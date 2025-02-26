package bin;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.Map.Entry;
import bloom.BloomFilter;
import bloom.KmerCountAbstract;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import shared.LineParser1;
import shared.LineParser2;
import shared.LineParserS1;
import shared.LineParserS4;
import shared.Parse;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.ConcurrentReadInputStream;
import stream.FastaReadInputStream;
import stream.Read;
import stream.SamLine;
import stream.SamLineStreamer;
import structures.ByteBuilder;
import structures.IntHashMap;
import structures.IntLongHashMap;
import structures.ListNum;
import tax.TaxTree;
import ukmer.Kmer;

public class DataLoader extends BinObject {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	DataLoader(PrintStream outstream_){outstream=outstream_;}
	
	boolean parse(String arg, String a, String b) {
	
		if(a.equals("tree")){
			if(b==null || b.equalsIgnoreCase("t") || b.equalsIgnoreCase("true")) {
				treePath="auto";
			}else if(b.equalsIgnoreCase("f") || b.equalsIgnoreCase("false")) {
				treePath=null;
			}else {
				treePath=b;
			}
		}else if(a.equals("maxdepths") || a.equals("maxdepthcount") || a.equals("samples")){
			MAX_DEPTH_COUNT=Parse.parseIntKMG(b);
		}else if(a.equals("maxreads")){
			maxReads=Parse.parseKMG(b);
		}else if(a.equals("mincontig") || a.equals("minsizetoload")){
			minContigToLoad=Parse.parseIntKMG(b);
		}else if(a.equals("mapq") || a.equals("minmapq")){
			minMapq=Integer.parseInt(b);
		}
		
		else if(a.equals("bloomk") || a.equals("kbloom")){
			bloomkbig=Integer.parseInt(b);
			bloomkbig=Kmer.getKbig(bloomkbig);
		}else if(a.equals("bits") || a.equals("cbits")){
			cbits=Integer.parseInt(b);
		}else if(a.equals("ignoremissingcontigs")){
			ignoreMissingContigs=Parse.parseBoolean(b);
		}else if(a.equals("hashes")){
			hashes=Integer.parseInt(b);
		}
		
		else if(a.equals("makegraph") || a.equals("makepairgraph")){
			makePairGraph=Parse.parseBoolean(b);
		}
		
		else if(a.equals("ref") || a.equals("contigs") || a.equals("scaffolds") || a.equals("in")){
			ref=b;
		}else if(a.equals("readsin") || a.equals("reads1") || a.equals("read1") || a.equals("reads") || a.equals("sam") || a.equals("bam")){
			readFiles.clear();
			Tools.getFileOrFiles(b, readFiles, false, true, true, false);
		}else if(a.equals("covstats")){
			covstats.clear();
			Tools.getFileOrFiles(b, covstats, true, true, true, true);
		}else if(a.equals("cov") || a.equals("covin")){
			covIn=b;
		}else if(new File(arg).isFile()) {
			FileFormat ff=FileFormat.testInput(arg, FileFormat.SAM, null, false, false);
			if(ff.samOrBam() || ff.fastq()) {
				System.err.println("Adding reads file "+ff.simpleName());
				readFiles.add(arg);
			}else if(ff.fasta()) {
				System.err.println("Adding ref file "+ff.simpleName());
				ref=arg;
			}else {
				System.err.println("Adding covstats file "+ff.simpleName());
				covstats.add(arg);
			}
		}else if(a.equals("samserial") || a.equals("loadsamserial")){
			loadSamSerial=Parse.parseBoolean(b);
		}else if(a.equals("samparallel") || a.equals("loadsamparallel")){
			loadSamSerial=!Parse.parseBoolean(b);
		}
		
		else {return false;}
		
		return true;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------          Validation          ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Adjust file-related static fields as needed for this program */
	private static void checkStatics(){
		//Adjust the number of threads for input file reading
		if(!ByteFile.FORCE_MODE_BF1 && !ByteFile.FORCE_MODE_BF2 && Shared.threads()>2){
			ByteFile.FORCE_MODE_BF2=true;
		}
		
		assert(FastaReadInputStream.settingsOK());
	}
	
	void checkInput() {
		//TODO; validate input files
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Wrapper           ----------------*/
	/*--------------------------------------------------------------*/
	
	ArrayList<Contig> loadData() {
		//Turn off read validation in the input threads to increase speed
		final boolean vic=Read.VALIDATE_IN_CONSTRUCTOR;
		Read.VALIDATE_IN_CONSTRUCTOR=Shared.threads()<4;
		
		final boolean parseTax=true;
		final boolean parseCov=(covIn==null && readFiles.isEmpty() && covstats.isEmpty());

		final ArrayList<Contig> contigs=loadContigs(ref, minContigToLoad);
		
		if(parseCov) {
			depthCalculated=true;
		}else {
			calculateDepth(contigs);
		}
		assert(depthCalculated);
		
		
		outstream.print("Calculating kmer frequency spectra: \t");
		phaseTimer.start();
		
		SpectraCounter sc=new SpectraCounter(outstream, parseTax, parseCov);
		sc.makeSpectra(contigs);
		sizeMap=sc.sizeMap;
		phaseTimer.stopAndPrint();
		
		if(validation && makePairGraph) {
			validateEdges(contigs);
		}
		
		if(sketchContigs) {
			System.err.println("Sketching contigs. ");
			sketcher.sketch(contigs, true);
		}
		
		//Reset read validation
		Read.VALIDATE_IN_CONSTRUCTOR=vic;
		numDepths=(contigs.isEmpty()) ? 0 : contigs.get(0).numDepths();
		assert(numDepths>0) : "At least one depth is required.";
		return contigs;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Data Loading         ----------------*/
	/*--------------------------------------------------------------*/
	
	void calculateDepth(ArrayList<Contig> contigs) {
		final boolean sam=(!readFiles.isEmpty() && FileFormat.isSamOrBamFile(readFiles.get(0)));

		if(covIn!=null) {
			loadCovFile(covIn, contigs, MAX_DEPTH_COUNT);
		}else if(covstats.size()>0) {
			HashMap<String, Contig> contigMap=toMap(contigs);
			for(int i=0; i<readFiles.size(); i++) {
				calcDepthFromCovStats(covstats.get(i), i, contigMap, minContigToLoad);
			}
		}else if(sam) {
			Timer t=new Timer(outstream);
			HashMap<String, Contig> contigMap=toMap(contigs);
			if(loadSamSerial) {
				for(int i=0; i<readFiles.size(); i++) {
					FileFormat ff=FileFormat.testInput(readFiles.get(i), FileFormat.SAM, null, true, false);
					calcDepthFromSam(ff, i, contigMap);
				}
			}else {
				if(makePairGraph) {graph=new IntHashMap[contigMap.size()];}
				SamLoader samLoader=new SamLoader(outstream);
				samLoader.minMapq=minMapq;
				samLoader.minMateq=minMapq;
				samLoader.load(readFiles, contigMap, contigs, graph);
				if(makePairGraph) {
					synchronized(graph) {
						for(int i=0; i<contigs.size(); i++) {
							Contig c=contigs.get(i);
							IntHashMap pairMap=graph[i];
							if(pairMap!=null) {
								synchronized(c) {
									synchronized(pairMap) {c.pairMap=pairMap;}
								}
							}
						}
					}
				}
				depthCalculated=true;
			}
			t.stop("Loaded "+Tools.plural(" read file", readFiles.size())+": ");
		}
		
		if(readFiles.isEmpty() && !depthCalculated) {
			calcDepthFromHeader(contigs);
		}else if(!readFiles.isEmpty() && !sam) {
			for(int i=0; i<readFiles.size(); i++) {
				BloomFilter bloomFilter=makeBloomFilter(readFiles.get(i), null);
				calcDepthFromBloomFilter(contigs, bloomFilter, 0);
			}
		}
	}
	
	void validateEdges(ArrayList<Contig> contigs) {
		//TODO: add cov file only flag to quickbin
		//TODO: have BBMap suppress ambig or low-Q alignments
		//TODO: track good/bad rates of first, second, 3rd, and nth edges
		//TODO: suppress edges not near contig ends
		//TODO: Track good/bad rate as a function of relative weight to top edge.
		//TODO: Track good/bad rate as a function of edge weight.
		outstream.println("Validating edges: ");
		phaseTimer.start();
		long goodEdges=0, badEdges=0, otherEdges=0;
		long goodWeight=0, badWeight=0, otherWeight=0;
		long reciprocalGood=0, reciprocalBad=0, reciprocalOther=0;
		long sizeSumGood=0;
		long sizeSumBad=0;
		final int max=4;
		long[] edgeCount=new long[max+1];
		long[] weightCount=new long[max+1];
		for(Contig c : contigs) {
			IntHashMap pairMap=c.pairMap;
			if(pairMap==null) {
				edgeCount[0]++;
			}else {
				final int[] keys=pairMap.keys(), values=pairMap.values();
				final int invalid=pairMap.invalid();
				edgeCount[Tools.min(max, pairMap.size())]++;
				for(int kpos=0; kpos<keys.length; kpos++) {
					int dest=keys[kpos], weight=values[kpos];
					if(dest!=invalid) {
						weightCount[Tools.min(max, weight)]++;
						Contig c2=contigs.get(dest);
						boolean reciprocal=(c2.pairMap!=null && c2.pairMap.contains(c.id()));
						if(c.labelTaxid<1 || c2.labelTaxid<1) {
							otherEdges++;
							otherWeight+=weight;
							if(reciprocal) {reciprocalOther++;}
//							assert(false) : "\n"+c+"\n"+c2+"\n";
						}else if(c.labelTaxid==c2.labelTaxid) {
							goodEdges++;
							goodWeight+=weight;
							if(reciprocal) {reciprocalGood++;}
							sizeSumGood+=Tools.min(c.size(), c2.size());
						}else {
							badEdges++;
							badWeight+=weight;
							if(reciprocal) {reciprocalBad++;}
							sizeSumBad+=Tools.min(c.size(), c2.size());
						}
					}
				}
			}
		}
		phaseTimer.stopAndPrint();
		outstream.println("Type:       \tGood\tBad\tOther");
		outstream.println("Edges:      \t"+goodEdges+"\t"+badEdges+"\t"+otherEdges);
		outstream.println("Reciprocal: \t"+reciprocalGood+"\t"+reciprocalBad+"\t"+reciprocalOther);
		float w1=goodWeight*1f/goodEdges;
		float w2=badWeight*1f/badEdges;
		float w3=otherWeight*1f/otherEdges;
		outstream.println(String.format("Weights:    \t%.2f\t%.2f\t%.2f", w1, w2, w3));
		float s1=sizeSumGood*1f/goodEdges;
		float s2=sizeSumBad*1f/badEdges;
		outstream.println(String.format("Sizes:       \t%.2f\t%.2f", s1, s2));
		outstream.println("Counts:     \t0\t1\t2\t3\t4+");
		outstream.println("Edges:      \t"+Arrays.toString(edgeCount));
		outstream.println("Weights:    \t"+Arrays.toString(weightCount));
	}
	
	TaxTree loadTree() {
		if("auto".equals(treePath)){treePath=TaxTree.defaultTreeFile();}
		if(treePath!=null) {
			tree=TaxTree.loadTaxTree(treePath, outstream, false, false);
		}
		return tree;
	}
	
	HashMap<String, Contig> toMap(Collection<Contig> list){
		HashMap<String, Contig> map=new HashMap<String, Contig>(list.size());
		for(Contig c : list) {map.put(c.shortName, c);}
		return map;
	}
	
	ArrayList<Contig> loadContigs(String fname, int minlen){
		assert(fname!=null) : "No contigs specified.";
		outstream.print("Loading contigs:\t");
		phaseTimer.start();
		FileFormat ff=FileFormat.testInput(fname, FileFormat.FASTA, null, true, true);
		ArrayList<Contig> list=new ArrayList<Contig>();

		ConcurrentReadInputStream cris=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ff, null);
		cris.start(); //Start the stream
		if(verbose){outstream.println("Started cris");}
		
		//Grab the first ListNum of reads
		ListNum<Read> ln=cris.nextList();

		//Check to ensure pairing is as expected
		if(ln!=null && !ln.isEmpty()){
			Read r=ln.get(0);
			assert(r.mate==null);
		}

		//As long as there is a nonempty read list...
		while(ln!=null && ln.size()>0){
			
			for(Read r : ln) {
				Contig c=loadContig(r, minlen);
				if(c!=null) {list.add(c);}
			}
			
			//Notify the input stream that the list was used
			cris.returnList(ln);
//			if(verbose){outstream.println("Returned a list.");} //Disabled due to non-static access
			
			//Fetch a new list
			ln=cris.nextList();
		}

		//Notify the input stream that the final list was used
		if(ln!=null){
			cris.returnList(ln.id, ln.list==null || ln.list.isEmpty());
		}
		phaseTimer.stopAndPrint();
		Collections.sort(list);
		for(int i=0; i<list.size(); i++) {list.get(i).setID(i);}
		return list;
	}
	
	Contig loadContig(Read r, int minlen) {
		contigsLoaded++;
		basesLoaded+=r.length();
//		final int tid=parseTaxID(r.name());
//		if(tid>0) {sizeMap.increment(tid, r.length());}
		if(r.length()<minlen) {return null;}
		contigsRetained++;
		basesRetained+=r.length();
		Contig c=new Contig(r.name(), r.bases, (int)contigsRetained);
//		if(tid>0) {
//			c.labelTaxid=tid;
//			if(!validation) {c.taxid=tid;}
//		}
		return c;
	}
	
	public static int parseTaxID(String line) {
		if(!parseTaxid) {return -1;}
		final String term="tid_";
		int pos=line.indexOf(term);
		if(pos<0) {return -1;}
		long id=0;
		for(int i=pos+4; i<line.length(); i++) {
			char c=line.charAt(i);
			if(c<'0' || c>'9') {break;}
			id=id*10+(c-'0');
		}
		assert(id>0 && id<Integer.MAX_VALUE) : id+"\n"+line+"\n";
		return (int)id;
	}
	
	BloomFilter makeBloomFilter(String in1, String in2) {
//		if(ffin1.samOrBam()) {return null;}
		outstream.print("Making Bloom filter: \t");
		KmerCountAbstract.CANONICAL=true;
		bloomkbig=Kmer.getKbig(bloomkbig);
		int bloomk=Kmer.getK(bloomkbig);
		BloomFilter filter=new BloomFilter(in1, in2, null, bloomk, bloomkbig, cbits, hashes, 1,
				true, false, false, 0.5f);
		phaseTimer.stopAndPrint();
		outstream.println(filter.filter.toShortString());
		return filter;
	}
	
	void calcDepthFromCovStats(String fname, int sample, HashMap<String, Contig> contigMap, int minlen) {
		ByteFile bf=ByteFile.makeByteFile(fname, true);
		outstream.print("Loading covstats file ("+fname+"): \t");
		byte[] line=null;
		LineParser2 lp=new LineParser2('\t');
		int found=0;
		for(line=bf.nextLine(); line!=null; line=bf.nextLine()){
			if(Tools.startsWith(line, '#')) {
				//header
				assert(Tools.startsWith(line, "#ID\tAvg_fold")) : new String(line);
			}else {
				lp.set(line);
				String name=lp.parseString();
				String shortName=ContigRenamer.toShortName(name);
				float depth=lp.parseFloat();
				Contig c=contigMap.get(shortName);
				assert(c!=null || minlen>0 || ignoreMissingContigs) : 
					"Can't find contig that is specified in covstats: "+shortName;
				if(c!=null) {
					if(c!=null) {
						c.setDepth(depth, sample);
						found++;
					}
				}
			}
		}
		assert(found<=contigMap.size()) : "Duplicate entries found in covstats file.";
		assert(found==contigMap.size() || minlen>0) : "Some contigs were not found in covstats file.";
		assert(found>0) : "No matching entries found in covstats file.";
		depthCalculated=true;
		phaseTimer.stopAndPrint();
	}
	
	@Deprecated
	void calcDepthFromSam(FileFormat ff, final int sample, HashMap<String, Contig> contigMap) {
		SamLineStreamer ss=null;
		outstream.print("Loading sam file: \t");
		phaseTimer.start();
		final int streamerThreads=Tools.min(4, Shared.threads());
		ss=new SamLineStreamer(ff, streamerThreads, false, maxReads);
		ss.start();
		processSam(ss, sample, contigMap);
		for(Entry<String, Contig> e : contigMap.entrySet()) {
			Contig c=e.getValue();
			c.setDepth(c.depthSum/(Tools.max(1f, c.bases.length)), sample);
			c.depthSum=0;
		}
		depthCalculated=true;
		phaseTimer.stopAndPrint();
	}

	@Deprecated
	private void processSam(SamLineStreamer ss, final int sample, HashMap<String, Contig> contigMap) {
		ListNum<SamLine> ln=ss.nextLines();
		ArrayList<SamLine> reads=(ln==null ? null : ln.list);

		while(ln!=null && reads!=null && reads.size()>0){

			for(int idx=0; idx<reads.size(); idx++){
				SamLine sl=reads.get(idx);
				if(sl.mapped()) {
					String rname=ContigRenamer.toShortName(sl.rname());
					Contig c=contigMap.get(rname);
					assert(c!=null || minContigToLoad>0 || ignoreMissingContigs) : 
						"Can't find contig for rname "+rname+"\nsize="+contigMap.size()
						+"\n[0]="+contigMap.keySet().toArray(new String[0])[0];
					if(c!=null) {
						c.depthSum+=sl.length();
					}
				}
			}
			ln=ss.nextLines();
			reads=(ln==null ? null : ln.list);
		}
	}
	
	public void calcDepthFromHeader(Collection<Contig> list) {
		outstream.print("Parsing depth from contig headers: \t");
		phaseTimer.start();
		LineParserS1 lps=new LineParserS1('_');
		LineParserS4 lpt=new LineParserS4(",,=,");
		for(Contig c : list) {
			boolean b=parseAndSetDepth(c, lps, lpt);
			assert(b) : "Could not parse depth from "+c.name;
//			c.depthepth=parseDepth(c.name, lps, lpt);
		}
		depthCalculated=true;
		phaseTimer.stopAndPrint();
	}
	
	public static boolean parseAndSetDepth(Contig c, LineParserS1 lps, LineParserS4 lpt) {
		String name=c.name;
		c.clearDepth();
		if(name.contains(" cov_")) {//Multiple coverage
			String[] array=Tools.whitespacePlus.split(name);
			for(String s : array) {
				lps.set(s);
				if(lps.terms()==2 && lps.termEquals("cov", 0)) {
					float depth=lps.parseFloat(1);
					c.appendDepth(depth);
					if(c.numDepths()>=MAX_DEPTH_COUNT) {break;}
				}
			}
		}else if(name.startsWith("NODE_") && name.contains("_cov_")) {//Spades
			lps.set(name);
			float depth=lps.parseFloat(5);
			c.setDepth(depth, 0);
		}else if(name.startsWith("contig_") && name.contains(",cov=")) {//Tadpole
			lpt.set(name);
			float depth=lpt.parseFloat(3);
			c.setDepth(depth, 0);
		}else {
			return false;
		}
		return true;
	}
	
	public void calcDepthFromBloomFilter(Collection<Contig> list, BloomFilter bf, int sample) {
		outstream.print("Calculating depth from Bloom filter: \t");
		for(Contig c : list) {c.setDepth(bf.averageCount(c.bases), sample);}
		depthCalculated=true;
		phaseTimer.stopAndPrint();
	}
	
	public void loadCovFile(String fname, ArrayList<Contig> contigs, final int maxSamples) {
		outstream.print("Loading coverage from "+fname+": \t");
		phaseTimer.start();
		LineParser1 lp=new LineParser1('\t');
		ByteFile bf=ByteFile.makeByteFile(fname, true);
		
		byte[] line=bf.nextLine();
		int numDepths=0;
		for(; Tools.startsWith(line, '#'); line=bf.nextLine()) {
			if(Tools.startsWith(line, "#Depths")) {
				numDepths=lp.set(line).parseInt(1);
			}
		}
		assert(numDepths>0) : numDepths;
		final int edgeStart=3+numDepths;
		final int samples=Tools.min(numDepths, maxSamples);
		int loaded=0;
		for(; line!=null; line=bf.nextLine()) {
			lp.set(line);
			int id=lp.parseInt(1);
			if(id>=contigs.size()) {continue;}//Presumably small contigs that were not loaded
			int size=lp.parseInt(2);
			int edges=(lp.terms()-edgeStart)/2;
			Contig c=contigs.get(id);
			assert(c.id()==id);
			assert(c.size()==size);
			assert(lp.termEquals(c.shortName, 0));
			assert(c.numDepths()==0);
			for(int i=0; i<samples; i++) {
				float f=lp.parseFloat(i+3);
				c.setDepth(f, i);
			}
			if(lp.terms()>edgeStart && makePairGraph) {
				c.pairMap=new IntHashMap((1+(edges*4)/3)|1);
				for(int term=edgeStart; term<lp.terms(); term+=2) {
					int dest=lp.parseInt(term);
					int weight=lp.parseInt(term+1);
					c.pairMap.put(dest, weight);
				}
			}
			loaded++;
		}
		assert(loaded==contigs.size()) : loaded+", "+contigs.size();
		depthCalculated=true;
		phaseTimer.stopAndPrint();
	}
	
	public static void writeCov(String fname, Collection<Contig> contigs, int numDepths, PrintStream outstream) {
		outstream.print("Writing coverage to "+fname+": ");
		Timer t=new Timer(outstream);
		ByteStreamWriter bsw=new ByteStreamWriter(fname, true, false, true);
		bsw.start();
		ByteBuilder bb=new ByteBuilder();
		bb.append("#Contigs\t").append(contigs.size()).nl();
		bb.append("#Depths\t").append(numDepths).nl();
		bb.append("#ShortName\tID\tSize");
		for(int i=0; i<numDepths; i++) {bb.tab().append("Cov_"+i);}
		bb.tab().append("Edge\tWeight");
		bsw.println(bb);
		for(Contig c : contigs) {
			c.toCov(bb.clear());
			bsw.print(bb.nl());
		}
		bsw.poisonAndWait();
		t.stopAndPrint();
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Assembly path */
	String ref=null;
	String covIn=null;
	
	private final ArrayList<String> covstats=new ArrayList<String>();
	
	/** Primary read input file path */
	private final ArrayList<String> readFiles=new ArrayList<String>();
	
//	/** Primary read input file */
//	FileFormat ffin1;
	
	boolean depthCalculated=false;
	boolean ignoreMissingContigs=false;
	
	long contigsLoaded=0;
	long basesLoaded=0;
	
	long contigsRetained=0;
	long basesRetained=0;
	
	private int bloomkbig=31;
	private int cbits=16;
	private int hashes=3;
	
	BinSketcher sketcher;

	int minContigToLoad=100;
	int minMapq=20;
	
	/** Quit after processing this many input reads; -1 means no limit */
	private long maxReads=-1;
	
	String treePath="auto";
	TaxTree tree;
	
	IntLongHashMap sizeMap;
	IntHashMap[] graph;
	
	static int MAX_DEPTH_COUNT=999;
	static boolean loadSamSerial=false;
	boolean makePairGraph=true;
	int numDepths=0;
	
	/*--------------------------------------------------------------*/
	/*----------------         Final Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	final Timer phaseTimer=new Timer();
	final PrintStream outstream;
//	static final Pattern covPattern=Pattern.compile(" cov_");
	
}
