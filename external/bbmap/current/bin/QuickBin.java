package bin;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.concurrent.locks.ReadWriteLock;
import java.util.concurrent.locks.ReentrantReadWriteLock;

import fileIO.ByteStreamWriter;
import fileIO.ReadWrite;
import shared.Parse;
import shared.Parser;
import shared.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import sketch.DisplayParams;
import sketch.Sketch;
import stream.ConcurrentReadInputStream;
import stream.ConcurrentReadOutputStream;
import stream.Read;
import structures.ByteBuilder;
import structures.IntHashMap;
import structures.IntLongHashMap;
import structures.ListNum;
import structures.LongList;
import tax.TaxTree;
import template.Accumulator;
import template.ThreadWaiter;
import tracker.ReadStats;

/**
 * Prototype for metagenome contig binning.
 * 
 * @author Brian Bushnell
 * @date December 6, 2024
 *
 */
public class QuickBin extends BinObject implements Accumulator<QuickBin.ProcessThread> {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Code entrance from the command line.
	 * @param args Command line arguments
	 */
	public static void main(String[] args){
		//Start a timer immediately upon code entrance.
		Timer t=new Timer();
		
		//Create an instance of this class
		QuickBin x=new QuickBin(args);
		
		//Run the object
		x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(BinObject.outstream);
	}
	
	/**
	 * Constructor.
	 * @param args Command line arguments
	 */
	public QuickBin(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		//Set shared static variables prior to parsing
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.setZipThreads(Shared.threads());
		
		{//Parse the arguments
			loader=new DataLoader(outstream);
			binner=new Binner(outstream);
			final Parser parser=parse(args);
			Parser.processQuality();
			
			overwrite=ReadStats.overwrite=parser.overwrite;
			append=ReadStats.append=parser.append;

			outPattern=parser.out1;
			extout=parser.extout;
		}

		validateParams();
		checkFileExistence(); //Ensure files can be read and written
		loader.checkInput();
		
		BinObject.tree=(sketchContigs || sketchClusters) ? loader.loadTree() : null;
		sketcher=(sketchContigs || sketchClusters || sketchOutput) ? new BinSketcher(16, 2000) : null;
		binner.sketcher=loader.sketcher=sketcher;
		Sketch.defaultParams.format=DisplayParams.FORMAT_JSON;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------    Initialization Helpers    ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Parse arguments from the command line */
	private Parser parse(String[] args){
		originalArgs=args.clone();
		
		//Create a parser object
		Parser parser=new Parser();
		
		//Set any necessary Parser defaults here
		//parser.foo=bar;
		
		//Parse each argument
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			
			//Break arguments into their constituent parts, in the form of "a=b"
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			if(b!=null && b.equalsIgnoreCase("null")){b=null;}
			
			if(a.equals("verbose")){
				verbose=Parse.parseBoolean(b);
			}
			
			else if(a.equals("clusterbytax") || a.equals("clusterbytaxid")){
				clusterByTaxid=Parse.parseBoolean(b);
			}else if(a.equals("clusterbytet") || a.equals("clusterbytetramer")){
				clusterByTetramer=Parse.parseBoolean(b);
			}else if(a.equals("refine") || a.equals("refineclusters")){
				refineClusters=Parse.parseBoolean(b);
			}else if(a.equals("residue") || a.equals("processresidue")){
				processResidue=Parse.parseBoolean(b);
			}else if(a.equals("entropy")){
				SpectraCounter.calcEntropy=Parse.parseBoolean(b);
			}
			
			else if(a.equalsIgnoreCase("sketchcontigs")){
				sketchContigs=Parse.parseBoolean(b);
			}else if(a.equalsIgnoreCase("sketchclusters") || a.equalsIgnoreCase("sketchbins")){
				sketchClusters=Parse.parseBoolean(b);
			}else if(a.equalsIgnoreCase("sketchoutput")){
				sketchOutput=Parse.parseBoolean(b);
			}else if(a.equalsIgnoreCase("sketch")){
				sketchClusters=sketchContigs=sketchOutput=Parse.parseBoolean(b);
			}else if(a.equalsIgnoreCase("density")){
				float f=Float.parseFloat(b);
				if(f<1) {sketchDensity=f;}
				else {sketchDensity=1/f;}
				
			}else if(a.equals("sketchbulk") || a.equalsIgnoreCase("sketchinbulk") || a.equals("bulk")){
				sketchInBulk=Parse.parseBoolean(b);
			}else if(a.equals("sketchsectionsize")){
				BinSketcher.sectionSize=Tools.mid(1, Integer.parseInt(b), 10000);
			}else if(a.equals("quant") || a.equals("quantize") || a.equals("quantizer")){
				BinObject.setQuant(Tools.max(1, Integer.parseInt(b)));
			}
			
			else if(a.equals("mincluster") || a.equals("minclustersize")){
				minBasesPerCluster=Parse.parseIntKMG(b);
			}else if(a.equals("mincontigs")){
				minContigsPerCluster=Parse.parseIntKMG(b);
			}else if(a.equals("validate") || a.equals("validation")){
				validation=Parse.parseBoolean(b);
			}else if(a.equalsIgnoreCase("followEdges1") || a.equals("e1")){
				followEdge1Passes=(Tools.startsWithDigit(b) ? Integer.parseInt(b) :
					Parse.parseBoolean(b) ? 4 : 0);
			}else if(a.equalsIgnoreCase("followEdges2") || a.equals("e2") || a.equals("followEdges")){
				followEdge2Passes=(Tools.startsWithDigit(b) ? Integer.parseInt(b) :
					Parse.parseBoolean(b) ? 4 : 0);
			}else if(a.equalsIgnoreCase("edgestringency1")){
				edgeStringency1=Float.parseFloat(b);
			}else if(a.equalsIgnoreCase("edgestringency2") || a.equalsIgnoreCase("edgestringency")){
				edgeStringency2=Float.parseFloat(b);
			}
			
			else if(a.equals("covout") || a.equals("outcov")) {
				covOut=b;
			}
			
			else if(a.equals("parse_flag_goes_here")){
				long fake_variable=Parse.parseKMG(b);
				//Set a variable here
			}else if(binner.parse(arg, a, b)){
				//do nothing
			}else if(SimilarityMeasures.parse(arg, a, b)){
				//do nothing
			}else if(loader.parse(arg, a, b)){
				//do nothing
			}else if(parser.parse(arg, a, b)){//Parse standard flags in the parser
				//do nothing
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
			}
		}
		
		return parser;
	}
	
	private void reprocessArgs() {
		for(int i=0; i<originalArgs.length; i++){
			String arg=originalArgs[i];
			
			//Break arguments into their constituent parts, in the form of "a=b"
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			if(b!=null){binner.parse(arg, a, b);}
		}
	}
	
	/** Ensure parameter ranges are within bounds and required parameters are set */
	private boolean validateParams(){
//		assert(minfoo>0 && minfoo<=maxfoo) : minfoo+", "+maxfoo;
//		assert(false) : "TODO";
		return true;
	}
	
	/** Ensure files can be read and written */
	private void checkFileExistence(){
		//Ensure output files can be written
		String o=(outPattern!=null && outPattern.indexOf('%')<0) ? outPattern : null;
		if(!Tools.testOutputFiles(overwrite, append, false, o)){
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output file "+o+"\n");
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	void process(Timer t) {
		process_mapVersion(t);
	}

	/** Create read streams and process all data */
	void process_mapVersion(Timer t){
		
		Timer t2=new Timer();
		
		contigList=loader.loadData();
		assert(isValid(contigList, false));
		if(covOut!=null) {DataLoader.writeCov(covOut, contigList, loader.numDepths, outstream);}
//		System.err.println("Thresholds:");
//		binner.printThresholds();
		binner.setSamples(loader.numDepths);
//		System.err.println("Thresholds:");
//		binner.printThresholds();
		reprocessArgs();
//		System.err.println("Thresholds:");
//		binner.printThresholds();
//		assert(false);
		
		IntLongHashMap sizeMap=loader.sizeMap;
//		@SuppressWarnings("unchecked")
//		ArrayList<? extends Bin> bins=(ArrayList<Contig>) contigList.clone();
//		if(clusterByTaxid) {
//			bins=binner.clusterByTaxid(bins);
//		}
		
		ArrayList<Bin> binList=null;
		if(followEdge1Passes>0 && loader.makePairGraph) {
			t2.start();
			outstream.println("Following Edges:");
			assert(isValid(contigList, false));
			int sum=0;
			int merges=binner.followEdges(contigList, edgeStringency1);
			sum+=merges;
			assert(isValid(contigList, true));
			for(int i=1; i<followEdge1Passes && merges>0; i++) {
				binList=Binner.toBinList(contigList);
				assert(isValid(binList, false));
				merges=binner.followEdges(contigList, binList, edgeStringency1);
				sum+=merges;
			}
			if(merges>0) {
				binList=Binner.toBinList(contigList);
			}
			t2.stop("Merged "+sum+" bins:   ");
		}
		
		BinMap binMap=binner.makeBinMap(contigList, binList);
		assert(binMap.isValid());
		binList=null;
//		comparisonsCreate=binMap.comparisons;
//		slowComparisonsCreate=binMap.slowComparisons;
		
		if(refineClusters) {
			binner.refineBinMap(binMap);
			assert(binMap.isValid());
		}
		
		if(followEdge2Passes>0 && loader.makePairGraph) {
			t2.start();
			outstream.println("Following Edges:");
			int merges=100, sum=0;
			for(int i=0; i<followEdge2Passes && merges>1; i++) {
				assert(isValid(contigList, true));
				binList=Binner.toBinList(contigList);
				assert(isValid(binList, false));
				merges=binner.followEdges(contigList, binList, edgeStringency2);
				sum+=merges;
			}
			t2.stop("Merged "+sum+" bins:   ");
			if(sum>0) {
				binMap.clear(true);
				binMap.addAll(Binner.toBinList(contigList), Binner.minSizeToMerge);
			}
		}
		
		if(processResidue) { 
			int removedResidue=binner.processResidue(binMap, Binner.residueStringency, 
					TaxTree.SPECIES, true, true, binner.residueRange);
			if(removedResidue>1) {
				//Honestly it would be best to have a bigger compare size limit for residue prior to this
				int removedThisPhase=binner.refinePhase(binMap, 
						"f", 1.0f, -1, true, true, binner.baseRange+2, 2);
			}
		}
//		comparisonsRefine=binMap.comparisons-comparisonsCreate;
//		slowComparisonsRefine=binMap.slowComparisons-slowComparisonsCreate;
		
		t2.stop();
		ArrayList<Cluster> bins=binMap.toList(true);
		Collections.sort(bins);
//		assert(isValid(bins, false));
		
		if(sketchOutput) {
			System.err.println("Sketching output.");
			sketcher.sketch(bins, true);
		}
		
		long cleanBins=0, cleanSize=0, partialCleanSize=0;
		long contamBins=0, contamSize=0, partialContamSize=0;
		double contamScore=0;
		double compltScore=0;
		long sizeOverLimit=0, contigsOverLimit=0, binsOverLimit=0;
		IntHashMap tidBins=new IntHashMap();
		
//		t2.start();
//		System.err.print("Calculating statistics.");
		for(Bin b : bins) {
			if(b.size()<minBasesPerCluster) {continue;}
			if(validation) {b.calcContam(sizeMap);}
			long contam=(int)(b.size()*b.contam);
			compltScore+=((b.size()-contam)*b.completeness);
			contamScore+=(b.size()*b.contam);
			if(contam>0) {
				contamBins++;
				contamSize+=b.size();
				partialContamSize+=contam;
				partialCleanSize+=(b.size()-contam);
			}else {
				cleanBins++;
				cleanSize+=b.size();
			}
			if(b.size()>=minBasesPerCluster) {
				sizeOverLimit+=b.size();
				binsOverLimit++;
				contigsOverLimit+=b.numContigs();
				if(b.labelTaxid>0) {tidBins.increment(b.labelTaxid);}
			}
		}
//		t2.stopAndPrint();
		
//		System.err.println("\nFinal clusters:");
//		for(int i=0; i<bins.size(); i++) {
//			Bin a=bins.get(i);
//			if(a.size()>=4000 || a.numContigs()>1){
//				System.out.println(a.toBytes());
//			}
//			if(a.size()<1000 || i>=4) {break;}
//		}
		
		outputClusters(outPattern, bins, minBasesPerCluster, minContigsPerCluster);
//		
//		//Create a read input stream
//		final ConcurrentReadInputStream cris=makeCris(ffin1, ffin2);
//		
//		//Optionally create a read output stream
//		final ConcurrentReadOutputStream ros=makeCros(outPattern);
//		
//		//Reset counters
//		readsProcessed=readsOut=0;
//		basesProcessed=basesOut=0;
//		
//		//Process the reads in separate threads
//		spawnThreads(cris, ros);
		
		if(verbose){outstream.println("Finished; closing streams.");}
		
//		//Write anything that was accumulated by ReadStats
//		errorState|=ReadStats.writeAll();
//		//Close the read streams
//		errorState|=ReadWrite.closeStreams(cris, ros);
		
		//Report timing and results
		outstream.println();
		t.stopAndPrint();
//		outstream.println(Tools.timeReadsBasesProcessed(t, readsProcessed, basesProcessed, 8));
//		outstream.println(Tools.readsBasesOut(readsProcessed, basesProcessed, readsOut, basesOut, 8, false));

//		outstream.println("Initial Fast Comparisons:    \t"+comparisonsCreate);
//		outstream.println("Initial Slow Comparisons:    \t"+slowComparisonsCreate);
//		outstream.println("Refinement Fast Comparisons: \t"+comparisonsRefine);
//		outstream.println("Refinement Slow Comparisons: \t"+slowComparisonsRefine);
//		float cps=((binMap.slowComparisons+binner.refinementComparisonsSlow)/(float)t2.elapsed)*1000000000;
//		outstream.println("Comparisons Per Second:      \t"+Tools.padKMB((long)cps, 0));
		outstream.println();

		if(validation) {
			outstream.println(formatString("Good Merges Follow", 29, binner.goodMergesFollow, binner.badMergesFollow));
			outstream.println(formatString("Bad Merges Follow", 29, binner.badMergesFollow, binner.goodMergesFollow));
			outstream.println(formatString("Good Merges Create", 29, binner.goodMergesCreate, binner.badMergesCreate));
			outstream.println(formatString("Bad Merges Create", 29, binner.badMergesCreate, binner.goodMergesCreate));
			outstream.println(formatString("Good Merges Refine", 29, binner.goodMergesRefine, binner.badMergesRefine));
			outstream.println(formatString("Bad Merges Refine", 29, binner.badMergesRefine, binner.goodMergesRefine));
			outstream.println(formatString("Good Merges Residue", 29, binner.goodMergesResidue, binner.badMergesResidue));
			outstream.println(formatString("Bad Merges Residue", 29, binner.badMergesResidue, binner.goodMergesResidue));
			outstream.println();
			outstream.println(formatString("Clean Bins", 29, cleanBins, contamBins));
			outstream.println(formatString("Dirty Bins", 29, contamBins, cleanBins));
			outstream.println(formatString("Clean Bin Bases", 29, cleanSize, contamSize));
			outstream.println(formatString("Dirty Bin Bases", 29, contamSize, cleanSize));
			outstream.println(formatString("Tainted Bases", 29, partialCleanSize, cleanSize+contamSize-partialCleanSize));
			outstream.println(formatString("Contam Bases", 29, partialContamSize, cleanSize+contamSize-partialContamSize));
			outstream.println();
		}

		outstream.println("Sequence Recovery:           \t"+
				String.format("%.3f", sizeOverLimit*100.0/loader.basesLoaded));
		outstream.println("Contig Recovery:             \t"+
				String.format("%.3f", contigsOverLimit*100.0/loader.contigsLoaded));

		if(validation) {
			outstream.println("Genomes Represented:         \t"+
					String.format("%.3f", (tidBins.size())*100.0/sizeMap.size()));
			outstream.println("Completeness Score:          \t"+
					String.format("%.3f", 100*compltScore/loader.basesLoaded));
			outstream.println("Contamination Score:         \t"+
					String.format("%.4f", 100*contamScore/loader.basesLoaded));
		}
		outstream.println();
		
		printL90(bins, loader.basesLoaded);
		
		//Throw an exception of there was an error in a thread
		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}
	
	static String formatString(String term, int len, long a, long b) {
		float pct=a*100f/(a+b);
		while(term.length()<len) {term=term+" ";}
		String apad=Long.toString(a);
		while(apad.length()<12) {apad+=" ";}
		return(term+"\t"+apad+"\t"+String.format("%.3f%%", pct));
	}
	
	private void printL90(ArrayList<? extends Bin> bins, long basesLoaded) {
		LongList sizes=new LongList(bins.size());
		for(Bin b : bins) {
			sizes.add(b.size());
		}
		GradeBins.printL90(sizes, basesLoaded);
	}
	
	private void outputClusters(String pattern, ArrayList<? extends Bin> clusters, long minBases, int minContigs) {
//		if(pattern==null) {return;}
		if(pattern!=null) {outstream.println("Writing clusters to "+pattern);}
		long sizeOverLimit=0;
		if(pattern!=null && pattern.indexOf('%')>=0) {
			final ByteStreamWriter chaff=ByteStreamWriter.makeBSW(pattern.replace("%", "chaff"), overwrite, append, true);
			final ByteBuilder bb=new ByteBuilder(8192);
			for(int i=0; i<clusters.size(); i++) {
				Bin a=clusters.get(i);
				if(a.size()>=minBases && a.numContigs()>=minContigs) {
					String fname=pattern.replace("%", Integer.toString(i));
					final ByteStreamWriter bsw=ByteStreamWriter.makeBSW(fname, overwrite, append, true);
					printBin(a, bsw, bb, -1);
					bsw.poison();
					clustersWritten++;
					contigsWritten+=a.numContigs();
					basesWritten+=a.size();
				}else {
					printBin(a, chaff, bb, i+1);
				}
			}
			chaff.poisonAndWait();
		}else {
			final ByteBuilder bb=new ByteBuilder(8192);
			final ByteStreamWriter bsw=ByteStreamWriter.makeBSW(pattern, overwrite, append, true);
			for(int i=0; i<clusters.size(); i++) {
				Bin a=clusters.get(i);
				printBin(a, bsw, bb, i+1);
				clustersWritten++;
				contigsWritten+=a.numContigs();
				basesWritten+=a.size();
			}
			if(bsw!=null) {bsw.poisonAndWait();}
		}
		float cpct=contigsWritten*100f/loader.contigsLoaded;
		float bpct=basesWritten*100f/loader.basesLoaded;
		outstream.println("\nMetric   \t        In\t       Out\tPercent");
		outstream.println("Clusters \t"+Tools.padLeft(0, 10)+"\t"+Tools.padLeft(clustersWritten, 10));
		outstream.println("Contigs  \t"+Tools.padLeft(loader.contigsLoaded, 10)+"\t"+
				Tools.padLeft(contigsWritten, 10)+"\t"+Tools.format("%.2f%%", cpct));
		outstream.println("Bases    \t"+Tools.padLeft(loader.basesLoaded, 10)+"\t"+
				Tools.padLeft(basesWritten, 10)+"\t"+Tools.format("%.2f%%", bpct));
	}
	
	private void printBin(Bin a, ByteStreamWriter bsw, ByteBuilder bb, int id) {
		if(a.isCluster()) {printCluster((Cluster)a, bsw, bb, id);}
		else {printContig((Contig)a, bsw, bb, id);}
	}
	
	private void printCluster(Cluster a, ByteStreamWriter bsw, ByteBuilder bb, int id) {
		ArrayList<Contig> contigs=a.contigs;
		Collections.sort(contigs);
		for(Contig c : contigs) {
			c.appendTo(bb, id);
			if(bb.length>4096) {
				if(bsw!=null) {bsw.print(bb);}
				bb.clear();
			}
		}
		if(bsw!=null && !bb.isEmpty()) {bsw.print(bb);}
		bb.clear();
	}
	
	private void printContig(Contig c, ByteStreamWriter bsw, ByteBuilder bb, int id) {
		c.appendTo(bb, id);
		if(bsw!=null && !bb.isEmpty()) {bsw.print(bb);}
		bb.clear();
	}
	
//	private ConcurrentReadInputStream makeCris(FileFormat ff1, FileFormat ff2){
//		ConcurrentReadInputStream cris=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ff1, ff2);
//		cris.start(); //Start the stream
//		if(verbose){outstream.println("Started cris");}
//		boolean paired=cris.paired();
////		if(!ffin1.samOrBam()){outstream.println("Input is being processed as "+(paired ? "paired" : "unpaired"));}
//		return cris;
//	}
//	
//	private ConcurrentReadOutputStream makeCros(String fname){
//		if(fname==null) {return null;}
//		FileFormat ff=FileFormat.testOutput(fname, FileFormat.FASTA, null, true, overwrite, append, false);
//
//		//Select output buffer size based on whether it needs to be ordered
//		final int buff=8;
//
//		final ConcurrentReadOutputStream ros=ConcurrentReadOutputStream.getStream(ff, null, buff, null, false);
//		ros.start(); //Start the stream
//		return ros;
//	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	
	
	/*--------------------------------------------------------------*/
	/*----------------       Thread Management      ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Spawn process threads */
	private void spawnThreads(final ConcurrentReadInputStream cris, final ConcurrentReadOutputStream ros){
		
		//Do anything necessary prior to processing
		
		//Determine how many threads may be used
		final int threads=Shared.threads();
		
		//Fill a list with ProcessThreads
		ArrayList<ProcessThread> alpt=new ArrayList<ProcessThread>(threads);
		for(int i=0; i<threads; i++){
			alpt.add(new ProcessThread(cris, ros, i));
		}
		
		//Start the threads and wait for them to finish
		boolean success=ThreadWaiter.startAndWait(alpt, this);
		errorState&=!success;
		
		//Do anything necessary after processing
		
	}
	
	@Override
	public final void accumulate(ProcessThread pt){
		synchronized(pt) {
//			readsProcessed+=pt.readsProcessedT;
//			basesProcessed+=pt.basesProcessedT;
//			readsOut+=pt.readsOutT;
//			basesOut+=pt.basesOutT;
			errorState|=(!pt.success);
		}
	}
	
	@Override
	public final boolean success(){return !errorState;}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** This class is static to prevent accidental writing to shared variables.
	 * It is safe to remove the static modifier. */
	static class ProcessThread extends Thread {
		
		//Constructor
		ProcessThread(final ConcurrentReadInputStream cris_, final ConcurrentReadOutputStream ros_, final int tid_){
			cris=cris_;
			ros=ros_;
			tid=tid_;
		}
		
		//Called by start()
		@Override
		public void run(){
			//Do anything necessary prior to processing
			
			//Process the reads
			processInner();
			
			//Do anything necessary after processing
			
			//Indicate successful exit status
			success=true;
		}
		
		/** Iterate through the reads */
		void processInner(){
			
			//Grab the first ListNum of reads
			ListNum<Read> ln=cris.nextList();

			//Check to ensure pairing is as expected
			if(ln!=null && !ln.isEmpty()){
				Read r=ln.get(0);
//				assert(ffin1.samOrBam() || (r.mate!=null)==cris.paired()); //Disabled due to non-static access
			}

			//As long as there is a nonempty read list...
			while(ln!=null && ln.size()>0){
//				if(verbose){outstream.println("Fetched "+reads.size()+" reads.");} //Disabled due to non-static access
				
				processList(ln);
				
				//Notify the input stream that the list was used
				cris.returnList(ln);
//				if(verbose){outstream.println("Returned a list.");} //Disabled due to non-static access
				
				//Fetch a new list
				ln=cris.nextList();
			}

			//Notify the input stream that the final list was used
			if(ln!=null){
				cris.returnList(ln.id, ln.list==null || ln.list.isEmpty());
			}
		}
		
		void processList(ListNum<Read> ln){

			//Grab the actual read list from the ListNum
			final ArrayList<Read> reads=ln.list;
			
			//Loop through each read in the list
			for(int idx=0; idx<reads.size(); idx++){
				final Read r1=reads.get(idx);
				final Read r2=r1.mate;
				
				//Validate reads in worker threads
				if(!r1.validated()){r1.validate(true);}
				if(r2!=null && !r2.validated()){r2.validate(true);}

				//Track the initial length for statistics
				final int initialLength1=r1.length();
				final int initialLength2=r1.mateLength();

				//Increment counters
				readsProcessedT+=r1.pairCount();
				basesProcessedT+=initialLength1+initialLength2;
				
				{
					//Reads are processed in this block.
					boolean keep=processReadPair(r1, r2);
					
					if(!keep){reads.set(idx, null);}
					else{
						readsOutT+=r1.pairCount();
						basesOutT+=r1.pairLength();
					}
				}
			}

			//Output reads to the output stream
			if(ros!=null){ros.add(reads, ln.id);}
		}
		
		/**
		 * Process a read or a read pair.
		 * @param r1 Read 1
		 * @param r2 Read 2 (may be null)
		 * @return True if the reads should be kept, false if they should be discarded.
		 */
		boolean processReadPair(final Read r1, final Read r2){
			throw new RuntimeException("TODO: Implement this method."); //TODO
//			return true;
		}

		/** Number of reads processed by this thread */
		protected long readsProcessedT=0;
		/** Number of bases processed by this thread */
		protected long basesProcessedT=0;
		
		/** Number of reads retained by this thread */
		protected long readsOutT=0;
		/** Number of bases retained by this thread */
		protected long basesOutT=0;
		
		/** True only if this thread has completed successfully */
		boolean success=false;
		
		/** Shared input stream */
		private final ConcurrentReadInputStream cris;
		/** Shared output stream */
		private final ConcurrentReadOutputStream ros;
		/** Thread ID */
		final int tid;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Primary output file path; should contain % symbol for bins */
	private String outPattern=null;
	
	/** Override output file extension */
	private String extout=null;
	
	private String covOut=null;
	
	private ArrayList<Contig> contigList;
	
	boolean clusterByTaxid=false;
	boolean clusterByTetramer=true;
	boolean refineClusters=true;
	boolean processResidue=true;
	int followEdge1Passes=0;
	int followEdge2Passes=4;
	float edgeStringency1=0.25f;
	float edgeStringency2=2f;
	
	/*--------------------------------------------------------------*/

//	/** Number of reads processed */
//	protected long readsProcessed=0;
//	/** Number of bases processed */
//	protected long basesProcessed=0;
//
//	/** Number of reads retained */
//	protected long readsOut=0;
//	/** Number of bases retained */
//	protected long basesOut=0;

	private long clustersWritten=0;
	private long contigsWritten=0;
	private long basesWritten=0;
	
	private DataLoader loader;
	private Binner binner;
	

	private long comparisonsCreate=0, slowComparisonsCreate=0;
	private long comparisonsRefine=0, slowComparisonsRefine=0;
	
	/*--------------------------------------------------------------*/
	/*----------------         Final Fields         ----------------*/
	/*--------------------------------------------------------------*/

//	private final Timer phaseTimer=new Timer();
	
	private final BinSketcher sketcher;
	
	@Override
	public final ReadWriteLock rwlock() {return rwlock;}
	private final ReadWriteLock rwlock=new ReentrantReadWriteLock();
	private String[] originalArgs;
	
	/*--------------------------------------------------------------*/
	/*----------------        Common Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Print verbose messages */
	public static boolean verbose=false;
	/** True if an error was encountered */
	public boolean errorState=false;
	/** Overwrite existing output files */
	private boolean overwrite=true;
	/** Append to existing output files */
	private boolean append=false;
	
}
