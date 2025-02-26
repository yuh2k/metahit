package bin;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Random;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.atomic.AtomicLong;

import dna.AminoAcid;
import fileIO.ByteFile;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import shared.Parse;
import shared.Parser;
import shared.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.ConcurrentReadInputStream;
import stream.ConcurrentReadOutputStream;
import stream.FASTQ;
import stream.FastaReadInputStream;
import stream.Read;
import structures.ByteBuilder;
import structures.ListNum;
import tax.TaxNode;
import tax.TaxTree;
import tracker.ReadStats;

/**
 * This class does nothing.
 * It serves as a template for generating reads,
 * potentially based on some input sequence.
 * 
 * @author Brian Bushnell
 * @date June 8, 2019
 *
 */
public class RandomReadsMG {
	
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
		RandomReadsMG x=new RandomReadsMG(args);
		
		//Run the object
		x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	/**
	 * Constructor.
	 * @param args Command line arguments
	 */
	public RandomReadsMG(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		//Set shared static variables prior to parsing
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.setZipThreads(Shared.threads());
		FASTQ.TEST_INTERLEAVED=FASTQ.FORCE_INTERLEAVED=false;
		
		{//Parse the arguments
			final Parser parser=parse(args);
			Parser.processQuality();
			
			overwrite=ReadStats.overwrite=parser.overwrite;
			append=ReadStats.append=parser.append;

			out1=parser.out1;
			out2=parser.out2;
			qfout1=parser.qfout1;
			qfout2=parser.qfout2;
			extout=parser.extout;
		}

		validateParams();
		doPoundReplacement(); //Replace # with 1 and 2
		checkFileExistence(); //Ensure files can be read and written
		checkStatics(); //Adjust file-related static fields as needed for this program 
		
		//Create output FileFormat objects
		ffout1=FileFormat.testOutput(out1, FileFormat.FASTQ, extout, true, overwrite, append, false);
		ffout2=FileFormat.testOutput(out2, FileFormat.FASTQ, extout, true, overwrite, append, false);
		if("auto".equalsIgnoreCase(taxTreeFile)){taxTreeFile=TaxTree.defaultTreeFile();}
		tree=TaxTree.loadTaxTree(taxTreeFile, outstream, true, false);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------    Initialization Helpers    ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Parse arguments from the command line */
	private Parser parse(String[] args){
		
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
			}else if(a.equals("tree")){
				taxTreeFile=b;
			}else if(a.equals("in") || a.equals("ref")){
				Tools.getFileOrFiles(b, inputFiles, true, false, false, false);
			}else if(a.equals("depth") || a.equals("cov")){
				minDepth=maxDepth=Float.parseFloat(b);
			}else if(a.equals("mindepth") || a.equals("mincov")){
				minDepth=Float.parseFloat(b);
			}else if(a.equals("maxdepth") || a.equals("maxcov")){
				maxDepth=Float.parseFloat(b);
			}else if(a.equals("depthvariance")){
				depthVariance=Float.parseFloat(b);
			}else if(a.equals("insert") || a.equals("avginsert")){
				avgInsert=Float.parseFloat(b);
			}else if(a.equals("len") || a.equals("length") || a.equals("readlen") || a.equals("readlength")){
				readlen=Integer.parseInt(b);
			}else if(a.equals("paired") || a.equals("int") || a.equals("interleaved")){
				paired=Parse.parseBoolean(b);
			}else if(a.equals("depthvariance")){
				depthVariance=Float.parseFloat(b);
			}else if(a.equals("seed")){
				seed=Long.parseLong(b);
			}else if(a.equals("parse_flag_goes_here")){
				long fake_variable=Parse.parseKMG(b);
				//Set a variable here
			}else if(parser.parse(arg, a, b)){//Parse standard flags in the parser
				//do nothing
			}else{
				File f=new File(arg);
				if(f.exists() && f.canRead()) {
					Tools.getFileOrFiles(arg, inputFiles, true, false, false, false);
				}else {
					outstream.println("Unknown parameter "+args[i]);
					assert(false) : "Unknown parameter "+args[i];
				}
			}
		}
		
		return parser;
	}
	
	/** Replace # with 1 and 2 in headers */
	private void doPoundReplacement(){

		//Do output file # replacement
		if(out1!=null && out2==null && out1.indexOf('#')>-1){
			out2=out1.replace("#", "2");
			out1=out1.replace("#", "1");
		}

		//Ensure out2 is not set without out1
		if(out1==null && out2!=null){throw new RuntimeException("Error - cannot define out2 without defining out1.");}
	}
	
	/** Ensure files can be read and written */
	private void checkFileExistence(){
		//Ensure output files can be written
		if(!Tools.testOutputFiles(overwrite, append, false, out1, out2)){
			outstream.println((out1==null)+", "+(out2==null)+", "+out1+", "+out2);
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output files "+out1+", "+out2+"\n");
		}
		
		//Ensure input files can be read
		if(!Tools.testInputFiles(false, true, inputFiles.toArray(new String[0]))){
			throw new RuntimeException("\nCan't read some input files.\n");  
		}
		
		//Ensure that no file was specified multiple times
		if(!Tools.testForDuplicateFiles(true, out1, out2)){
			throw new RuntimeException("\nSome file names were specified multiple times.\n");
		}
	}
	
	/** Adjust file-related static fields as needed for this program */
	private static void checkStatics(){
		//Adjust the number of threads for input file reading
		if(!ByteFile.FORCE_MODE_BF1 && !ByteFile.FORCE_MODE_BF2 && Shared.threads()>2){
			ByteFile.FORCE_MODE_BF2=true;
		}
		
		assert(FastaReadInputStream.settingsOK());
	}
	
	/** Ensure parameter ranges are within bounds and required parameters are set */
	private boolean validateParams(){
//		assert(false) : "TODO";
		return true;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Create read streams and process all data */
	void process(Timer t){
		
		//Turn off read validation in the input threads to increase speed
		final boolean vic=Read.VALIDATE_IN_CONSTRUCTOR;
		Read.VALIDATE_IN_CONSTRUCTOR=true;
		
		//Optionally create a read output stream
		final ConcurrentReadOutputStream ros=makeCros();
		
		//Reset counters
		readsProcessed=readsOut=0;
		basesProcessed=basesOut=0;
		
		//Process the reads in separate threads
		spawnThreads(inputFiles, ros);
		
		if(verbose){outstream.println("Finished; closing streams.");}
		
		//Write anything that was accumulated by ReadStats
		errorState|=ReadStats.writeAll();
		//Close the read streams
		errorState|=ReadWrite.closeStream(ros);
		
		//Reset read validation
		Read.VALIDATE_IN_CONSTRUCTOR=vic;
		
		//Report timing and results
		t.stop();
		outstream.println(Tools.timeReadsBasesProcessed(t, readsProcessed, basesProcessed, 8));
		outstream.println(Tools.readsBasesOut(t.elapsed, readsOut, basesOut, 8));
		
		//Throw an exception of there was an error in a thread
		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}
	
	/** Create a Read Input Stream */
	private ConcurrentReadInputStream makeCris(FileFormat ff){
		ConcurrentReadInputStream cris=ConcurrentReadInputStream.getReadInputStream(-1, true, ff, null);
		cris.start(); //Start the stream
		if(verbose){outstream.println("Started cris");}
		return cris;
	}
	
	/** Create a Read Output Stream */
	private ConcurrentReadOutputStream makeCros(){
		if(ffout1==null){return null;}

		//Set output buffer size
		final int buff=4;

		//Notify user of output mode
		if(paired && out2==null){
			outstream.println("Writing interleaved.");
		}

		final ConcurrentReadOutputStream ros=ConcurrentReadOutputStream.getStream(
				ffout1, ffout2, qfout1, qfout2, buff, null, false);
		ros.start(); //Start the stream
		return ros;
	}
	
	/** Spawn process threads */
	private void spawnThreads(final ArrayList<String> files, final ConcurrentReadOutputStream ros){
		
		//Do anything necessary prior to processing
		
		//Determine how many threads may be used
		final int threads=Shared.threads();
		
		//Fill a list with ProcessThreads
		ArrayList<ProcessThread> alpt=new ArrayList<ProcessThread>(threads);
		AtomicInteger atom=new AtomicInteger(0);
		for(int i=0; i<threads; i++){
			alpt.add(new ProcessThread(files, ros, i, atom));
		}
		
		//Start the threads
		for(ProcessThread pt : alpt){
			pt.start();
		}
		
		//Wait for threads to finish
		waitForThreads(alpt);
		
		//Do anything necessary after processing
		
	}
	
	/** Wait until all worker threads are finished, then return */
	private void waitForThreads(ArrayList<ProcessThread> alpt){
		
		//Wait for completion of all threads
		boolean success=true;
		for(ProcessThread pt : alpt){
			
			//Wait until this thread has terminated
			while(pt.getState()!=Thread.State.TERMINATED){
				try {
					//Attempt a join operation
					pt.join();
				} catch (InterruptedException e) {
					//Potentially handle this, if it is expected to occur
					e.printStackTrace();
				}
			}
			
			//Accumulate per-thread statistics
			readsProcessed+=pt.readsInT;
			basesProcessed+=pt.basesInT;
			readsOut+=pt.readsOutT;
			basesOut+=pt.basesOutT;
			success&=pt.success;
		}
		
		//Track whether any threads failed
		if(!success){errorState=true;}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	private ConcurrentReadInputStream makeCris(String fname){
		FileFormat ff=FileFormat.testInput(fname, FileFormat.FASTA, null, true, true);
		ConcurrentReadInputStream cris=ConcurrentReadInputStream.getReadInputStream(-1, false, ff, null);
		cris.start(); //Start the stream
		if(verbose){outstream.println("Started cris");}
		return cris;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/
	
	private class ProcessThread extends Thread {
		
		//Constructor
		ProcessThread(final ArrayList<String> files_, final ConcurrentReadOutputStream ros_, 
				int tid_, final AtomicInteger nextFile_){
			files=files_;
			ros=ros_;
			tid=tid_;
			nextFile=nextFile_;
			
			minRoot=(float)Math.sqrt(minDepth);
			range=(float)(Math.sqrt(maxDepth)-minRoot);
		}
		
		//Called by start()
		@Override
		public void run(){
			//Do anything necessary prior to processing
			randy=Shared.threadLocalRandom(seed>=0 ? seed+tid : -1);
			
			//Process the files
			for(int i=nextFile.getAndIncrement(); i<files.size(); i=nextFile.getAndIncrement()) {
				String fname=files.get(i);
				processFile(fname, i);
			}
			
			//Do anything necessary after processing
			
			//Indicate successful exit status
			success=true;
		}
		
		/** Iterate through the reads */
		void processFile(String fname, int fnum){
//			System.err.println("Thread "+tid+" processing file "+fnum+"; next="+nextFile);
			ConcurrentReadInputStream cris=makeCris(fname);
				
			//Grab the first ListNum of reads
			ListNum<Read> ln=cris.nextList();
			
			final float rootDepth=minRoot+(Tools.min(randy.nextFloat(), randy.nextFloat(), 
					randy.nextFloat(), randy.nextFloat()))*range;
			final float fileDepth=rootDepth*rootDepth;
			int taxID=TaxTree.parseHeaderStatic2(fname, tree);
			if(taxID<0 && ln.size()>0) {
				Read c0=ln.get(0);
				taxID=TaxTree.parseHeaderStatic2(c0.id, tree);
			}
			System.err.println("File "+fnum+", tid "+taxID+": depth="+String.format("%.2f", fileDepth));
			assert(taxID>0) : "Can't parse taxID from "+fname;

			//As long as there is a nonempty read list...
			while(ln!=null && ln.size()>0){
				
				for(Read c : ln) {
					float depthRatio=1f+(depthVariance*(randy.nextFloat()+randy.nextFloat()))-depthVariance;
					float contigDepth=depthRatio*fileDepth;
//					System.err.println("depthRatio = "+depthRatio+"; contigDepth="+contigDepth);
					processContig(c, contigDepth, taxID, fnum);
				}

				//Notify the input stream that the list was used
				cris.returnList(ln);

				//Fetch a new list
				ln=cris.nextList();
			}

			//Notify the input stream that the final list was used
			if(ln!=null){
				cris.returnList(ln.id, ln.list==null || ln.list.isEmpty());
			}
			
		}
		
		private void processContig(Read contig, float depth, int taxID, int fnum) {
			final int basesPerRead=(paired ? 2*readlen : readlen);
			readsInT++;
			basesInT+=contig.length();
			if(contig.length()<basesPerRead+10) {return;}
			if(paired && contig.length()<avgInsert) {return;}
			
			long basesToGenerate=(long)(depth*contig.length());
			long readsGenerated=0;
			long basesGenerated=0;
			ArrayList<Read> list=new ArrayList<Read>(200);
			
//			System.err.println("Generating "+basesToGenerate+" for depth-"+depth+" contig "+contig.id);
			
			for(long i=0; basesGenerated<basesToGenerate; i++) {
				Read r=generateRead(contig, i, taxID, fnum, contig.numericID);
				if(r!=null) {
					list.add(r);
					readsGenerated+=r.pairCount();
					basesGenerated+=r.pairLength();
				}
				if(list.size()>=200) {
					if(ros!=null) {ros.add(list, 0);}
					list=new ArrayList<Read>(200);
				}
			}
			if(list.size()>0) {if(ros!=null) {ros.add(list, 0);}}
//			System.err.println("Generated "+basesGenerated+" for depth-"+depth+" contig "+contig.id);
			
			readsOutT+=readsGenerated;
			basesOutT+=basesGenerated;
		}
		
		private Read generateRead(Read contig, long rnum, int taxID, int fnum, long cnum) {
			if(paired) {return generateReadPair(contig, rnum, taxID, fnum, cnum);}
			else {return generateReadSingle(contig, rnum, taxID, fnum, cnum);}
		}
		
		private Read generateReadSingle(Read contig, long rnum, int taxID, int fnum, long cnum) {
			int insert=readlen;
			int start=randy.nextInt(contig.length()-insert);
			int strand=randy.nextInt(2);
			byte[] bases=Arrays.copyOfRange(contig.bases, start, start+readlen);
			if(strand==1) {AminoAcid.reverseComplementBasesInPlace(bases);}
			String header=makeHeader(start, strand, insert, taxID, fnum, cnum, 0);
			Read r=new Read(bases, null, header, rnum);
			return r;
		}
		
		private Read generateReadPair(Read contig, long rnum, int taxID, int fnum, long cnum) {
			double g=randy.nextGaussian()*0.25f;
			int insert=(int)((1+g)*avgInsert);
			while(insert>=contig.length() || insert<readlen) {
				g=randy.nextGaussian()*0.25f;
				insert=(int)((1+g)*avgInsert);
			}
			int start1=randy.nextInt(contig.length()-insert);
			int strand=randy.nextInt(2);
			int start2=start1+insert-readlen;

			byte[] bases1=Arrays.copyOfRange(contig.bases, start1, start1+readlen);
			byte[] bases2=Arrays.copyOfRange(contig.bases, start2, start2+readlen);
			AminoAcid.reverseComplementBasesInPlace(bases2);
			if(strand==1) {
				byte[] temp=bases1;
				bases1=bases2;
				bases2=temp;
			}
			String header1=makeHeader(start1, strand, insert, taxID, fnum, cnum, 0);
			String header2=makeHeader(start1, strand, insert, taxID, fnum, cnum, 1);
			Read r1=new Read(bases1, null, header1, rnum);
			Read r2=new Read(bases2, null, header2, rnum);
			r2.setPairnum(1);
			r1.mate=r2;
			r2.mate=r1;
			return r1;
		}
		
		private String makeHeader(int start, int strand, int insert, int taxID, int fnum, long cnum, int pnum) {
			bb.clear().append('f').under().append(fnum).under().append('c').under().append(cnum);
			bb.under().append('s').under().append(strand).under().append('p').under().append(start);
			bb.under().append('i').under().append(insert).under().append("tid");
			bb.under().append(taxID).space().append(pnum+1).colon();
			return bb.toString();
		}
		
		/** Number of reads retained by this thread */
		protected long readsOutT=0;
		/** Number of bases retained by this thread */
		protected long basesOutT=0;
		
		/** Number of input contigs processed by this thread */
		protected long readsInT=0;
		/** Number of input bases processed by this thread */
		protected long basesInT=0;
		
		/** True only if this thread has completed successfully */
		boolean success=false;
		
		/** Shared output stream */
		private final ConcurrentReadOutputStream ros;
		/** Thread ID */
		final int tid;
		final AtomicInteger nextFile;
		private final ArrayList<String> files;
		private Random randy;
		
		private final float minRoot;
		private final float range;
		private ByteBuilder bb=new ByteBuilder(128);
	}
	
	private class Genome {
		
		public Genome(String fname_) {
			fname=fname_;
		}
		
		private void load() {
			FileFormat ff=FileFormat.testInput(fname, FileFormat.FASTA, null, true, true);
			contigs=ConcurrentReadInputStream.getReads(-1, false, ff, null, null, null);
			taxID=TaxTree.parseHeaderStatic2(fname, tree);
			if(taxID<0) {taxID=TaxTree.parseHeaderStatic2(contigs.get(0).id, tree);}
		}
		
		ArrayList<Read> contigs;
		String fname;
		int taxID=-1;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Primary input file path */
	private ArrayList<String> inputFiles=new ArrayList<String>();

	/** Primary output file path */
	private String out1=null;
	/** Secondary output file path */
	private String out2=null;

	private String qfout1=null;
	private String qfout2=null;
	
	/** Override output file extension */
	private String extout=null;
	
	private String taxTreeFile=null;
	
	/*--------------------------------------------------------------*/

	/** Number of reads processed */
	protected long readsProcessed=0;
	/** Number of bases processed */
	protected long basesProcessed=0;

	/** Number of reads retained */
	protected long readsOut=0;
	/** Number of bases retained */
	protected long basesOut=0;
	
	private AtomicLong nextReadID=new AtomicLong(0);
	
	private float minDepth=1;
	private float maxDepth=256;
	private float depthVariance=0.3f;
	private float avgInsert=300;
	private int readlen=150;
	private boolean paired=true;
	private long seed=-1;
	
	/*--------------------------------------------------------------*/
	/*----------------         Final Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Primary output file */
	private final FileFormat ffout1;
	/** Secondary output file */
	private final FileFormat ffout2;
	
	private final TaxTree tree;
	
	/*--------------------------------------------------------------*/
	/*----------------        Common Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Print status messages to this output stream */
	private PrintStream outstream=System.err;
	/** Print verbose messages */
	public static boolean verbose=false;
	/** True if an error was encountered */
	public boolean errorState=false;
	/** Overwrite existing output files */
	private boolean overwrite=true;
	/** Append to existing output files */
	private boolean append=false;
	
}
