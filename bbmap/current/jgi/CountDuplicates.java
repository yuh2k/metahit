package jgi;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Locale;
import java.util.Random;

import fileIO.ByteFile;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import kmer.HashArray1D;
import kmer.HashBuffer;
import kmer.KmerTableSet;
import shared.Parse;
import shared.Parser;
import shared.PreParser;
import tracker.ReadStats;
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
import structures.LongList;
import template.Accumulator;
import template.ThreadWaiter;

/**
 * Counts duplicate sequences probabilistically,
 * using around 20 bytes per unique read.  Read pairs are treated
 * as a single read.  Reads are converted to a hashcode and only
 * the hashcode is stored when tracking duplicates, so (rare) hash
 * collisions will result in false positive duplicate detection.
 * Optionally outputs the deduplicated and/or duplicate reads.
 * 
 * @author Brian Bushnell
 * @date November 19, 2022
 *
 */
public class CountDuplicates implements Accumulator<CountDuplicates.ProcessThread> {
	
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
		CountDuplicates x=new CountDuplicates(args);
		
		//Run the object
		x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	/**
	 * Constructor.
	 * @param args Command line arguments
	 */
	public CountDuplicates(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		//Set shared static variables prior to parsing
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.MAX_ZIP_THREADS=Shared.threads();
		
		{//Parse the arguments
			final Parser parser=parse(args);
			Parser.processQuality();
			
			maxReads=parser.maxReads;
			overwrite=ReadStats.overwrite=parser.overwrite;
			append=ReadStats.append=parser.append;
			setInterleaved=parser.setInterleaved;
			
			in1=parser.in1;
			in2=parser.in2;
			qfin1=parser.qfin1;
			qfin2=parser.qfin2;
			extin=parser.extin;

			out1=parser.out1;
			out2=parser.out2;
			qfout1=parser.qfout1;
			qfout2=parser.qfout2;
			extout=parser.extout;
		}
		
		Shared.capThreads(20); //Does not need more than ~6 threads.  So this accounts for 4-way SMT.

		validateParams();
		doPoundReplacement(); //Replace # with 1 and 2
		adjustInterleaving(); //Make sure interleaving agrees with number of input and output files
		fixExtensions(); //Add or remove .gz or .bz2 as needed
		checkFileExistence(); //Ensure files can be read and written
		checkStatics(); //Adjust file-related static fields as needed for this program 
		
		//Create output FileFormat objects
		ffout1=FileFormat.testOutput(out1, FileFormat.FASTQ, extout, true, overwrite, append, ordered);
		ffout2=FileFormat.testOutput(out2, FileFormat.FASTQ, extout, true, overwrite, append, ordered);
		
		//Create output FileFormat objects
		ffoutd1=FileFormat.testOutput(outd1, FileFormat.HEADER, null, true, overwrite, append, ordered);
		ffoutd2=FileFormat.testOutput(outd2, FileFormat.HEADER, null, true, overwrite, append, ordered);

		//Create input FileFormat objects
		ffin1=FileFormat.testInput(in1, FileFormat.FASTQ, extin, true, true);
		ffin2=FileFormat.testInput(in2, FileFormat.FASTQ, extin, true, true);
		
		sampleThresh=(int)(samplerate*(sampleMask+1));
		
		Parser.silent=PreParser.silent=true;
		tables=new KmerTableSet(args, 12);
		tables.allocateTables();
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
			}else if(a.equals("ordered")){
				ordered=Parse.parseBoolean(b);
			}else if(a.equals("bases") || a.equals("hashbases")){
				hashBases=Parse.parseBoolean(b);
			}else if(a.equals("names") || a.equals("headers") || a.equals("hashnames")){
				hashNames=Parse.parseBoolean(b);
			}else if(a.equals("qualities") || a.equals("hashqualities") || a.equals("quals")){
				hashQualities=Parse.parseBoolean(b);
			}else if(a.equals("stats") || a.equals("outstats")){
				stats=b;
			}else if(a.equals("maxfraction")){
				maxFraction=Double.parseDouble(b);
			}else if(a.equals("maxrate")){
				maxRate=Double.parseDouble(b);
			}else if(a.equals("exitcode") || a.equals("failcode")){
				failCode=Integer.parseInt(b);
			}else if(a.equals("samplerate")){
				samplerate=Double.parseDouble(b);
			}else if(a.equals("outd") || a.equals("dupes") || a.equals("outd1") || a.equals("dupes1")){
				outd1=b;
			}else if(a.equals("outd2") || a.equals("dupes2")){
				outd2=b;
			}else if(a.equals("parse_flag_goes_here")){
				long fake_variable=Parse.parseKMG(b);
				//Set a variable here
			}else if(parser.parse(arg, a, b)){//Parse standard flags in the parser
				//do nothing
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
			}
		}
		
		return parser;
	}
	
	/** Replace # with 1 and 2 in headers */
	private void doPoundReplacement(){
		//Do input file # replacement
		if(in1!=null && in2==null && in1.indexOf('#')>-1 && !new File(in1).exists()){
			in2=in1.replace("#", "2");
			in1=in1.replace("#", "1");
		}

		//Do output file # replacement
		if(out1!=null && out2==null && out1.indexOf('#')>-1){
			out2=out1.replace("#", "2");
			out1=out1.replace("#", "1");
		}

		//Do output file # replacement
		if(outd1!=null && outd2==null && outd1.indexOf('#')>-1){
			outd2=outd1.replace("#", "2");
			outd1=outd1.replace("#", "1");
		}
		
		//Ensure there is an input file
		if(in1==null){throw new RuntimeException("Error - at least one input file is required.");}

		//Ensure out2 is not set without out1
		if(out1==null && out2!=null){throw new RuntimeException("Error - cannot define out2 without defining out1.");}
	}
	
	/** Add or remove .gz or .bz2 as needed */
	private void fixExtensions(){
		in1=Tools.fixExtension(in1);
		in2=Tools.fixExtension(in2);
		qfin1=Tools.fixExtension(qfin1);
		qfin2=Tools.fixExtension(qfin2);
	}
	
	/** Ensure files can be read and written */
	private void checkFileExistence(){
		//Ensure output files can be written
		if(!Tools.testOutputFiles(overwrite, append, false, out1, out2, outd1, outd2)){
			outstream.println((out1==null)+", "+(out2==null)+", "+out1+", "+out2);
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output files "+out1+", "+out2+"\n");
		}
		
		//Ensure input files can be read
		if(!Tools.testInputFiles(false, true, in1, in2)){
			throw new RuntimeException("\nCan't read some input files.\n");  
		}
		
		//Ensure that no file was specified multiple times
		if(!Tools.testForDuplicateFiles(true, in1, in2, out1, out2, outd1, outd2)){
			throw new RuntimeException("\nSome file names were specified multiple times.\n");
		}
	}
	
	/** Make sure interleaving agrees with number of input and output files */
	private void adjustInterleaving(){
		//Adjust interleaved detection based on the number of input files
		if(in2!=null){
			if(FASTQ.FORCE_INTERLEAVED){outstream.println("Reset INTERLEAVED to false because paired input files were specified.");}
			FASTQ.FORCE_INTERLEAVED=FASTQ.TEST_INTERLEAVED=false;
		}

		//Adjust interleaved settings based on number of output files
		if(!setInterleaved){
			assert(in1!=null && (out1!=null || out2==null)) : "\nin1="+in1+"\nin2="+in2+"\nout1="+out1+"\nout2="+out2+"\n";
			if(in2!=null){ //If there are 2 input streams.
				FASTQ.FORCE_INTERLEAVED=FASTQ.TEST_INTERLEAVED=false;
				outstream.println("Set INTERLEAVED to "+FASTQ.FORCE_INTERLEAVED);
			}else{ //There is one input stream.
				if(out2!=null){
					FASTQ.FORCE_INTERLEAVED=true;
					FASTQ.TEST_INTERLEAVED=false;
					outstream.println("Set INTERLEAVED to "+FASTQ.FORCE_INTERLEAVED);
				}
			}
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
		assert(samplerate>=0 && samplerate<=1) : "Invalid samplerate: Range is 0-1.";
		assert(maxFraction<0 || maxFraction<=1) : "Invalid maxFraction: Range is negative, or 0-1.";
		assert(maxRate<0 || maxRate>=1) : "Invalid maxRate: Range is negative, or 1+.";
		return true;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Create read streams and process all data */
	void process(Timer t){
		
		//Turn off read validation in the input threads to increase speed
		final boolean vic=Read.VALIDATE_IN_CONSTRUCTOR;
		Read.VALIDATE_IN_CONSTRUCTOR=Shared.threads()<4;
		
		//Create a read input stream
		final ConcurrentReadInputStream cris=makeCris();
		
		//Optionally create a read output stream
		final ConcurrentReadOutputStream ros=makeCros(cris.paired());
		
		//Optionally create a duplicate read output stream
		final ConcurrentReadOutputStream rosd=makeCrosD(cris.paired());
		
		//Reset counters
		readsProcessed=readsOut=0;
		basesProcessed=basesOut=0;
		
		//Process the reads in separate threads
		spawnThreads(cris, ros, rosd);
		
		if(verbose){outstream.println("Finished; closing streams.");}
		
		System.err.println();
		calcStatistics(stats);
		
		//Write anything that was accumulated by ReadStats
		errorState|=ReadStats.writeAll();
		//Close the read streams
		errorState|=ReadWrite.closeStreams(cris, ros, rosd);
		
		//Reset read validation
		Read.VALIDATE_IN_CONSTRUCTOR=vic;
		
		//Report timing and results
		t.stop();
		outstream.println(Tools.timeReadsBasesProcessed(t, readsProcessed, basesProcessed, 8));
		outstream.println(Tools.readsBasesOut(readsProcessed, basesProcessed, readsOut, basesOut, 8, false));
		
		//Throw an exception of there was an error in a thread
		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
		
		if(failed){
			System.err.println("FAIL");
			System.exit(failCode);
		}else if(maxFraction>=0 || maxRate>=0){
			System.err.println("PASS");
		}
	}
	
	private ConcurrentReadInputStream makeCris(){
		ConcurrentReadInputStream cris=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ffin1, ffin2, qfin1, qfin2);
		cris.start(); //Start the stream
		if(verbose){outstream.println("Started cris");}
		boolean paired=cris.paired();
		if(!ffin1.samOrBam()){outstream.println("Input is being processed as "+(paired ? "paired" : "unpaired"));}
		return cris;
	}
	
	private ConcurrentReadOutputStream makeCros(boolean pairedInput){
		if(ffout1==null){return null;}

		//Select output buffer size based on whether it needs to be ordered
		final int buff=(ordered ? Tools.mid(16, 128, (Shared.threads()*2)/3) : 8);

		//Notify user of output mode
		if(pairedInput && out2==null && (in1!=null && !ffin1.samOrBam() && !ffout1.samOrBam())){
			outstream.println("Writing interleaved.");
		}

		final ConcurrentReadOutputStream ros=ConcurrentReadOutputStream.getStream(ffout1, ffout2, qfout1, qfout2, buff, null, false);
		ros.start(); //Start the stream
		return ros;
	}
	
	private ConcurrentReadOutputStream makeCrosD(boolean pairedInput){
		if(ffoutd1==null){return null;}

		//Select output buffer size based on whether it needs to be ordered
		final int buff=(ordered ? Tools.mid(16, 128, (Shared.threads()*2)/3) : 8);

		final ConcurrentReadOutputStream rosd=ConcurrentReadOutputStream.getStream(ffoutd1, ffoutd2, null, null, buff, null, false);
		rosd.start(); //Start the stream
		return rosd;
	}
	
	private void calcStatistics(String fname){
		ByteBuilder bb=new ByteBuilder();
		
		long[] counts=tables.fillHistogram(4000);
		long sumTotal=Tools.sumHistogram(counts);
		long sumUnique=Tools.sum(counts);
		
		long[] multiplicity=new long[10];
		multiplicity[1]=counts[1];
		for(int i=1; i<counts.length; i++){
			long count=counts[i];
			for(int j=2; j<multiplicity.length; j++){
				if(i%j==0){multiplicity[j]+=count*i;}
			}
		}
		final double mult=(sumTotal<1 ? 0 : 1.0/sumTotal);
		final long x1=counts[1];
		double fraction=mult*(sumTotal-x1);
		double rate=sumTotal*1.0/sumUnique;
		bb.appendln("#Fraction of reads with duplicates:"+
				String.format(Locale.ROOT, "\t%.4f", fraction));
		bb.appendln("#Average duplication rate:         "+
				String.format(Locale.ROOT, "\t%.4f", rate));
		bb.appendln("#Read Copy Count Distribution");
		bb.appendln("#Copies\tCount\tFraction");
		bb.appendln("1\t"+x1+String.format(Locale.ROOT, "\t%.4f", mult*x1));
		for(int i=2; i<multiplicity.length; i++){
			long x=multiplicity[i];
			bb.appendln(i+"N\t"+x+String.format(Locale.ROOT, "\t%.4f", mult*x));
		}
		if(fname!=null){
			ReadWrite.writeString(bb, fname);
		}

		if(maxFraction>=0 && fraction>maxFraction){failed=true;}
		if(maxRate>=0 && rate>maxRate){failed=true;}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------       Thread Management      ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Spawn process threads */
	private void spawnThreads(final ConcurrentReadInputStream cris, 
			final ConcurrentReadOutputStream ros, final ConcurrentReadOutputStream rosd){
		
		//Do anything necessary prior to processing
		
		//Determine how many threads may be used
		final int threads=Shared.threads();
		
		//Fill a list with ProcessThreads
		ArrayList<ProcessThread> alpt=new ArrayList<ProcessThread>(threads);
		for(int i=0; i<threads; i++){
			alpt.add(new ProcessThread(cris, ros, rosd, i));
		}
		
		//Start the threads and wait for them to finish
		boolean success=ThreadWaiter.startAndWait(alpt, this);
		errorState&=!success;
		
		//Do anything necessary after processing
		
	}
	
	@Override
	public final void accumulate(ProcessThread pt){
		readsProcessed+=pt.readsProcessedT;
		basesProcessed+=pt.basesProcessedT;
		readsOut+=pt.readsOutT;
		basesOut+=pt.basesOutT;
		errorState|=(!pt.success);
	}
	
	@Override
	public final boolean success(){return !errorState;}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/
	
	class ProcessThread extends Thread {

		//Constructor
		ProcessThread(final ConcurrentReadInputStream cris_, 
				final ConcurrentReadOutputStream ros_, final ConcurrentReadOutputStream rosd_, final int tid_){
			cris=cris_;
			ros=ros_;
			rosd=rosd_;
			tid=tid_;
		}

		//Called by start()
		@Override
		public void run(){
			//Do anything necessary prior to processing

			//Process the reads
			processInner();

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
			final ArrayList<Read> dupes=new ArrayList<Read>(reads.size());
			
			//Generate hashcodes and place them in counts array.
			//This is the slow bit so it goes outside the lock.
			hashList(reads);
			
			//This involves the lock
			int numDupes=processCodes();
			
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
				
				long code=codes.get(idx);
				if(code>=0){
					readsOutT+=r1.pairCount();
					basesOutT+=r1.pairLength();
				}else{
					reads.set(idx, null);
					if(code==DUPLICATE){dupes.add(r1);}
				}
			}

			//Output reads to the output stream
			if(ros!=null){ros.add(reads, ln.id);}
			if(rosd!=null){rosd.add(dupes, ln.id);}
		}
		
		//Critical section
		private int processCodes(){
			int dupes=0;
			synchronized(tables){
				for(int i=0; i<codes.size; i++){
					long code=codes.get(i);
					if(code>=0){
						HashArray1D table=tables.getTableForKey(code);
						int count=table.increment(code, 1);
						if(count>maxCount){
							codes.set(i, DUPLICATE);
							dupes++;
						}
					}
				}
			}
			return dupes;
		}

		private void hashList(ArrayList<Read> reads){
			codes.clear();
			for(final Read r1 : reads){
				codes.add(r1==null ? 0 : hashRead(r1, r1.mate));
			}
		}
		
		long hashRead(final Read r1, final Read r2){
			long code=(hashBases ? hashBases(r1, r2) : 0);
			code^=(hashNames ? hashNames(r1, r2) : 0);
			code^=(hashQualities ? hashQualities(r1, r2) : 0);
			if(samplerate<1 && (code&sampleMask)<sampleThresh){return UNSAMPLED;}
			return code&Long.MAX_VALUE;
		}
		
		private long hashBases(final Read r1, final Read r2){
			long code=hash(r1.bases, 0);
			if(r2!=null){code=code^Long.rotateRight(hash(r2.bases, 1), 3);}
			return code;
		}
		
		private long hashNames(final Read r1, final Read r2){
			long code=hash(r1.id.getBytes(), 7);
			return code;
		}
		
		private long hashQualities(final Read r1, final Read r2){
			long code=hash(r1.quality, 0);
			if(r2!=null){code=code^Long.rotateRight(hash(r2.quality, 1), 3);}
			return code;
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
		/** Shared output stream for duplicate reads */
		private final ConcurrentReadOutputStream rosd;
		/** Thread ID */
		final int tid;
		
		final LongList codes=new LongList(Shared.bufferLen());
	}
	
	public static long hash(byte[] bytes, int pairnum){
		if(bytes==null){return 0;}
		long code=(bytes.length+pairnum)^salt;
		for(int i=0; i<bytes.length; i++){
			byte b=bytes[i];
			int mode=(int)(code&31);
			assert(hashcodes[b]!=null) : "Invalid sequence character: '"+(char)b+"'";
			code=code^hashcodes[b][mode];
			code=Long.rotateLeft(code, 1);
		}
		return code;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Primary input file path */
	private String in1=null;
	/** Secondary input file path */
	private String in2=null;
	
	private String qfin1=null;
	private String qfin2=null;

	/** Primary output file path */
	private String out1=null;
	/** Secondary output file path */
	private String out2=null;

	private String qfout1=null;
	private String qfout2=null;

	/** Primary dupe output file path */
	private String outd1=null;
	/** Secondary dupe output file path */
	private String outd2=null;
	
	/** Override input file extension */
	private String extin=null;
	/** Override output file extension */
	private String extout=null;

	/** Stats output */
	private String stats="stdout.txt";
	
	private double maxFraction=-1;
	private double maxRate=-1;
	private int failCode=0;
	boolean failed=false;
	
	/** Whether interleaved was explicitly set. */
	private boolean setInterleaved=false;

	private KmerTableSet tables;

	private boolean hashBases=true;
	private boolean hashNames=false;
	private boolean hashQualities=false;
	
	/*--------------------------------------------------------------*/

	/** Number of reads processed */
	protected long readsProcessed=0;
	/** Number of bases processed */
	protected long basesProcessed=0;

	/** Number of reads retained */
	protected long readsOut=0;
	/** Number of bases retained */
	protected long basesOut=0;

	/** Quit after processing this many input reads; -1 means no limit */
	private long maxReads=-1;
	private double samplerate=1;
	private int maxCount=1;
	
	/*--------------------------------------------------------------*/
	/*----------------         Final Fields         ----------------*/
	/*--------------------------------------------------------------*/

	/** Primary input file */
	private final FileFormat ffin1;
	/** Secondary input file */
	private final FileFormat ffin2;
	
	/** Primary output file */
	private final FileFormat ffout1;
	/** Secondary output file */
	private final FileFormat ffout2;
	
	/** Primary dupe output file */
	private final FileFormat ffoutd1;
	/** Secondary dupe output file */
	private final FileFormat ffoutd2;
	
	private static final long[][] hashcodes=Dedupe.makeCodes2(32);
	private static final long salt=new Random(173).nextLong();
	private final int sampleMask=1023;
	private final int sampleThresh;
	
	private final int UNSAMPLED=-1;
	private final int DUPLICATE=-2;
	
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
	private boolean overwrite=false;
	/** Append to existing output files */
	private boolean append=false;
	/** Reads are output in input order */
	private boolean ordered=false;
	
}
