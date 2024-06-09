package barcode;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;

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
import stream.FASTQ;
import stream.FastaReadInputStream;
import stream.Read;
import structures.ListNum;
import tracker.ReadStats;

/**
 * Counts barcodes in a fastq file and produces a summary.
 * 
 * @author Brian Bushnell
 * @date June 20, 2014
 *
 */
public class CountBarcodes2 {
	
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
		CountBarcodes2 x=new CountBarcodes2(args);
		
		//Run the object
		x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	/**
	 * Constructor.
	 * @param args Command line arguments
	 */
	public CountBarcodes2(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		//Set shared static variables prior to parsing
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.MAX_ZIP_THREADS=Shared.threads();
		Shared.capBuffers(4); //Only for singlethreaded programs
		
		{//Parse the arguments
			final Parser parser=parse(args);
			Parser.processQuality();
			
			maxReads=parser.maxReads;
			overwrite=ReadStats.overwrite=parser.overwrite;
			append=ReadStats.append=parser.append;
			setInterleaved=parser.setInterleaved;

			if(parser.in1!=null) {
				for(String s : parser.in1.split(",")) {
					in1.add(s);
				}
			}
			if(parser.in2!=null) {
				for(String s : parser.in2.split(",")) {
					in2.add(s);
				}
			}
			extin=parser.extin;

			out1=parser.out1;
			extout=parser.extout;
		}

		doPoundReplacement(); //Replace # with 1 and 2
		adjustInterleaving(); //Make sure interleaving agrees with number of input and output files
		fixExtensions(); //Add or remove .gz or .bz2 as needed
		checkFileExistence(); //Ensure files can be read and written
		checkStatics(); //Adjust file-related static fields as needed for this program 
		
		//Create output FileFormat objects
		ffout1=FileFormat.testOutput(out1, FileFormat.TXT, extout, true, overwrite, append, false);

		//Create input FileFormat objects
		for(String s : in1) {
			ffin1.add(FileFormat.testInput(s, FileFormat.FASTQ, extin, true, true));
		}
		for(String s : in2) {
			ffin2.add(FileFormat.testInput(s, FileFormat.FASTQ, extin, true, true));
		}
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
			}else if(a.equals("expected") || a.equals("valid") || 
					a.equals("barcodes") || a.equals("expectedbarcodes")){
				expectedBarcodesFile=b;
			}else if(a.equals("delimiter")){
				delimiter=(b==null ? -1 : (byte)b.charAt(0));
				barcodesPerRead=(delimiter<0 ? 1 : 2);
			}else if(a.equals("parse_flag_goes_here")){
				long fake_variable=Parse.parseKMG(b);
				//Set a variable here
			}else if(parser.parse(arg, a, b)){//Parse standard flags in the parser
				//do nothing
			}else if(b==null && new File(arg).exists()){
				in1.add(arg);
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
		if(in2.isEmpty()) {
			for(int i=0; i<in1.size(); i++) {
				String s=in1.get(i);
				if(s.indexOf('#')>-1 && !new File(s).exists()){
					String a=s.replace("#", "1");
					String b=s.replace("#", "2");
					in1.set(i, a);
					in2.add(b);
				}
			}
		}
		
		//Ensure there is an input file
		if(in1.isEmpty()){throw new RuntimeException("Error - at least one input file is required.");}
	}
	
	/** Add or remove .gz or .bz2 as needed */
	private void fixExtensions(){
		in1=Tools.fixExtension(in1);
		in2=Tools.fixExtension(in2);
	}
	
	/** Ensure files can be read and written */
	private void checkFileExistence(){
		//Ensure output files can be written
		if(!Tools.testOutputFiles(overwrite, append, false, out1)){
			outstream.println((out1==null)+", "+out1);
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output file "+out1+"\n");
		}
		
		ArrayList<String> temp=new ArrayList<String>();
		temp.addAll(in1);
		temp.addAll(in2);
		if(expectedBarcodesFile!=null) {temp.add(expectedBarcodesFile);}
		
		
		//Ensure input files can be read
		if(!Tools.testInputFiles(false, true, temp.toArray(new String[0]))){
			throw new RuntimeException("\nCan't read some input files.\n");  
		}
		
		if(out1!=null) {temp.add(out1);}
		//Ensure that no file was specified multiple times
		if(!Tools.testForDuplicateFiles(true, temp.toArray(new String[0]))){
			throw new RuntimeException("\nSome file names were specified multiple times.\n");
		}
	}
	
	/** Make sure interleaving agrees with number of input and output files */
	private void adjustInterleaving(){
		//Adjust interleaved detection based on the number of input files
		if(!in2.isEmpty()){
			if(FASTQ.FORCE_INTERLEAVED){outstream.println("Reset INTERLEAVED to false because paired input files were specified.");}
			FASTQ.FORCE_INTERLEAVED=false; //FASTQ.TEST_INTERLEAVED=false; interferes with barcode detection
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
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Create read streams and process all data */
	void process(Timer t){
		
		//Reset counters
		readsProcessed=0;
		basesProcessed=0;
		
		//Process the read stream
		for(int i=0; i<ffin1.size(); i++){
			processInner(ffin1.get(i), (ffin2.size()<=i ? null : ffin2.get(i)));
		}

		if(calcStats) {
			if(bs.barcodesPerRead>1) {
				bs.leftStats=bs.makeLeft();
				bs.rightStats=bs.makeRight();
				bs.leftStats.calcStats();
				bs.rightStats.calcStats();
			}
			bs.calcStats();

			String s=bs.toStats("Barcodes");
			System.err.println(s);
			if(bs.leftStats!=null){
				System.err.println(bs.leftStats.toStats("\nBarcode1"));
			}
			if(bs.rightStats!=null){
				System.err.println(bs.rightStats.toStats("\nBarcode2"));
			}
		}
		
		if(ffout1!=null){
			bs.printToFile(ffout1);
		}
		
		if(verbose){outstream.println("Finished; closing streams.");}
		
		//Write anything that was accumulated by ReadStats
		errorState|=ReadStats.writeAll();
		
		//Report timing and results
		t.stop();
		outstream.println();
		outstream.println(Tools.timeReadsBasesProcessed(t, readsProcessed, basesProcessed, 8));
		
		//Throw an exception of there was an error in a thread
		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	private ConcurrentReadInputStream makeCris(FileFormat ff1, FileFormat ff2){
		ConcurrentReadInputStream cris=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ff1, ff2, null, null);
		cris.start(); //Start the stream
		if(verbose){outstream.println("Started cris");}
		boolean paired=cris.paired();
		if(!ff1.samOrBam()){outstream.println("Input is being processed as "+(paired ? "paired" : "unpaired"));}
		return cris;
	}
	
	/** Iterate through the reads */
	void processInner(final FileFormat ff1, final FileFormat ff2){
		
		//Do anything necessary prior to processing
		
		//Create a read input stream
		final ConcurrentReadInputStream cris=makeCris(ff1, ff2);
		
		if(bs==null) {
			if(delimiter<0) {
				delimiter=(byte)ff1.barcodeDelimiter();
				barcodesPerRead=ff1.barcodesPerRead();
			}
			bs=new BarcodeStats(delimiter, barcodesPerRead, extin);
			if(expectedBarcodesFile!=null) {bs.loadBarcodeList(expectedBarcodesFile);}
		}
		readFile(cris);
		
		//Close the read streams
		errorState|=ReadWrite.closeStreams(cris);
		
		//Do anything necessary after processing
		
	}
	
	void readFile(final ConcurrentReadInputStream cris) {
		{
			//Grab the first ListNum of reads
			ListNum<Read> ln=cris.nextList();

			//Check to ensure pairing is as expected
			if(ln!=null && !ln.isEmpty()){
				Read r=ln.get(0);
				assert(r.samline!=null || (r.mate!=null)==cris.paired());
			}

			//As long as there is a nonempty read list...
			while(ln!=null && ln.size()>0){
//				if(verbose){outstream.println("Fetched "+reads.size()+" reads.");} //Disabled due to non-static access
				
				processList(ln, cris);

				//Fetch a new list
				ln=cris.nextList();
			}
			
			//Notify the input stream that the final list was used
			if(ln!=null){
				cris.returnList(ln.id, ln.list==null || ln.list.isEmpty());
			}
		}
	}
	
	/**
	 * Process a list of Reads.
	 * @param ln The list.
	 * @param cris Read Input Stream
	 * @param ros Read Output Stream for reads that will be retained
	 */
	void processList(ListNum<Read> ln, final ConcurrentReadInputStream cris){

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
			readsProcessed+=r1.pairCount();
			basesProcessed+=initialLength1+initialLength2;
			
			{
				//Reads are processed in this block.
				processReadPair(r1, r2);
			}
		}

		//Notify the input stream that the list was used
		cris.returnList(ln);
//		if(verbose){outstream.println("Returned a list.");} //Disabled due to non-static access
	}
	
	
	/**
	 * Process a single read pair.
	 * @param r1 Read 1
	 * @param r2 Read 2 (may be null)
	 * @return True if the reads should be kept, false if they should be discarded.
	 */
	void processReadPair(final Read r1, final Read r2){
		String code=r1.barcode(true);
		if(code==null){return;}
		bs.increment(code, 1);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Primary input file path */
	private ArrayList<String> in1=new ArrayList<String>();
	/** Secondary input file path */
	private ArrayList<String> in2=new ArrayList<String>();
	
	/** Expected Barcodes */
	private String expectedBarcodesFile=null;

	/** Primary output file path */
	private String out1=null;
	
	/** Override input file extension */
	private String extin=null;
	/** Override output file extension */
	private String extout=null;
	
	/** Whether interleaved was explicitly set. */
	private boolean setInterleaved=false;
	
	private BarcodeStats bs;
	private boolean calcStats=true;
	
	byte delimiter=-1;
	int barcodesPerRead=-1;
	
	/*--------------------------------------------------------------*/

	/** Number of reads processed */
	protected long readsProcessed=0;
	/** Number of bases processed */
	protected long basesProcessed=0;

	/** Quit after processing this many input reads; -1 means no limit */
	private long maxReads=-1;
	
	/*--------------------------------------------------------------*/
	/*----------------         Final Fields         ----------------*/
	/*--------------------------------------------------------------*/

	/** Primary input file */
	private final ArrayList<FileFormat> ffin1=new ArrayList<FileFormat>();
	/** Secondary input file */
	private final ArrayList<FileFormat> ffin2=new ArrayList<FileFormat>();
	
	/** Primary output file */
	private final FileFormat ffout1;
	
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
	
}
