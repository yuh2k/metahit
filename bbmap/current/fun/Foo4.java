package fun;

import java.io.PrintStream;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.TimeZone;

import fileIO.ByteFile;
import fileIO.ByteFile1;
import fileIO.ByteFile2;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import shared.Parse;
import shared.Parser;
import shared.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import structures.LongList;

/**
 * Reads a text file.
 * Prints it to another text file.
 * Filters out invalid lines and prints them to an optional third file.
 * @author Brian Bushnell
 * @date May 9, 2016
 *
 */
public class Foo4 {
	
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
		Foo4 x=new Foo4(args);
		
		//Run the object
		x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	/**
	 * Constructor.
	 * @param args Command line arguments
	 */
	public Foo4(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, /*getClass()*/null, false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		//Set shared static variables prior to parsing
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.MAX_ZIP_THREADS=Shared.threads();
		
		{//Parse the arguments
			final Parser parser=parse(args);
			overwrite=parser.overwrite;
			append=parser.append;
			
			in1=parser.in1;

			out1=parser.out1;
		}

		validateParams();
		fixExtensions(); //Add or remove .gz or .bz2 as needed
		checkFileExistence(); //Ensure files can be read and written
		checkStatics(); //Adjust file-related static fields as needed for this program

		ffout1=FileFormat.testOutput(out1, FileFormat.TXT, null, true, overwrite, append, false);
		ffoutInvalid=FileFormat.testOutput(outInvalid, FileFormat.TXT, null, true, overwrite, append, false);
		ffin1=FileFormat.testInput(in1, FileFormat.TXT, null, true, true);
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
//		parser.out1="stdout";
		
		//Parse each argument
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			
			//Break arguments into their constituent parts, in the form of "a=b"
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			if(b!=null && b.equalsIgnoreCase("null")){b=null;}

			if(a.equals("invalid")){
				outInvalid=b;
			}else if(a.equals("lines")){
				maxLines=Long.parseLong(b);
				if(maxLines<0){maxLines=Long.MAX_VALUE;}
			}else if(a.equals("verbose")){
				verbose=Parse.parseBoolean(b);
				ByteFile1.verbose=verbose;
				ByteFile2.verbose=verbose;
				ReadWrite.verbose=verbose;
			}else if(a.equals("parse_flag_goes_here")){
				long fake_variable=Parse.parseKMG(b);
				//Set a variable here
			}else if(parser.parse(arg, a, b)){
				//do nothing
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}
		
		return parser;
	}
	
	/** Add or remove .gz or .bz2 as needed */
	private void fixExtensions(){
		in1=Tools.fixExtension(in1);
		if(in1==null){throw new RuntimeException("Error - at least one input file is required.");}
	}
	
	/** Ensure files can be read and written */
	private void checkFileExistence(){
		//Ensure output files can be written
		if(!Tools.testOutputFiles(overwrite, append, false, out1)){
			outstream.println((out1==null)+", "+out1);
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output file "+out1+"\n");
		}
		
		//Ensure input files can be read
		if(!Tools.testInputFiles(false, true, in1)){
			throw new RuntimeException("\nCan't read some input files.\n");  
		}
		
		//Ensure that no file was specified multiple times
		if(!Tools.testForDuplicateFiles(true, in1, out1)){
			throw new RuntimeException("\nSome file names were specified multiple times.\n");
		}
	}
	
	/** Adjust file-related static fields as needed for this program */
	private static void checkStatics(){
		//Adjust the number of threads for input file reading
		if(!ByteFile.FORCE_MODE_BF1 && !ByteFile.FORCE_MODE_BF2 && Shared.threads()>2){
			ByteFile.FORCE_MODE_BF2=true;
		}
		
//		if(!ByteFile.FORCE_MODE_BF2){
//			ByteFile.FORCE_MODE_BF2=false;
//			ByteFile.FORCE_MODE_BF1=true;
//		}
	}
	
	/** Ensure parameter ranges are within bounds and required parameters are set */
	private boolean validateParams(){
//		assert(minfoo>0 && minfoo<=maxfoo) : minfoo+", "+maxfoo;
//		assert(false) : "TODO";
		return true;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Create streams and process all data */
	void process(Timer t){
		
		ByteFile bf=ByteFile.makeByteFile(ffin1);
		ByteStreamWriter bsw=makeBSW(ffout1);
		ByteStreamWriter bswInvalid=makeBSW(ffoutInvalid);
		
//		assert(false) : "Header goes here.";
		if(bsw!=null){
//			assert(false) : "Header goes here.";
		}
		
		processInner(bf, bsw, bswInvalid);
		
		errorState|=bf.close();
		if(bsw!=null){errorState|=bsw.poisonAndWait();}
		if(bswInvalid!=null){errorState|=bswInvalid.poisonAndWait();}
		
		t.stop();
		
//		outstream.println(Tools.timeLinesBytesProcessed(t, linesProcessed, bytesProcessed, 8));
//		
//		outstream.println();
//		outstream.println("Valid Lines:       \t"+linesOut);
//		outstream.println("Invalid Lines:     \t"+(linesProcessed-linesOut));
		
		//Throw an exception of there was an error in a thread
		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	private void processInner(ByteFile bf, ByteStreamWriter bsw, ByteStreamWriter bswInvalid){
		byte[] line=bf.nextLine();
		LongList sizes=new LongList(100000000);

		ArrayList<Pair> list=new ArrayList<Pair>(10000000);
		while(line!=null){
			linesProcessed++;
			bytesProcessed+=(line.length+1);
			processLine(line, list, sizes);
			if(linesProcessed>=maxLines) {break;}
			line=bf.nextLine();
		}
		sizes.sort();
		
		long total=sizes.sumLong();
		long mode=sizes.mode();
		long mean=sizes.sumLong()/sizes.size;
		long median=sizes.get((int)(sizes.size*0.5));
		
		Shared.sort(list);
		final long tebi=1024L*1024L*1024L*1024L;
		final long tera=1000L*1000L*1000L*1000L;
		final int[] idxArray=new int[] {10, 20, 30, 40, 50, 60, 70, 80, 90, 95};
		final int[] pairArray=idxArray.clone();
		for(int i=0; i<idxArray.length; i++) {
			pairArray[i]=(int)(idxArray[i]*.01*list.size());
		}
		long tsum=0;
		for(int i=0, nextIdx=0, nextPair=pairArray[0]; i<list.size() && nextIdx<idxArray.length; i++) {
			Pair p=list.get(i);
			tsum+=p.size;
			if(i>=nextPair) {
				String s=idxArray[nextIdx]+" percent of files have not been accessed since: "+
						timeString(p.time*1000)+" ("+(tsum/tebi)+" tebibytes)";
				System.out.println(s);
				nextIdx++;
				if(nextIdx<idxArray.length) {nextPair=pairArray[nextIdx];}
			}
		}

		System.out.println("total size: \t"+(total/tera)+" TB \t("+total+")"+"\t"+"("+((total/tebi))+" tebibytes)");
		System.out.println("mean size:  \t"+mean+" bytes");
		System.out.println("P50 size:   \t"+median+" bytes");
		System.out.println("P80 size:   \t"+sizes.get((int)(sizes.size*0.8))+" bytes");
		System.out.println("P90 size:   \t"+sizes.get((int)(sizes.size*0.9))+" bytes");
		System.out.println("P95 size:   \t"+sizes.get((int)(sizes.size*0.95))+" bytes");
	}

	static final byte delimiter=(byte)'|';
	
	void processLine(byte[] line, ArrayList<Pair> list, LongList sizes) {
		int a=0, b=0;

		while(b<line.length && line[b]!=delimiter){b++;}
		b++;
		a=b;

		while(b<line.length && line[b]!=delimiter){b++;}
		b++;
		a=b;

		while(b<line.length && line[b]!=delimiter){b++;}
		b++;
		a=b;

		while(b<line.length && line[b]!=delimiter){b++;}
		assert(b>a) : "Missing term : '"+new String(line)+"'";
		long size=Parse.parseLong(line, a, b);
		b++;
		a=b;

		while(b<line.length && line[b]!=delimiter){b++;}
		b++;
		a=b;

		while(b<line.length && line[b]!=delimiter){b++;}
		b++;
		a=b;

		while(b<line.length && line[b]!=delimiter){b++;}
		if(line[a]!='F') {return;}
		b++;
		a=b;

		while(b<line.length && line[b]!=delimiter){b++;}
		b++;
		a=b;

		while(b<line.length && line[b]!=delimiter){b++;}
		b++;
		a=b;

		while(b<line.length && line[b]!=delimiter){b++;}
		b++;
		a=b;

		while(b<line.length && line[b]!=delimiter){b++;}
		b++;
		a=b;

		while(b<line.length && line[b]!=delimiter){b++;}
		long time=Parse.parseLong(line, a, b);
		b++;
		a=b;
		
		linesOut++;
		sizes.add(size);
		list.add(new Pair(size, time));
	}
	
//	10 percent of files have not been accessed since: 2017-02-14 09:58:18 -0800 PST (130 TB)
//	20 percent of files have not been accessed since: 2017-05-26 10:52:24 -0700 PDT (204 TB)
//	30 percent of files have not been accessed since: 2019-03-13 14:54:27 -0700 PDT (381 TB)
//	40 percent of files have not been accessed since: 2021-01-05 03:06:55 -0800 PST (465 TB)
//	50 percent of files have not been accessed since: 2021-04-26 06:02:52 -0700 PDT (484 TB)
//	60 percent of files have not been accessed since: 2022-04-17 23:01:48 -0700 PDT (537 TB)
//	70 percent of files have not been accessed since: 2022-06-19 23:15:29 -0700 PDT (598 TB)
//	80 percent of files have not been accessed since: 2023-01-30 00:35:47 -0800 PST (772 TB)
//	90 percent of files have not been accessed since: 2023-04-25 16:49:28 -0700 PDT (1129 TB)
//	95 percent of files have not been accessed since: 2023-05-01 17:26:50 -0700 PDT (1182 TB)
//	Total Size: 1411 TB
//	P50 (Median) File Size: 10762 bytes
//	P80 File Size: 162376 bytes
//	P90 File Size: 1187502 bytes
//	P95 File Size: 10694780 bytes
	
	private static ByteStreamWriter makeBSW(FileFormat ff){
		if(ff==null){return null;}
		ByteStreamWriter bsw=new ByteStreamWriter(ff);
		bsw.start();
		return bsw;
	}
	
	private static class Pair implements Comparable<Pair> {
		
		Pair(long size_, long time_){
			size=size_;
			time=time_;
		}
		
		@Override
		public int compareTo(Pair b) {
			return time>b.time ? 1 : time<b.time ? -1 : 0;
		}
		
		final long size, time;
		
	}
	
	static String timeString(long time){
		SimpleDateFormat sdf=new SimpleDateFormat("yyyy-MM-dd HH:mm:ss");
		sdf.setTimeZone(TimeZone.getTimeZone("PST"));
//		sdf.setTimeZone(TimeZone.getDefault());
		return sdf.format(new Date(time));
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Primary input file path */
	private String in1=null;

	/** Primary output file path */
	private String out1=null;

	/** Junk output file path */
	private String outInvalid=null;
	
	/*--------------------------------------------------------------*/
	
	private long linesProcessed=0;
	private long linesOut=0;
	private long bytesProcessed=0;
	private long bytesOut=0;
	
	private long maxLines=Long.MAX_VALUE;
	
	/*--------------------------------------------------------------*/
	/*----------------         Final Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Input File */
	private final FileFormat ffin1;
	/** Output File */
	private final FileFormat ffout1;
	/** Optional Output File for Junk */
	private final FileFormat ffoutInvalid;
	
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
