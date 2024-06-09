package ml;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Locale;

import dna.AminoAcid;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import shared.Parse;
import shared.Parser;
import shared.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.ConcurrentReadInputStream;
import stream.Read;
import structures.ByteBuilder;
import structures.ListNum;
import tracker.EntropyTracker;

/**
 * @author Brian Bushnell
 * @date Oct 6, 2014
 *
 */
public class SequenceToVector {

	public static void main(String[] args){
		//Start a timer immediately upon code entrance.
		Timer t=new Timer();
		
		//Create an instance of this class
		SequenceToVector x=new SequenceToVector(args);
		
		//Run the object
		x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	public SequenceToVector(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		Parser parser=new Parser();
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			if(b!=null && b.equalsIgnoreCase("null")){b=null;}

			if(a.equals("parse_flag_goes_here")){
				//Set a variable here
			}else if(a.equals("result")){
				result0=Float.parseFloat(b);
			}else if(a.equals("width")){
				width=Integer.parseInt(b);
			}else if(a.equals("rcomp")){
				rcomp=Parse.parseBoolean(b);
			}else if(a.equals("parse")){
				parseHeader=Parse.parseBoolean(b);
			}else if(parser.parse(arg, a, b)){
				//do nothing
			}else{
				//				throw new RuntimeException("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
				outstream.println("Unknown parameter "+args[i]);
			}
		}
		
		{//Process parser fields
			Parser.processQuality();
			
			maxReads=parser.maxReads;
			in1=parser.in1;
			out1=parser.out1;
		}
		
		ffout1=FileFormat.testOutput(out1, FileFormat.TXT, null, true, true, false, false);
		ffin1=FileFormat.testInput(in1, FileFormat.FASTQ, null, true, true);
		assert(result0!=-1 || parseHeader);
	}
	
	void process(Timer t){
		
		final ConcurrentReadInputStream cris;
		{
			cris=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ffin1, null);
			cris.start();
		}
		boolean paired=cris.paired();
		
		ByteStreamWriter bsw=new ByteStreamWriter(ffout1);
		bsw.start();
		
		long readsProcessed=0, basesProcessed=0;
		{
			
			ListNum<Read> ln=cris.nextList();
			ArrayList<Read> reads=(ln!=null ? ln.list : null);
			
			if(reads!=null && !reads.isEmpty()){
				Read r=reads.get(0);
				assert((ffin1==null || ffin1.samOrBam()) || (r.mate!=null)==cris.paired());
			}
			
			int inputs=width*4+4;
			bsw.print("#dims\t"+inputs+"\t1\n");
			final ByteBuilder bb=new ByteBuilder();
			while(ln!=null && reads!=null && reads.size()>0){//ln!=null prevents a compiler potential null access warning
				if(verbose){outstream.println("Fetched "+reads.size()+" reads.");}
				
				for(int idx=0; idx<reads.size(); idx++){
					final Read r1=reads.get(idx);
					final Read r2=r1.mate;
					readsProcessed+=r1.pairCount();
					basesProcessed+=r1.pairLength();

					float result=(parseHeader ? Parse.parseFloat(r1.id, "result=", '\t') : result0);
					toVector(r1, bb, bsw, width, result, rcomp);
					toVector(r2, bb, bsw, width, result, rcomp);
				}

				cris.returnList(ln);
				if(verbose){outstream.println("Returned a list.");}
				ln=cris.nextList();
				reads=(ln!=null ? ln.list : null);
			}
			if(ln!=null){
				cris.returnList(ln.id, ln.list==null || ln.list.isEmpty());
			}
		}
		errorState=ReadWrite.closeStreams(cris) | errorState;
		if(verbose){outstream.println("Finished reading data.");}

		errorState=bsw.poisonAndWait() | errorState;
		
		t.stop();
		outstream.println("Time:                         \t"+t);
		outstream.println("Reads Processed:    "+readsProcessed+" \t"+String.format(Locale.ROOT, "%.2fk reads/sec", (readsProcessed/(double)(t.elapsed))*1000000));
		assert(!errorState) : "An error was encountered.";
	}
	
	/*--------------------------------------------------------------*/
	
	public static void toVector(Read r, ByteBuilder bb, ByteStreamWriter bsw, int width, 
			float result, boolean rcomp) {
		if(r==null) {return;}
		toVector(r.bases, bb, width, result);
		if(bsw!=null) {bsw.println(bb);}
		bb.clear();
		if(!rcomp) {return;}
		toVector(r.reverseComplement(), bb, bsw, width, result, false);
	}
	
	public static ByteBuilder toVector(byte[] bases, ByteBuilder bb, int width, float result) {
		float len=bases.length/(float)(width+5);
		float gc=Tools.calcGC(bases);
		EntropyTracker[] eTrackers=localETrackers.get();
		final int tnum=Tools.min(bases.length, eTrackers.length-1);
		EntropyTracker eTracker=eTrackers[tnum];
		assert(eTracker!=null) : tnum+", "+eTrackers.length;
		float entropy=eTracker.averageEntropy(bases, true);
		float poly=Read.longestHomopolymer(bases);
		poly=poly/(poly+5);
		bb.append(len,4).tab().append(gc, 4).tab().append(entropy, 4).tab().append(poly, 4);
		for(int i=0; i<bases.length && i<width; i++) {
			byte b=bases[i];
			int n=AminoAcid.baseToNumber4[b];
			bb.tab().append(codes[n]);
		}
		for(int i=bases.length; i<width; i++) {//Pad with zeros
			bb.tab().append("0\t0\t0\t0");
		}
		if(result==(int)result) {
			bb.tab().append((int)result);
		}else {
			bb.tab().append(result, 4);
		}
		return bb;
	}

	public static float[] fillVector(byte[] bases, float[] vec) {
		assert(vec!=null);
		return fillVector(bases, vec, 0, bases.length-1);
	}
	public static float[] fillVector(byte[] bases, float[] vec, int from, int to) {
		assert(vec!=null);
		final int len=(to-from+1);
		final int width=(vec.length-4)/4;
		float flen=len/(float)(width+5);
		float gc=Tools.calcGC(bases, from, to);
		EntropyTracker[] eTrackers=localETrackers.get();
		EntropyTracker eTracker=eTrackers[Tools.min(len, eTrackers.length-1)];
		float entropy=eTracker.averageEntropy(bases, true, from, to);
		float poly=Read.longestHomopolymer(bases, from, to);
		poly=poly/(poly+5);
		Arrays.fill(vec, 0);
		vec[0]=flen;
		vec[1]=gc;
		vec[2]=entropy;
		vec[3]=poly;
//		System.err.println(Arrays.toString(vec));
		for(int i=from, j=4, lim=Tools.min(to+1, from+width, bases.length); i<lim; i++, j+=4) {
//		for(int i=0, j=4; i<bases.length && i<width; i++, j+=4) {
			byte b=bases[i];
			int n=AminoAcid.baseToNumber0[b];
			vec[j+n]=1;
		}
		assert(vec.length==width*4+4);
//		vec[width*4+4]=result;
		return vec;
	}
	
	public static float score(byte[] bases, float[] vec, CellNet net, boolean rcomp) {
		net.applyInput(vec);
		final float r, f=net.feedForwardFast();
		if(rcomp) {
			AminoAcid.reverseComplementBasesInPlace(bases);
			net.applyInput(vec);
			r=net.feedForwardFast();
			AminoAcid.reverseComplementBasesInPlace(bases);
		}else {r=f;}
		return Tools.max(r, f);
	}
	
	/*--------------------------------------------------------------*/
	
	private String in1=null;
	private String out1=null;
	
	private final FileFormat ffin1;
	private final FileFormat ffout1;

	private static final ThreadLocal<EntropyTracker[]> localETrackers=new ThreadLocal<EntropyTracker[]>(){
        @Override protected EntropyTracker[] initialValue() {
        	EntropyTracker[] array=new EntropyTracker[maxWindow+1];
        	for(int i=minWindow; i<array.length; i++) {
        		array[i]=new EntropyTracker(ke, i, false);
        	}
        	return array;
        }
    };
//	public EntropyTracker eTracker=new EntropyTracker(3, 20, false);
	
	static String[] codes=new String[] {"1\t0\t0\t0", "0\t1\t0\t0", "0\t0\t1\t0", "0\t0\t0\t1", "1\t0\t0\t0"};
	
	/*--------------------------------------------------------------*/

	private long maxReads=-1;
	private boolean errorState=false;

	private boolean rcomp=false;
	private boolean parseHeader=false;
	private int width=55;
	float result0=-1;

	public static final int minWindow=16;
	public static final int maxWindow=40;
	public static final int ke=3;
	
	
	/*--------------------------------------------------------------*/
	
	private java.io.PrintStream outstream=System.err;
	public static boolean verbose=false;
	
}
