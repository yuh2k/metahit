package bin;

import java.io.File;
import java.util.ArrayList;

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
import stream.Read;
import structures.IntHashMap;
import structures.IntLongHashMap;
import structures.ListNum;
import structures.LongList;

/**
 * Accepts multiple input files.
 * Reads them each sequentially, and outputs everything to a single output file.
 * Generically, it can be used to concatenate files while recompressing them
 * and avoiding the use of stdio.
 * @author Brian Bushnell
 * @date May 10, 2024
 *
 */
public class GradeBins {

	public static void main(String[] args){
		//Start a timer immediately upon code entrance.
		Timer t=new Timer();
		
		//Create an instance of this class
		GradeBins x=new GradeBins(args);
		
		//Run the object
		x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	public GradeBins(String[] args){
		
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
			}else if(a.equals("in")){
				Tools.getFileOrFiles(b, in, true, false, false, false);
			}else if(a.equals("size")){
				totalSize=Parse.parseKMG(b);
			}else if(a.equals("ref") || a.equals("contigs") || a.equals("assembly")){
				ref=b;
			}else if(a.equals("hist")){
				hist=b;
			}else if(b==null && new File(arg).isFile()){
				in.add(arg);
			}else if(b==null && new File(arg).isDirectory()){
				Tools.getFileOrFiles(arg, in, true, false, false, false);
			}else if(parser.parse(arg, a, b)){
				//do nothing
			}else{
				assert(false) : "Unknown parameter "+args[i];
				outstream.println("Unknown parameter "+args[i]);
			}
		}
		
		{//Process parser fields
			Parser.processQuality();
			
			maxReads=parser.maxReads;
			out1=parser.out1;
		}
	}
	
	void process(Timer t){
		
		BinObject.grading=true;
		IntLongHashMap sizeMap=makeSizeMap(ref);

		long cleanBins=0, contamBins=0;
		long cleanContigs=0, contamContigs=0;
		long cleanSize=0, contamSize=0;
		long partialCleanSize=0, partialContamSize=0;
		IntHashMap tidBins=new IntHashMap();
		for(String s : in) {
			Cluster c=loadCluster(s);
			c.calcContam(sizeMap);
			
			sizes.add(c.size);
			bins.add(new BinStats(c));
			if(c.taxid>0) {
				tidBins.increment(c.taxid);
			}
			long contam=Math.round(c.contam*c.size);
			contamScore+=contam;
			compltScore+=Math.round(c.completeness*(c.size-contam));
			if(contam<1) {
				cleanBins++;
				cleanSize+=c.size;
				cleanContigs+=c.numContigs();
			}else {
				contamBins++;
				contamSize+=c.size;
				contamContigs+=c.numContigs();
				partialCleanSize+=(c.size-contam);
				partialContamSize+=contam;
			}
		}
		
		if(verbose){outstream.println("Finished.");}

		outstream.println();
		outstream.println(QuickBin.formatString("Clean Bins", 29, cleanBins, contamBins));
		outstream.println(QuickBin.formatString("Dirty Bins", 29, contamBins, cleanBins));
		outstream.println(QuickBin.formatString("Clean Bin Bases", 29, cleanSize, contamSize));
		outstream.println(QuickBin.formatString("Dirty Bin Bases", 29, contamSize, cleanSize));
		outstream.println(QuickBin.formatString("Tainted Bases", 29, 
				partialCleanSize, cleanSize+contamSize-partialCleanSize));
		outstream.println(QuickBin.formatString("Contam Bases", 29, 
				partialContamSize, cleanSize+contamSize-partialContamSize));
		outstream.println();

		outstream.println("Sequence Recovery:           \t"+
				String.format("%.3f", (cleanSize+contamSize)*100.0/totalSize));
		outstream.println("Contig Recovery:             \t"+
				String.format("%.3f", (cleanContigs+contamContigs)*100.0/totalContigs));
		outstream.println("Genomes Represented:         \t"+
				String.format("%.3f", (tidBins.size())*100.0/taxIDsIn));
		outstream.println("Completeness Score:          \t"+
				String.format("%.3f", 100*compltScore/totalSize));
		outstream.println("Contamination Score:         \t"+
				String.format("%.4f", 100*contamScore/totalSize));
		outstream.println();
		
		printL90(sizes, totalSize);
		
		if(hist!=null) {
			ChartMaker.makeChartFromBinStats(hist, bins);
		}
		
		t.stop();
		outstream.println();
		outstream.println(Tools.timeReadsBasesProcessed(t, readsProcessed, basesProcessed, 8));
	}
	
	Cluster loadCluster(String fname) {
		FileFormat ffin=FileFormat.testInput(fname, FileFormat.FASTA, null, true, true);
		
		final ConcurrentReadInputStream cris;
		{
			cris=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ffin, null);
			cris.start();
		}
		Cluster c=new Cluster(0);		
		{
			ListNum<Read> ln=cris.nextList();
			ArrayList<Read> reads=(ln!=null ? ln.list : null);

			//ln!=null prevents a compiler potential null access warning
			while(ln!=null && reads!=null && reads.size()>0){
				if(verbose){outstream.println("Fetched "+reads.size()+" reads.");}
				
				for(int idx=0; idx<reads.size(); idx++){
					final Read r1=reads.get(idx);
					readsProcessed+=r1.pairCount();
					basesProcessed+=r1.pairLength();
					
					//  *********  Process reads here  *********
					Contig a=new Contig(r1.name(), r1.bases, (int)r1.numericID);
					c.counts=new int[0];
					int tid=DataLoader.parseTaxID(a.name);
					a.taxid=a.labelTaxid=tid;
					c.add(a);
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
		ReadWrite.closeStream(cris);
		return c;
	}
	
	static void printL90(LongList list, long basesLoaded) {
		long c99=(long)(0.99f*basesLoaded);
		long c95=(long)(0.95f*basesLoaded);
		long c90=(long)(0.90f*basesLoaded);
		long c80=(long)(0.80f*basesLoaded);
		long c75=(long)(0.75f*basesLoaded);
		long c50=(long)(0.50f*basesLoaded);
		long c40=(long)(0.40f*basesLoaded);
		long c30=(long)(0.30f*basesLoaded);
		long c25=(long)(0.25f*basesLoaded);
		long c20=(long)(0.20f*basesLoaded);
		long c10=(long)(0.10f*basesLoaded);
		long c05=(long)(0.05f*basesLoaded);
		long c01=(long)(0.01f*basesLoaded);
		
		list.sort();
		list.reverse();
		long prev=0, sum2=0;
		for(int i=0; i<list.size(); i++) {
			long size=list.get(i);
			prev=sum2;
			sum2+=size;
			int num=i+1;

			if(sum2>=c01 && prev<c01) {System.err.println("L01: "+size+"\t"+"N01: "+num);}
			if(sum2>=c05 && prev<c05) {System.err.println("L05: "+size+"\t"+"N05: "+num);}
			if(sum2>=c10 && prev<c10) {System.err.println("L10: "+size+"\t"+"N10: "+num);}
			if(sum2>=c20 && prev<c20) {System.err.println("L20: "+size+"\t"+"N20: "+num);}
//			if(sum2>=c25 && prev<c25) {System.err.println("L25: "+size+"\t"+"N25: "+num);}
			if(sum2>=c30 && prev<c30) {System.err.println("L30: "+size+"\t"+"N30: "+num);}
			if(sum2>=c40 && prev<c40) {System.err.println("L40: "+size+"\t"+"N40: "+num);}
			if(sum2>=c50 && prev<c50) {System.err.println("L50: "+size+"\t"+"N50: "+num);}
//			if(sum2>=c75 && prev<c75) {System.err.println("L75: "+size+"\t"+"N75: "+num);}
			if(sum2>=c80 && prev<c80) {System.err.println("L80: "+size+"\t"+"N80: "+num);}
			if(sum2>=c90 && prev<c90) {System.err.println("L90: "+size+"\t"+"N90: "+num);}
			if(sum2>=c95 && prev<c95) {System.err.println("L95: "+size+"\t"+"N95: "+num);}
			if(sum2>=c99 && prev<c99) {System.err.println("L99: "+size+"\t"+"N99: "+num);}
		}
	}
	
	IntLongHashMap makeSizeMap(String fname) {
		FileFormat ffin=FileFormat.testInput(fname, FileFormat.FASTA, null, true, true);
		
		final ConcurrentReadInputStream cris;
		{
			cris=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ffin, null);
			cris.start();
		}
		IntLongHashMap map=new IntLongHashMap();
		long sizeSum=0, contigSum=0;
		{
			ListNum<Read> ln=cris.nextList();
			ArrayList<Read> reads=(ln!=null ? ln.list : null);

			//ln!=null prevents a compiler potential null access warning
			while(ln!=null && reads!=null && reads.size()>0){
				if(verbose){outstream.println("Fetched "+reads.size()+" reads.");}
				
				for(int idx=0; idx<reads.size(); idx++){
					final Read r=reads.get(idx);
					readsProcessed++;
					basesProcessed+=r.length();
					sizeSum+=r.length();
					contigSum++;
					
					//  *********  Process reads here  *********
					int tid=DataLoader.parseTaxID(r.id);
					long ret=map.increment(tid, r.length());
					if(ret==r.length() && tid>0) {taxIDsIn++;}
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
		if(totalSize==0) {totalSize=sizeSum;}
		if(totalContigs==0) {totalContigs=contigSum;}
		ReadWrite.closeStream(cris);
		return map;
	}
	
	/*--------------------------------------------------------------*/
	
	/*--------------------------------------------------------------*/
	
	private ArrayList<String> in=new ArrayList<String>();
	private String ref=null;
	private String out1=null;
	private String hist=null;
	private LongList sizes=new LongList();
	private ArrayList<BinStats> bins=new ArrayList<BinStats>();
	private double contamScore=0;
	private double compltScore=0;
	
	/*--------------------------------------------------------------*/

	private long maxReads=-1;
	private long readsProcessed=0, basesProcessed=0;
	private long totalSize=0, totalContigs=0;
	private long taxIDsIn=0;
	
	/*--------------------------------------------------------------*/
	
	private java.io.PrintStream outstream=System.err;
	public static boolean verbose=false;
	
}
