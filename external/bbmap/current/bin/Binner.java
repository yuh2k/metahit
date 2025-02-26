package bin;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.concurrent.locks.ReadWriteLock;

import bin.SpectraCounter.LoadThread;
import shared.Parse;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import structures.IntHashMap;
import structures.IntHashSet;
import structures.IntLongHashMap;
import tax.TaxTree;
import template.Accumulator;
import template.ThreadWaiter;

public class Binner extends BinObject implements Accumulator<Binner.CompareThread> {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	Binner(PrintStream outstream_){outstream=outstream_;}
	
	boolean parse(String arg, String a, String b) {
	
		if(a.equalsIgnoreCase("productMult")){
			productMult=Float.parseFloat(b);
		}
		
		else if(a.equalsIgnoreCase("residueRange")){
			baseRange=Integer.parseInt(b);
		}else if(a.equalsIgnoreCase("residueRange")){
			residueRange=Integer.parseInt(b);
		}

		else if(a.equalsIgnoreCase("maxDif1") || a.equalsIgnoreCase("maxKmerDif1")){
			maxKmerDif1=Float.parseFloat(b);
		}else if(a.equalsIgnoreCase("maxRatio1") || a.equalsIgnoreCase("maxDepthRatio1")){
			maxDepthRatio1=Float.parseFloat(b);
		}else if(a.equalsIgnoreCase("maxGCDif1")){
			maxGCDif1=Float.parseFloat(b);
		}else if(a.equalsIgnoreCase("maxCovariance1")){
			maxCovariance1=Float.parseFloat(b);
		}
		
		else if(a.equalsIgnoreCase("maxDif2") || a.equalsIgnoreCase("maxKmerDif2")){
			maxKmerDif2=Float.parseFloat(b);
		}else if(a.equalsIgnoreCase("maxRatio2") || a.equalsIgnoreCase("maxDepthRatio2")){
			maxDepthRatio2=Float.parseFloat(b);
		}else if(a.equalsIgnoreCase("maxGCDif2")){
			maxGCDif2=Float.parseFloat(b);
		}else if(a.equalsIgnoreCase("maxCovariance2")){
			maxCovariance2=Float.parseFloat(b);
		}
		
//		else if(a.equalsIgnoreCase("minSizeToCluster")){
//			minSizeToCluster=Parse.parseIntKMG(b);
//		}else if(a.equalsIgnoreCase("minSizeToRefine")){
//			minSizeToRefine=Parse.parseIntKMG(b);
//		}
		else if(a.equalsIgnoreCase("minseedsize") || a.equals("minsizeseed") || a.equals("minseed")){
			minSizeToCompare=minSizeToMerge=Parse.parseIntKMG(b);
		}else if(a.equalsIgnoreCase("minSizeToCompare")){
			minSizeToCompare=Parse.parseIntKMG(b);
		}else if(a.equalsIgnoreCase("minSizeToMerge")){
			minSizeToMerge=Parse.parseIntKMG(b);
//		}else if(a.equalsIgnoreCase("minSizeToAdd")){
//			minSizeToAdd=Parse.parseIntKMG(b);
		}else if(a.equalsIgnoreCase("minSizeResidue") || a.equalsIgnoreCase("minResidue")){
			minSizeResidue=Parse.parseIntKMG(b);
		}
		
		else if(a.equalsIgnoreCase("residueStringency")){
			residueStringency=Float.parseFloat(b);
		}else if(a.equalsIgnoreCase("smallthresh")){
			smallThresh=Parse.parseIntKMG(b);
		}else if(a.equalsIgnoreCase("smallmult")){
			smallMult=Float.parseFloat(b);
		}else if(a.equalsIgnoreCase("bigthresh")){
			bigThresh=Parse.parseIntKMG(b);
		}else if(a.equalsIgnoreCase("bigmult")){
			bigMult=Float.parseFloat(b);
		}
		
		else if(a.equalsIgnoreCase("maxEdges")){
			maxEdges=Integer.parseInt(b);
		}else if(a.equalsIgnoreCase("minEdgeWeight")){
			minEdgeWeight=Integer.parseInt(b);
		}else if(a.equalsIgnoreCase("minEdgeWeight0")){
			minEdgeWeight0=Integer.parseInt(b);
		}else if(a.equalsIgnoreCase("minEdgeRatio")){
			minEdgeRatio=Float.parseFloat(b);
		}
		
//		else if(a.equalsIgnoreCase("mtcompare") || a.equalsIgnoreCase("comparemt")){
//			multiThreadedCompare=Parse.parseBoolean(b);
//		}
		
		else {return false;}
		
		return true;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------             Graph            ----------------*/
	/*--------------------------------------------------------------*/
	
	public int followEdges(ArrayList<Contig> contigs, float stringency){
		return followEdges(contigs, contigs, stringency);
	}
	
	public int followEdges(ArrayList<Contig> contigs, ArrayList<? extends Bin> input, float stringency){
//		outstream.print("Following Edges:  \t");
		phaseTimer.start();
		int merges=0;
		for(Bin a : input) {
			if(a.pairMap!=null) {
				int dest=followEdges(a, contigs, stringency);
				if(dest>=0) {
					merges++;
					a.dest=dest;
				}
			}
		}
		
		if(merges>0) {merges=mergeWithDest(contigs, input);}
		
//		phaseTimer.stopAndPrint();
		if(merges>0) {outstream.println("Merged "+merges+"/"+input.size()+" bins.");}
		return merges;
	}
	
	private int mergeWithDest(ArrayList<Contig> contigs, ArrayList<? extends Bin> input){
		int x=0, y=0, z=0;
		
		for(Bin a : input) {
			if(!a.isCluster() && a.cluster()!=null) {a.dest=-1;continue;}
			if(a.isCluster() && a.numContigs()==0) {a.dest=-1;continue;}
			assert(a.isValid());
			if(a.dest<0) {
				//ignore
				x++;
			}else {
				y++;
				Bin b=contigs.get(a.dest);
				if(b.cluster()!=null) {b=b.cluster();}
				if(!b.isEmpty() && !a.sameCluster(b)) {
					if(a.labelTaxid>=0 && b.labelTaxid>=0) {
						if(a.labelTaxid==b.labelTaxid) {goodMergesFollow++;}
						else {badMergesFollow++;}
					}
					assert(b.isValid());
					if(a.isCluster()) {
						try {
							((Cluster)a).add(b);
						} catch (Throwable e) {
							System.err.println(a.numContigs()+", "+((Cluster)a).contigs);
							// TODO Auto-generated catch block
							e.printStackTrace();
							throw new RuntimeException(e);
						}
					}else if(b.isCluster()) {
						((Cluster)b).add(a);
					}else if(a.cluster()!=null){
						assert(false);
						a.cluster().add(b);
					}else {
						Cluster c=a.toCluster();
						c.add(b);
					}
					z++;
				}
			}
			a.dest=-1;
		}
		return z;
	}
	
	private static final int countClustered(ArrayList<? extends Bin> list) {
		int clustered=0;
		for(Bin b : list) {
			clustered+=(b.cluster()!=null ? b.numContigs() : 0);
		}
		return clustered;
	}
	
	public static final ArrayList<Bin> toBinList(ArrayList<? extends Bin> list){
		ArrayList<Bin> bins=new ArrayList<Bin>();
		IntHashSet clusterSet=new IntHashSet(255);
		for(Bin a : list) {
			Cluster c=a.cluster();
			if(c==null) {//Contig
				bins.add(a);
			}else if(!clusterSet.contains(c.id())){
				bins.add(c);
				clusterSet.add(c.id());
			}
		}
		return bins;
	}
	
	private int followEdges(Bin a, ArrayList<Contig> contigs, float stringency) {
		ArrayList<KeyValue> edges=KeyValue.toList(a.pairMap);//TODO: Slow
		float bestScore=0;
		Bin target=null;
		assert(a.isCluster() || a.cluster()==null);
		
		int max=(a.isCluster() ? maxEdges : maxEdges+Tools.min(2, maxEdges)*Tools.min(8, a.numContigs()-1));
		max=Tools.min(max, edges.size());
		int minWeight=(int)Math.ceil(minEdgeRatio*edges.get(0).value);
		minWeight=Tools.max(minWeight, minEdgeWeight);
		
//		{System.err.print("^");}
		for(int i=0; i<max; i++) {
//			{System.err.print("_");}
			KeyValue kv=edges.get(i);
			if(kv.value<minWeight) {break;}
			Contig c=contigs.get(kv.key);
//			Cluster cc=c.cluster();
			Bin b=(c.cluster()!=null ? c.cluster() : c);
			if(a==b || target==b) {
//				{System.err.print("["+a.id()+","+b.id()+"]");}
				continue;
			}
//			{System.err.print("("+kv.value+","+c.countEdgesTo(a)+")");}
			int min=Tools.min((int)c.countEdgesTo(a), kv.value);
//			if(min>1) {System.err.print(".");}
			if(min>=minWeight && a.canMergeWith(b, stringency)) {
//				System.err.print("@");
				float f=a.similarityTo(b);
				if(f>bestScore) {
					target=b;
					bestScore=f;
				}
			}
		}
		return target==null ? -1 : target.id();
	}
	
	public Cluster findLinkedCluster(Bin a, ArrayList<Contig> contigs) {
		if(a.pairMap==null) {return null;}
		assert(a.isCluster() || a.cluster()==null);
		int bestValue=0;
		Cluster bestCluster=null;
		int[] keys=a.pairMap.keys(), values=a.pairMap.values();
//		final int max=Tools.max(values);
		for(int i=0; i<keys.length; i++) {
			int key=keys[i];
			if(key!=a.pairMap.invalid()) {
				int v=values[i];
				Contig a2=contigs.get(key);
				assert(key==a2.id());
				Cluster c=contigs.get(key).cluster();
				if(c!=null && v>bestValue && c!=a) {
					assert(c.contigSet.contains(key));
					if(!a.isCluster() || a.id()<=c.id) {
						if(c.countEdgesTo(a)>0) {
							bestValue=v;
							bestCluster=c;
						}
					}
				}
			}
		}
		return bestCluster;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Binning           ----------------*/
	/*--------------------------------------------------------------*/
	
	public BinMap makeBinMap(ArrayList<Contig> contigList, ArrayList<? extends Bin> input) {
		outstream.print("Making BinMap:    \t");
		phaseTimer.start();
		
		if(input==null) {
			input=contigList;
		}else {
			Collections.sort(input);
		}
		BinMap map=new BinMap(contigList);
		Key key=new Key();
		float stringency=1;
		long contigsAdded=0;
		long clustersCreated=0;
		
		for(int i=1; i<input.size(); i++) {
			assert(input.get(i).size()<=input.get(i-1).size());
		}
		
		for(Bin a : input) {
			int initialContigs=a.numContigs();
			Cluster c=map.addOrMerge(a, minSizeToCompare*8, minSizeToMerge*4, minSizeToCompare,	
					maxKmerDif1, maxDepthRatio1, maxGCDif1, maxCovariance1, stringency, 
					TaxTree.SPECIES, true, true, key, 0);
			if(c==null) {
				//residual.add(b);//Automatic
			}else {
				contigsAdded+=initialContigs;
				if(c.numContigs()==1) {
					clustersCreated++;
				}else if(c.labelTaxid>0 && a.labelTaxid>0) {
					goodMergesCreate+=(c.labelTaxid==a.labelTaxid) ? 1 : 0;
					badMergesCreate+=(c.labelTaxid==a.labelTaxid) ? 0 : 1;
				}
			}
		}
		
		phaseTimer.stopAndPrint();
		outstream.println("Made "+map.map.size()+" lists containing "+clustersCreated+
				" clusters and "+contigsAdded+" contigs from "+contigList.size()+" elements.");
		return map;
	}
	
	public int refineBinMapPass(BinMap map, float stringency, 
			int taxlevel, boolean allowNoTaxID, boolean allowHalfTaxID, int range, int minSize) {
//		System.err.println("Merging clusters pass.");
		
		float maxKmerDif=maxKmerDif2*stringency;
		float maxDepthRatio=1+((maxDepthRatio2-1)*stringency);
		float maxGCDif=maxGCDif2*stringency;
		float maxProduct=maxKmerDif*maxDepthRatio*Binner.productMult;
		float maxCovariance=maxCovariance2*stringency;
		
		ArrayList<Cluster> clusters=map.toList(false);
		Collections.sort(clusters);
		Key key=new Key();
		
		for(int i=1; i<clusters.size(); i++) {
			assert(clusters.get(i).size()<=clusters.get(i-1).size());
		}
		
//		System.err.println("maxKmerDif="+maxKmerDif+", maxDepthRatio="+maxDepthRatio+", maxGCDif="+maxGCDif
//				+", maxProduct="+maxProduct+", allowNoTaxID="+allowNoTaxID+", allowHalfTaxID="+allowHalfTaxID
//				+", range="+range);
		
		int merged=0;
		
		launchThreads(clusters, map, REFINE_MODE, 
				maxKmerDif, maxDepthRatio, maxGCDif, maxProduct, maxCovariance, 
				taxlevel, allowNoTaxID, allowHalfTaxID, range, minSize);
		for(Cluster c : clusters) {
			synchronized(c) {
				if(c.dest>=0) {
					merged++;
				}
			}
		}
//		assert(false);
		assert(map.isValid());
		if(merged<1) {return 0;}
		
		for(int i=clusters.size()-1; i>=0; i--) {
			Cluster a=clusters.get(i);
			synchronized(a) {
				final int dest=a.dest;
				if(dest>=0 && dest!=a.id()) {
					Cluster b=map.contigList.get(dest).cluster;
					synchronized(b) {
						if(b!=a) {
							assert(!b.contigSet.contains(a.id));
							assert(!a.contigSet.contains(b.id));
							assert(a.id()!=dest && a.id()!=b.id()) : a.id+", "+dest+", "+b.id+", "+(a.id()==dest)+", "+(b.id()==dest)+", "+(a.id()!=b.id());

							if(a.labelTaxid>0 && b.labelTaxid>0) {
								goodMergesRefine+=(a.labelTaxid==b.labelTaxid) ? 1 : 0;
								badMergesRefine+=(a.labelTaxid==b.labelTaxid) ? 0 : 1;
							}
							b.add(a);
							clusters.set(i, null);
						}
					}
				}
			}
		}

		map.clear(false);
		assert(map.isValid());
		Tools.condense(clusters);
		Collections.sort(clusters);
		for(Cluster c : clusters) {
			map.add(c, key);
		}
		assert(map.isValid());
		return merged;
	}
	
	public int processResidue(BinMap map, float stringency, 
			int taxlevel, boolean allowNoTaxID, boolean allowHalfTaxID, int range) {
		assert(map.isValid());
		System.err.println("Processing "+map.residual.size()+" residual contigs.");
		Timer t=new Timer(outstream, true);
		if(map.residual.isEmpty()) {return 0;}
		
		float maxKmerDif=maxKmerDif2*stringency;
		float maxDepthRatio=1+((maxDepthRatio2-1)*stringency);
		float maxGCDif=maxGCDif2*stringency;
		float maxProduct=maxKmerDif*maxDepthRatio*Binner.productMult;
		float maxCovariance=maxCovariance2*stringency;

//		System.err.println("maxKmerDif="+maxKmerDif+", maxDepthRatio="+maxDepthRatio+", maxGCDif="+maxGCDif
//				+", maxProduct="+maxProduct+", allowNoTaxID="+allowNoTaxID+", allowHalfTaxID="+allowHalfTaxID
//				+", range="+range);
		int merged=0;
		map.residual=toBinList(map.residual);
		assert(map.isValid());
		
		int minSize=Tools.max(minSizeToMerge, minBasesPerCluster/5);
		launchThreads(map.residual, map, RESIDUE_MODE, maxKmerDif, maxDepthRatio, maxGCDif, maxProduct, maxCovariance, 
				taxlevel, allowNoTaxID, allowHalfTaxID, range, minSize);
		for(Bin c : map.residual) {
			synchronized(c) {
				if(c.dest>=0) {merged++;}
			}
		}
		t.stopAndStart("Found "+merged+" merge targets.");
		
		if(merged<1) {return 0;}
		for(int i=0; i<map.residual.size(); i++) {
			Bin a=map.residual.get(i);
			if(a.dest>0) {
				Cluster b=map.contigList.get(a.dest).cluster;
				assert(a!=b) : "\n"+a.id()+"\n"+b.id()+"\n"+a.dest+"\n"+
						a.getClass()+"\n"+b.getClass()+"\n"+b.contigSet.contains(a.id());
				if(a.labelTaxid>0 && b.labelTaxid>0) {
					goodMergesResidue+=(a.labelTaxid==b.labelTaxid) ? 1 : 0;
					badMergesResidue+=(a.labelTaxid==b.labelTaxid) ? 0 : 1;
				}
				b.add(a);
				map.residual.set(i, null);
			}
		}
		
		Tools.condenseStrict(map.residual);
		for(ArrayList<Cluster> list : map.map.values()) {
			Collections.sort(list);
		}
		assert(map.isValid());

		t.stop("Merged "+merged+" contigs into clusters.");
		return merged;
	}
	
	public Cluster findBestResidualCluster(Bin a, BinMap map, float maxKmerDif, 
			float maxDepthRatio, float maxGCDif, float maxProduct, float maxCovariance,
			Key key, int range, int minSize) {
		if(a==null || a.size()<minSizeResidue) {return null;}
		int minSize2=(int)Tools.max(minSizeToCompare, a.size(), minSize);
		float mult=sizeAdjustMult(a.size());
		Cluster b=map.findBestCluster(a, minSize2, maxKmerDif*mult, maxDepthRatio*mult, maxProduct*mult, 
				maxGCDif*mult, maxCovariance, -1, true, true, key, range);
		return b;
	}
		
	public int refineBinMap(BinMap map) {
		System.err.println("Merging clusters.");
		if(sketchClusters) {sketcher.sketch(map.toList(false), false);}
		else {
			for(Cluster c : map) {
				if(c.sketchedSize()>=2*c.size()) {c.clearTax();}
			}
		}
		phaseTimer.start();
		
		int removedThisPhase=0;
		int removedTotal=0;
		if(sketchContigs || sketchClusters) {
			removedThisPhase=refinePhase(map, "aa", 2.5f, TaxTree.SPECIES, false, false, baseRange, 8);
			removedTotal+=removedThisPhase;

			removedThisPhase=refinePhase(map, "bb", 1.5f, TaxTree.SPECIES, true, false, baseRange, 8);
			removedTotal+=removedThisPhase;

			removedThisPhase=refinePhase(map, "cc", 1.0f, TaxTree.GENUS, true, true, baseRange, 8);
			removedTotal+=removedThisPhase;
		}else {
			removedThisPhase=refinePhase(map, "a", 0.75f, -1, true, true, baseRange, 2);
			removedTotal+=removedThisPhase;
			removedThisPhase=refinePhase(map, "b", 0.8f, -1, true, true, baseRange+2, 3);
			removedTotal+=removedThisPhase;
		}

		removedThisPhase=refinePhase(map, "d", 0.9f, -1, true, true, baseRange, 1);
		removedTotal+=removedThisPhase;
		removedThisPhase=refinePhase(map, "e", 1.0f, -1, true, true, baseRange+2, 4);
		removedTotal+=removedThisPhase;
		
		phaseTimer.stop("Refinement merged "+removedTotal+" clusters. ");
		return removedTotal;
	}
	
	int refinePhase(BinMap map, String phase,
			float stringency, int taxLevel, boolean noTax, boolean halfTax, int range, int passes) {
		Timer t=new Timer(outstream, true);
		int removedThisPhase=0;
		for(int pass=1; pass<=passes; pass++) {
			int initial=map.countClusters();
			int removed=refineBinMapPass(map, stringency, taxLevel, noTax, halfTax, range, minSizeToMerge*pass);
			removedThisPhase+=removed;
			System.err.print("Refinement Pass "+pass+phase+": Merged "+removed+"/"+initial+" clusters. ");
			t.stopAndStart("\t");
			if(removed<2) {break;}
		}
		if(sketchClusters && removedThisPhase>0) {sketcher.sketch(map.toList(false), false);}
		return removedThisPhase;
	}
	
	public ArrayList<Bin> clusterByTaxid(ArrayList<? extends Bin> bins){
		outstream.print("Clustering by Taxid: \t");
		phaseTimer.start();
		Collections.sort(bins);
		HashMap<Integer, Bin> map=new HashMap<Integer, Bin>();
		int clustersMade=0;
		int contigsClustered=0;

		ArrayList<Bin> out=new ArrayList<Bin>();
		for(int i=0; i<bins.size(); i++) {
			Bin b=bins.get(i);
			if(b.taxid()>0) {
				Integer key=Integer.valueOf(b.taxid());
				Bin old=map.get(key);
				if(old==null) {
					map.put(key, b);
					b=null;
				}else if(old.getClass()==Cluster.class) {
					//Todo: Write "similar" function.
					((Cluster)old).add(b);
					contigsClustered++;
					b=null;
				}else {
					Cluster a=new Cluster(clustersMade);
					a.add(old);
					a.add(b);
					map.put(key, a);
					clustersMade++;
					contigsClustered++;
					b=null;
				}
			}
			if(b!=null) {out.add(b);}
		}
		
		out.addAll(map.values());
		Collections.sort(out);
		
		phaseTimer.stopAndPrint();
		outstream.println("Made "+clustersMade+" clusters containing "+contigsClustered+"/"+bins.size()+" elements.");
		return out;
	}
	
	/*--------------------------------------------------------------*/
	
	void setSamples(int samples) {
		if(samples<2) {//Single mode
			maxKmerDif2=0.005f; //.005 best for 1 depth; .01 best for 4 depths.
			maxDepthRatio2=1.2f;
			maxGCDif2=0.025f; //0.02f for 1 depth, 0.05 for 4 depths
			smallThresh=10000;
			smallMult=2.2f;
		}else if(samples<3){//Two mode
			maxKmerDif2=0.0065f;
			maxDepthRatio2=1.25f;
			maxGCDif2=0.03f;
			smallThresh=10500;
			smallMult=2.25f;
		}else if(samples<4){//Three mode
			maxKmerDif2=0.01f;
			maxDepthRatio2=1.3f;
			maxGCDif2=0.04f;
			smallThresh=11000;
			smallMult=2.3f;
		}else if(samples<5){//Four mode
			maxKmerDif2=0.02f;
			maxDepthRatio2=1.6f;
			maxGCDif2=0.05f;
			smallThresh=12000;
			smallMult=2.4f;
			
			maxKmerDif1=0.004f;
			maxDepthRatio1=1.1f;
			maxGCDif1=0.02f;
			maxCovariance2=0.0003f;
		}else if(samples<7){//5-6
			maxKmerDif2=0.022f;
			maxDepthRatio2=1.65f;
			maxGCDif2=0.055f;
			smallThresh=12500;
			smallMult=2.5f;
			
			maxKmerDif1=0.004f;
			maxDepthRatio1=1.1f;
			maxGCDif1=0.02f;
			maxCovariance1=0.0001f;
			maxCovariance2=0.0013f;
		}else if(samples<9){//7-8
			maxKmerDif2=0.022f;
			maxDepthRatio2=1.65f;
			maxGCDif2=0.055f;
			smallThresh=12500;
			smallMult=2.5f;
			
			maxKmerDif1=0.004f;
			maxDepthRatio1=1.1f;
			maxGCDif1=0.02f;
			maxCovariance1=0.0001f;
			maxCovariance2=0.0018f;
		}else {//9+
			maxKmerDif2=0.023f;
			maxDepthRatio2=1.7f;
			maxGCDif2=0.055f;
			smallThresh=13000;
			smallMult=2.5f;
			
			maxKmerDif1=0.004f;
			maxDepthRatio1=1.1f;
			maxGCDif1=0.02f;
			maxCovariance1=0.0002f;
			maxCovariance2=0.0022f;
		}
	}
	
	void printThresholds(){
		System.err.println("maxKmerDif2 =    "+maxKmerDif2);
		System.err.println("maxKmerDif2 =    "+maxKmerDif2);
		System.err.println("maxDepthRatio2 = "+maxDepthRatio2);
		System.err.println("maxDepthRatio2 = "+maxDepthRatio2);
		System.err.println("maxGCDif2 =      "+maxGCDif2);
		System.err.println("maxGCDif2 =      "+maxGCDif2);
		System.err.println("maxCovariance1 = "+maxCovariance1);
		System.err.println("maxCovariance2 = "+maxCovariance2);
		System.err.println("smallThresh =    "+smallThresh);
		System.err.println("smallMult =      "+smallMult);
		System.err.println("bigThresh =      "+bigThresh);
		System.err.println("bigMult =        "+bigMult);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------     Classes and Threading    ----------------*/
	/*--------------------------------------------------------------*/
	
	private void launchThreads(ArrayList<? extends Bin> list, BinMap map, int mode,
				float maxKmerDif, float maxDepthRatio, float maxGCDif, float maxProduct, float maxCovariance,
				int taxLevel, boolean allowNoTaxID, boolean allowHalfTaxID, int range, int minSize) {

		//Do anything necessary prior to processing
		
		//Determine how many threads may be used
		final int threads=Tools.mid(1, list.size()/64, Shared.threads());
		
		//Fill a list with LoadThreads
		ArrayList<CompareThread> alpt=new ArrayList<CompareThread>(threads);
		for(int i=0; i<threads; i++){
			CompareThread lt=new CompareThread(list, map, i, threads, mode,
					maxKmerDif, maxDepthRatio, maxGCDif, maxProduct, maxCovariance,
					taxLevel, allowNoTaxID, allowHalfTaxID, range, minSize);
			alpt.add(lt);
		}
		
		//Start the threads and wait for them to finish
		boolean success=ThreadWaiter.startAndWait(alpt, this);
		
		//Do anything necessary after processing
	}
	
	@Override
	public void accumulate(CompareThread t) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public ReadWriteLock rwlock() {return null;}

	@Override
	public boolean success() {return true;}
	
	class CompareThread extends Thread {
		
		
		public CompareThread(ArrayList<? extends Bin> list, BinMap map_, int tid_, int threads_, int mode_,
				float maxKmerDif_, float maxDepthRatio_, float maxGCDif_, float maxProduct_, float maxCovariance_,
				int taxLevel_, boolean allowNoTaxID_, boolean allowHalfTaxID_, int range_, int minSize_) {
			input=list;
			map=map_;
			tid=tid_;
			threads=threads_;
			mode=mode_;
			maxKmerDif=maxKmerDif_;
			maxDepthRatio=maxDepthRatio_;
			maxGCDif=maxGCDif_;
			maxProduct=maxProduct_;
			maxCovariance=maxCovariance_;
			taxLevel=taxLevel_;
			allowNoTaxID=allowNoTaxID_;
			allowHalfTaxID=allowHalfTaxID_;
			range=range_;
			minSize=minSize_;
		}
		
		
		public void run() {
			synchronized(this) {
				if(mode==REFINE_MODE) {
					refine();
				}else if(mode==RESIDUE_MODE) {
					residue();
				}
			}
		}
		
		private void refine() {
			for(int i=tid; i<input.size(); i+=threads) {
				Bin a=input.get(i);
				synchronized(a) {
					a.dest=-1;
					if(a.size()>=minSizeToMerge) {
						processBin(a);
					}
				}
			}
		}
		
		private void residue() {
			for(int i=tid; i<input.size(); i+=threads) {
				Bin a=input.get(i);
				synchronized(a) {
					a.dest=-1;
					
					Cluster b=findBestResidualCluster(a, map, maxKmerDif, 
							maxDepthRatio, maxGCDif, maxProduct, maxCovariance, key, range, minSize);
					assert(a!=b);
					if(b==null) {
						b=findLinkedCluster(a, map.contigList);
						assert(a!=b);
					}
					
					if(b!=null) {
						a.dest=b.id;
						assert(a!=b);
						assert(a.cluster()!=b);
						assert(map.contigList.get(b.id).cluster==b);
					}
				}
			}
		}
		
		private void processBin(Bin a) {
			a.dest=-1;
			
			float mult=sizeAdjustMult(a.size());
			float maxDepthRatio2=1f+(maxDepthRatio-1f)*mult;
			
			int minSize2=(int)Tools.mid(minSize, a.size(), Integer.MAX_VALUE);
			Cluster b=map.findBestCluster(a, minSize2, maxKmerDif*mult, maxDepthRatio2, maxProduct*mult, 
					maxGCDif*mult, maxCovariance*mult, taxLevel, allowNoTaxID, allowHalfTaxID, key, range);
			if(b==null) {
				//do nothing
			}else {
//				assert(a.size()<=b.size());
				assert(a!=b);
				assert(a.id()!=b.id);
				assert(!b.contigSet.contains(a.id()));
				assert(a.cluster()==null || !a.cluster().contigSet.contains(b.id));
				a.dest=b.id;
			}
		}
		
		final ArrayList<? extends Bin> input;
		final BinMap map;
		final int tid;
		final int threads;
		final int mode;
		
		final float maxKmerDif;
		final float maxDepthRatio;
		final float maxGCDif;
		final float maxProduct;
		final float maxCovariance;
		final int taxLevel;
		final boolean allowNoTaxID;
		final boolean allowHalfTaxID;
		final int range;
		final int minSize;
		final Key key=new Key();
		
	}
	
	static float sizeAdjustMult(long size) {
		if(size<smallThresh) {return 1f+smallMult*(smallThresh-size)/(float)smallThresh;}
		if(size>2*hugeThresh) {return hugeMult;}
		if(size>hugeThresh) {
			float range=1f-hugeMult;
			return Tools.min(bigThresh, 1f-(size-hugeThresh)*range/hugeThresh);
		}
		if(size>2*bigThresh) {return bigMult;}
		if(size>bigThresh) {
			float range=1f-bigMult;
			return 1f-(size-bigThresh)*range/bigThresh;
		}
		return 1f;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	boolean multiThreadedCompare=true;

//	long refinementComparisons=0;
//	long refinementComparisonsSlow=0;
	BinSketcher sketcher;
	int baseRange=1;
	int residueRange=3;

	long goodMergesFollow=0;
	long badMergesFollow=0;
	long goodMergesCreate=0;
	long badMergesCreate=0;
	long goodMergesRefine=0;
	long badMergesRefine=0;
	long goodMergesResidue=0;
	long badMergesResidue=0;
	
	BinMap binMap;
	
	/*--------------------------------------------------------------*/
	/*----------------         Final Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	final Timer phaseTimer=new Timer();
	final PrintStream outstream;
	
	/*--------------------------------------------------------------*/
	/*----------------            Static Fields          ----------------*/
	/*--------------------------------------------------------------*/
	
	static final int REFINE_MODE=0;
	static final int RESIDUE_MODE=1;
	
	static float residueStringency=1.25f;
	static float productMult=0.70f;
	
	static int minEdgeWeight0=0;
	static int maxEdges=4;
	static int minEdgeWeight=2;
	static float minEdgeRatio=0.4f;

	static int hugeThresh=1500000;
	static float hugeMult=0.65f;
	static int bigThresh=100000;
	static float bigMult=0.8f;
	static int smallThresh=12000; //10000 better for 1 depth
	static float smallMult=2.4f; //2.2f is better for 1 depth
	
//	static int minSizeToCluster=0;
//	static int minSizeToRefine=500;
	/** Size of the bigger one */
	static int minSizeToCompare=1000;
	/** Size of the smaller one being compared */
	static int minSizeToMerge=1000;
//	static int minSizeToAdd=4000;//Should always be the same as minSizeToCompare...
	static int minSizeResidue=200;
	
	//Optimal selection when forming clusters
	static float maxKmerDif1=0.004f;
	static float maxDepthRatio1=1.1f;
	static float maxGCDif1=0.02f;
	static float maxCovariance1=0.0001f;
	
	//When merging clusters
	static float maxKmerDif2=0.005f; //.005 best for 1 depth; .01 best for 4 depths.
	static float maxDepthRatio2=1.25f;
	static float maxGCDif2=0.025f; //0.02f for 1 depth, 0.05 for 4 depths
	static float maxCovariance2=0.0002f;
	
}
