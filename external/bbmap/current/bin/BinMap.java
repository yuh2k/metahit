package bin;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.Map.Entry;
import java.util.concurrent.ConcurrentHashMap;

import shared.Tools;

public class BinMap extends BinObject implements Iterable<Cluster> {
	
	public BinMap(ArrayList<Contig> contigs) {contigList=contigs;}
	
	/** Higher is less stringent; 1.0 is neutral, 0 is exact match */
	public Cluster addOrMerge(Bin a, int minSizeToCompare, int minSizeToMerge, int minSizeToAdd, 
			float maxKmerDif, float maxDepthRatio, float maxGCDif, float maxCovariance, float stringency, 
			int taxlevel, boolean allowNoTaxID, boolean allowHalfTaxID, Key key, int matrixRange) {
		if(minSizeToMerge>=0 && a.size()>=minSizeToMerge) {
			maxKmerDif=maxKmerDif*stringency;
			maxDepthRatio=1+((maxDepthRatio-1)*stringency);
			assert(maxDepthRatio>=1) : maxDepthRatio+", "+stringency;
			maxGCDif*=stringency;
			float maxProduct=maxKmerDif*maxDepthRatio*Binner.productMult;
			Cluster best=findBestCluster(a, minSizeToCompare, maxKmerDif, maxDepthRatio, maxProduct, maxGCDif,
					maxCovariance, taxlevel, allowNoTaxID, allowHalfTaxID, key, matrixRange);
			if(best!=null) {
				best.add(a);
				return best;
			}
		}
		if(minSizeToAdd>=0 && a.size()>=minSizeToAdd) {return add(a, key.set(a));}
		residual.add(a);
		return null;
	}
	
	void addAll(Collection<? extends Bin> bins, int minSize) {
		Key key=new Key();
		for(Bin b : bins) {
			if(b.size()<minSize) {
				residual.add(b);
			}else {
				add(b, key);
			}
		}
	}
	
	Cluster add(Bin a, Key key) {
		Cluster c=a.toCluster();
		if(key==null) {key=new Key();}
		key.set(a);
		ArrayList<Cluster> list=getOrMakeList(key);
		synchronized(list) {list.add(c);}
		return c;
	}
	
	public ArrayList<Cluster> getOrMakeList(Key key){
		ArrayList<Cluster> list=map.get(key);
		if(list==null) {
			map.putIfAbsent((Key)(key.clone()), new ArrayList<Cluster>(8));
			list=map.get(key);
		}
		return list;
	}
	
	public Cluster findBestCluster(Bin a, int minSizeToCompare, float maxKmerDif, float maxDepthRatio, 
			float maxProduct, float maxGCDif, float maxCovariance, int taxlevel, boolean allowNoTaxID, 
			boolean allowHalfTaxID, Key key, int matrixRange) {
		if(key==null) {key=new Key();}
		final float gc=a.gc();
		final float depth=a.depth(0);
		final float depth2=a.depth(1);
		key.set(a);
		assert(maxDepthRatio>=1) : maxDepthRatio;
		final int minDepthLevel=Tools.max(0, key.covLevel-matrixRange, Key.quantizeDepth(depth/maxDepthRatio));
		final int maxDepthLevel=Tools.min(key.covLevel+matrixRange, Key.quantizeDepth(depth*maxDepthRatio));
		final int minGCLevel=Tools.max(0, key.gcLevel-matrixRange, Key.quantizeGC(gc-maxGCDif));
		final int maxGCLevel=Tools.min(key.gcLevel+matrixRange, Key.quantizeGC(gc+maxGCDif));
		assert(minGCLevel<=maxGCLevel && minDepthLevel<=maxDepthLevel) : "mingc="+minGCLevel+", maxgc="+maxGCLevel+
			", mind="+minDepthLevel+", maxd="+maxDepthLevel+"\nrange="+matrixRange+", "+key;
		
		int minDepthLevel2=0;
		int maxDepthLevel2=0;
		if(a.numDepths()>1) {
			minDepthLevel2=Tools.max(0, key.covLevel2-matrixRange, Key.quantizeDepth(depth2/maxDepthRatio));
			maxDepthLevel2=Tools.min(key.covLevel2+matrixRange, Key.quantizeDepth(depth2*maxDepthRatio));
		}
		
		Cluster best=null;
		float bestSimilarity=-1;
		for(int depthLevel=minDepthLevel; depthLevel<=maxDepthLevel; depthLevel++) {
			for(int depthLevel2=minDepthLevel2; depthLevel2<=maxDepthLevel2; depthLevel2++) {
				for(int gcLevel=minGCLevel; gcLevel<=maxGCLevel; gcLevel++) {
					if(verbose) {System.err.println("Looking at depth "+depthLevel+", gc "+gcLevel);}
					key.setLevel(gcLevel, depthLevel, depthLevel2);
					ArrayList<Cluster> list=map.get(key);
					int idx=-1;
					if(list!=null && !list.isEmpty()) {
						idx=findBestBinIndex(a, list, minSizeToCompare, maxKmerDif, maxDepthRatio, 
								maxProduct, maxGCDif, maxCovariance, taxlevel, allowNoTaxID, allowHalfTaxID);}

					if(idx>=0) {
						Cluster b=list.get(idx);
						if(verbose) {System.err.println("Found idx "+idx+"; b="+b.id()+" ("+b.size+"bp)");}
						float similarity=0;
						if(a.size()<b.size) {similarity=a.similarityTo(b);}
						else {similarity=b.similarityTo(a);}
						if(similarity>bestSimilarity) {
							bestSimilarity=similarity;
							best=b;
							if(verbose) {System.err.println("Set best to "+b.id());}
						}
					}
				}
			}
		}
		return best;
	}
	
	private int findBestBinIndex(Bin a, ArrayList<? extends Bin> clusters, 
			int minSizeToCompare, float maxKmerDif, float maxDepthRatio, float maxProduct, float maxGCDif, 
			float maxCovariance, int taxlevel, boolean allowNoTaxID, boolean allowHalfTaxID) {
		
		int bestIdx=-1;
		float bestSimilarity=0;
		for(int i=0; i<clusters.size(); i++) {
//			comparisons++;
			Bin b=clusters.get(i);
			if(b==null || a==b) {continue;}
			if(b.size()<minSizeToCompare) {break;}
			if(!allowHalfTaxID && (a.taxid<1 || b.taxid<1)) {continue;}
			if(!allowNoTaxID && a.taxid<1 && b.taxid<1) {continue;}
			if(taxlevel>=0 && tree!=null && a.taxid!=b.taxid && a.taxid>0 && b.taxid>0) {
				int commonAncestorLevel=tree.commonAncestorLevel(a.taxid, b.taxid);
				if(commonAncestorLevel>taxlevel) {continue;}
			}
			
			if(!a.canMergeWith(b, maxGCDif, maxDepthRatio, maxKmerDif, maxProduct, maxCovariance)) {
				continue;
			}
			
//			slowComparisons++;
			float similarity=a.similarityTo(b);
			if(similarity>bestSimilarity) {
				bestIdx=i;
				bestSimilarity=similarity;
			}
		}
		return bestIdx;
	}
	
	ArrayList<Cluster> toList(boolean addResidue){
		ArrayList<Cluster> list=new ArrayList<Cluster>();
		for(Entry<Key, ArrayList<Cluster>> e : map.entrySet()) {
			list.addAll(e.getValue());
		}
		if(addResidue) {
			for(Bin b : residual) {
				list.add(b.cluster()==null ? b.toCluster() : b.cluster());
			}
		}
		return list;
	}
	
	@Override
	public Iterator<Cluster> iterator() {
		return toList(false).iterator();
	}
	
	public void clear(boolean clearResidual) {
		map.clear();
		if(clearResidual) {residual.clear();}
	}
	
	public boolean isValid() {
		assert(isValid(contigList, true));
		assert(isValid(residual, false));
		for(ArrayList<Cluster> list : map.values()) {
			assert(isValid(list, false));
		}
		return true;
	}
	
	public int countClusters() {
		int sum=0;
		for(ArrayList<Cluster> list : map.values()) {sum+=list.size();}
		return sum;
	}
	
	public ConcurrentHashMap<Key, ArrayList<Cluster>> map=new ConcurrentHashMap<Key, ArrayList<Cluster>>(2000, 0.7f, 32);
	public ArrayList<Bin> residual=new ArrayList<Bin>();//Should really be contigs
	public ArrayList<Contig> contigList;
	
//	public long comparisons=0;
//	public long slowComparisons=0;
	
}
