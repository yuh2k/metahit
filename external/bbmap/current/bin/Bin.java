package bin;

import json.JsonObject;
import shared.Tools;
import structures.ByteBuilder;
import structures.FloatList;
import structures.IntHashMap;
import structures.IntLongHashMap;

public abstract class Bin extends BinObject implements Sketchable, Iterable<Contig> {
	
	@Override
	public final int taxid() {return taxid;}

	@Override
	public final float gc() {return gcSum/(float)size();}
	
	public final void clearDepth() {
		depth.clear();
		normDepth=null;
	}
	
	public float depthRatio(Bin b) {
		float max=1;
		for(int i=0; i<depth.size; i++) {
			float d1=depth.get(i)+0.5f;
			float d2=b.depth.get(i)+0.5f;
			float ratio=Tools.max(d1,d2)/Tools.min(d1,d2);
			max=Tools.max(max, ratio);
		}
		return max;
	}
	
	public final void setDepth(float d, int sample) {
		depth.set(sample, d);
	}
	
	public final void appendDepth(float d) {
		depth.add(d);
	}
	
	public final int numDepths() {
		return depth.size;
	}
	
	public final float depth(int sample) {
		return (sample==1 && depth.size==1 ? 0 : depth.get(sample));
	}
	
	public float[] normDepth() {
		if(depth.size()<2) {return null;}
		if(normDepth==null) {fillNormDepth();}
		return normDepth;
	}
	
	void fillNormDepth() {
		assert(normDepth==null || (normDepth.length>1 && normDepth.length==numDepths()));
		if(normDepth==null) {normDepth=new float[depth.size];}
		float sum=0;
		for(int i=0; i<depth.size; i++) {
			float f=depth.get(i);
			f=(float)Math.log(f+1);
			sum+=f;
			normDepth[i]=f;
		}
		float inv=1/Tools.max(sum, 0.1f);
		for(int i=0; i<normDepth.length; i++) {
			normDepth[i]*=inv;
		}
	}
	
	/** Uses a weighted sum of linear and geometric means */
	public final float depth() {
		if(depthZeroProxy) {return depth.get(0);}
		if(avgDepthValid) {return avgDepth;}
		synchronized(this) {
			if(depth.size()==1) {avgDepth=depth.get(0);}
			else {
				double product=1;
				double sum=0;
				for(int i=0; i<depth.size; i++) {
					float d=depth.get(i);
					product*=(d+0.25f);
					sum+=d;
				}
				double inv=1.0/depth.size;
				float geo=(float)(Math.pow(product, inv)-0.25);
				float linear=(float)(sum*inv);
				avgDepth=geo*0.75f+linear*0.25f;
			}
			avgDepthValid=true;
		}
		return avgDepth;
	}
	
	@Override
	/** Biggest first */
	public final int compareTo(Sketchable o) {
		if(size()!=o.size()) {return size()>o.size() ? -1 : 1;}//Biggest first
		return o.id()-id();
	}
	
	@Override
	public final void setFrom(JsonObject all) {
		assert(sketchedSize<size());
		clearTax();
		JsonObject top=null, second=null;
		if(all!=null && all.jmapSize()>0) {
			for(String key : all.jmap.keySet()){
				JsonObject hit=all.jmap.get(key);
				if(top==null) {top=hit;}
				else {
					if(hit.getLong("TaxID")!=1806490) {//Achromobacter sp. ATCC35328; messes with E.coli.
						second=hit;
						break;
					}
				}
			}
		}
		topHit=(top==null ? null : new SketchRecord(top));
		secondHit=(second==null ? null : new SketchRecord(second));
		taxid=(topHit==null ? -1 : topHit.taxid);
		genusTaxid=(topHit==null ? -1 : topHit.genusTaxid);
		sketchedSize=size();
	}
	
	@Override
	public final void clearTax() {
		taxid=genusTaxid=-1;
		topHit=secondHit=null;
		sketchedSize=0;
	}
	
	@Override
	public final String toString() {
		return toBytes().toString();
	}
	
	public final ByteBuilder toBytes() {
		ByteBuilder bb=new ByteBuilder();
		bb.append(isCluster() ? "Cluster " : "Contig ").append(id()).append(":");
		bb.tab().append("Size ").append(size());
		bb.tab().append("Contigs ").append(numContigs());
		bb.tab().append("GC ").append(gc(), 3);
		bb.tab().append("Depth ").append(depth(), 1);
		if(depth.size()>1) {
			for(int i=0; i<depth.size; i++) {bb.comma().append(depth(i), 1);}
		}
		bb.tab().append("TaxID ").append(taxid);
		if(validation) {
			bb.tab().append("TaxID0 ").append(labelTaxid);
			if(completeness>=0) {
				bb.tab().append("Complt ").append(completeness*100, 2);
				bb.tab().append("Contam ").append(contam*100, 2);
			}
		}
//		if(labelTaxid>0) {bb.tab().append("TaxID0 ").append(labelTaxid);}
		if(topHit!=null) {topHit.appendTo(bb.nl().tab().tab());}
		if(secondHit!=null) {secondHit.appendTo(bb.nl().tab().tab());}
		return bb;
	}
	
	/** Higher is more similar */
	public final float similarityTo(Bin b) {
		final float ratio=depthRatio(b);
		final float gc=gc(), gc2=b.gc();
		final float gcDif=Math.abs(gc-gc2)+1f;
		final float simDif=SimilarityMeasures.calculateDifferenceAverage(counts, b.counts)*0.5f+1f;
		final float covariance=1+covariance(b)*32;
		float product=simDif*ratio*gcDif*covariance;
		return 1f/product;
	}
	
	public final boolean canMergeWith(Bin b, float stringency) {
		long size=Tools.min(size(), b.size());
		float mult=Binner.sizeAdjustMult(size);
		stringency*=mult;
		
		float maxKmerDif=Binner.maxKmerDif2*stringency;
		float maxDepthRatio=1+((Binner.maxDepthRatio2-1)*stringency);
		float maxGCDif=Binner.maxGCDif2*stringency;
		float maxProduct=maxKmerDif*maxDepthRatio*Binner.productMult;
		float maxCovariance=Binner.maxCovariance2*stringency;
		return canMergeWith(b, maxGCDif, maxDepthRatio, maxKmerDif, maxProduct, maxCovariance);
	}
	
	/** Higher is more similar */
	public final boolean canMergeWith(Bin b, float maxGCDif, float maxDepthRatio, 
			float maxKmerDif, float maxProduct, float maxCovariance) {
		long edges1=countEdgesTo(b);
		long edges2=b.countEdgesTo(this);
//		if(Tools.min(edges1, edges2)<Binner.minEdgeWeight0) {return false;}//Too fragmented
		float mult=(edges1>1 ? 1.4f : 1f)*(edges2>1 ? 1.4f : 1f);
		float gcDif=Math.abs(gc()-b.gc());
		if(gcDif>maxGCDif*mult) {return false;}
		final float depthRatio=depthRatio(b);
		final float covariance=covariance(b);
		if(depthRatio>maxDepthRatio*mult || covariance>maxCovariance*mult) {return false;}
		final float kmerDif=SimilarityMeasures.calculateDifferenceAverage(counts, b.counts);
		final float product=kmerDif*depthRatio;
		if(kmerDif>maxKmerDif*mult || product>maxProduct*mult) {return false;}
		return true;
	}
	
	public float covariance(Bin b) {
		if(depth.size()<2) {return 0;}
		float f=SimilarityMeasures.cosineDifference(normDepth(), b.normDepth());
		return f;
	}
	
	public final long sketchedSize() {return sketchedSize;}
	
	public final void calcContam(IntLongHashMap sizeMap) {
		IntLongHashMap taxmap=new IntLongHashMap(7);
		long sum=0;
		for(Contig c : this) {
			int tid=c.labelTaxid;
			taxmap.increment(tid, c.size());
			sum+=c.size();
		}
		assert(sum==size());
		int[] keys=taxmap.keys();
		long[] values=taxmap.values();
		final int invalid=taxmap.invalid();
		int tid=-1;
		long maxSize=-1;
		for(int i=0; i<keys.length; i++) {
			int key=keys[i];
			long value=values[i];
			if(key!=invalid && value>maxSize) {
				tid=key;
				maxSize=value;
			}
		}
		taxid=tid;
		long targetSize=sizeMap.get(tid);
		if(targetSize==sizeMap.invalid()) {targetSize=sum;}//unknown...
		completeness=maxSize/(float)targetSize;
		contam=(sum-maxSize)/(float)sum;
	}
	
	abstract boolean sameCluster(Bin b);
	
	public abstract boolean isCluster();
	
	public abstract Cluster toCluster();
	
	public abstract Cluster cluster();
	
	public abstract boolean isValid();
	
	public final boolean isEmpty() {return numContigs()<1;}
	
	public long countEdgesTo(Bin b) {
		if(!b.isCluster()) {return countEdgesTo((Contig)b);}
		else {return countEdgesTo((Cluster)b);}
	}
	
	public long countEdgesTo(Contig b) {
		return pairMap==null ? 0 : Tools.max(0, pairMap.get(b.id()));
	}
	
	public long countEdgesTo(Cluster b) {
		if(pairMap==null) {return 0;}
		int[] keys=pairMap.keys(), values=pairMap.values();
		long sum=0;
		for(int i=0, invalid=pairMap.invalid(); i<keys.length; i++) {
			int key=keys[i];
			if(key!=invalid && b.contigSet.contains(key)) {
				sum+=values[i];
			}
		}
		return sum;
	}
	
	public int kmers;
	public float invKmers;
	
	public int[] counts;
	public long gcSum;
	public long sketchedSize;
	
	public long depthSum=0;
	private float avgDepth=-1;
	boolean avgDepthValid=false;
	private FloatList depth=new FloatList(1);
	private float[] normDepth;
	public IntHashMap pairMap;
	float completeness=0, contam=0;
	float entropy;
	
	int dest=-1;
	
	public int taxid;
	public int genusTaxid;
	public int labelTaxid;//For validation on labeled data
	SketchRecord topHit;
	SketchRecord secondHit;
	
}
