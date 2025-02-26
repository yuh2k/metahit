package bin;

import java.io.PrintStream;
import java.util.Arrays;
import java.util.Collection;

import dna.AminoAcid;
import shared.Tools;
import structures.IntHashMap;
import tax.TaxTree;

/**
 * Superclass for binner classes
 * @author Brian Bushnell
 * @date Dec 4, 2014
 *
 */
public class BinObject {
	
	public static void setQuant(int x) {
		quant=x;
		assert(quant>0);
		setK(k);
	}
	
	public static void setK(int x) {
		assert(k>0 && k<16);
		k=x;
		remap=makeRemap(k);
		canonicalKmers=Tools.max(remap)+1;
		invCanonicalKmers=1f/canonicalKmers;
	}
	
	public static int[] makeRemap(int k){
		final int bits=2*k;
		final int max=(int)((1L<<bits)-1);
		int count=0;
		IntHashMap map=new IntHashMap();
		for(int kmer=0; kmer<=max; kmer++){
			int canon=Tools.min(kmer, AminoAcid.reverseComplementBinaryFast(kmer, k));
			if(canon%quant==0 && !map.containsKey(canon)) {
				map.put(canon, count);
				count++;
			}
			int idx=map.get(canon);
			map.put(kmer, idx);
		}
		int[] remap=new int[max+1];
		Arrays.fill(remap, -1);
		for(int kmer=0; kmer<=max; kmer++){
			remap[kmer]=map.get(kmer);
		}
		return remap;
	}
	
	public static int countKmers(final byte[] bases, final int[] counts){
		if(quant>1) {return countKmers_quantized(bases, counts);}
		if(bases==null || bases.length<k){return 0;}
//		counts=(counts!=null ? counts : new int[canonicalKmers]);
		
		final int shift=2*k;
//		final int shift2=shift-2;
		final int mask=~((-1)<<shift);
		
		int kmer=0;
//		int rkmer=0;
		int len=0;
		int valid=0;
		
		for(int i=0; i<bases.length; i++){
			byte b=bases[i];
			int x=AminoAcid.baseToNumber[b];
//			int x2=AminoAcid.baseToComplementNumber[b];
			kmer=((kmer<<2)|x)&mask;
//			rkmer=((rkmer>>>2)|(x2<<shift2))&mask;
			if(x>=0){
				len++;
				if(len>=k) {
					valid++;
					counts[remap[kmer]]++;
				}
			}else{len=kmer=0;}
		}
		return valid;
	}
	
	public static int countKmers_quantized(final byte[] bases, final int[] counts){
		if(bases==null || bases.length<k){return 0;}
//		counts=(counts!=null ? counts : new int[canonicalKmers]);
		
		final int shift=2*k;
		final int mask=~((-1)<<shift);
		
		int kmer=0;
		int len=0;
		int valid=0;
		
		for(int i=0; i<bases.length; i++){
			byte b=bases[i];
			int x=AminoAcid.baseToNumber[b];
			kmer=((kmer<<2)|x)&mask;
			if(x>=0){
				len++;
				if(len>=k) {
					valid++;
					int pos=remap[kmer];
					if(pos>=0) {counts[pos]++;}
				}
			}else{len=kmer=0;}
		}
		return valid;
	}
	
	/**
	 * @param a Contig kmer frequencies
	 * @param b Cluster kmer frequencies
	 * @return Score
	 */
	static final float absDif(float[] a, float[] b){
		assert(a.length==b.length);
		double sum=0;
		for(int i=0; i<a.length; i++){
			sum+=Math.abs(a[i]-b[i]);
		}

		return (float)sum;
	}
	
	static final float absDif(int[] a, int[] b){
		return absDif(a, b, 1f/Tools.sum(a), 1f/Tools.sum(b));
	}
	
	/**
	 * @param a Contig kmer counts
	 * @param b Cluster kmer counts
	 * @return Score
	 */
	static final float absDif(int[] a, int[] b, float inva, float invb){
		assert(a.length==b.length);
		float sum=0;
		for(int i=0; i<a.length; i++){
			float ai=a[i]*inva, bi=b[i]*invb;
			sum+=Math.abs(ai-bi);
		}
		return sum;
	}
	
	/**
	 * @param a Contig kmer frequencies
	 * @param b Cluster kmer frequencies
	 * @return Score
	 */
	static final float rmsDif(float[] a, float[] b){
		assert(a.length==b.length);
		double sum=0;
		for(int i=0; i<a.length; i++){
//			double d=Tools.absdif((double)a[i], (double)b[i]);
			double d=(a[i])-(b[i]);
			sum+=d*d;
		}

		return (float)Math.sqrt(sum/a.length);
	}
	
	static final float rmsDif(int[] a, int[] b){
		return rmsDif(a, b, 1f/Tools.sum(a), 1f/Tools.sum(b));
	}
	
	/**
	 * @param a Contig kmer counts
	 * @param b Cluster kmer counts
	 * @return Score
	 */
	static final float rmsDif(int[] a, int[] b, float inva, float invb){
		assert(a.length==b.length);
		long sum=0;
		for(int i=0; i<a.length; i++){
			float ai=a[i]*inva, bi=b[i]*invb;
			float d=(ai-bi);
			sum+=d*d;
		}
		return (float)Math.sqrt(sum/a.length);
	}
	
	/**
	 * @param a Contig kmer frequencies
	 * @param b Cluster kmer frequencies
	 * @return Score
	 */
	static final float ksFunction(float[] a, float[] b){
		assert(a.length==b.length);
		double sum=0;
		for(int i=0; i<a.length; i++){
			double ai=a[i]+0.0005;
			double bi=b[i]+0.0005;
			double d=(double)ai*Math.log(ai/bi);
			sum+=d;
		}
		
		return (float)sum;
	}
	
	static boolean isValid(Collection<? extends Bin> list, boolean allowImproperContigs) {
		for(Bin b : list) {
			if(b.isCluster()) {
				Cluster c=(Cluster)b;
				assert(c.isValid());
				for(Contig x : c.contigs) {assert(x.isValid());}
			}else {
				Contig c=(Contig)b;
				assert(c.isValid());
				assert(allowImproperContigs || c.cluster()==null);
				assert(c.cluster()==null || c.cluster().isValid());
			}
		}
		return true;
	}

	/** Kmer length for frequencies */
	private static int k=4;
	/** Number of canonical kmers; frequency array length */
	public static int canonicalKmers=-1;
	private static int quant=1;//Determines how many tetramers to use for comparisons
	public static float invCanonicalKmers=-1;
	/** Maps a kmer to index in frequency array */
	public static int[] remap;
	
	/** Print status messages to this output stream */
	static PrintStream outstream=System.err;
	static TaxTree tree=null;
	
	static int minBasesPerCluster=50000;
	static int minContigsPerCluster=1;
	
	static boolean verbose;
	static float sketchDensity=1/100f;
	static boolean sketchContigs=false;
	static boolean sketchClusters=false;
	static boolean sketchOutput=false;
	static boolean sketchInBulk=true;
	static int sketchSize=20000;
	
	static boolean validation=false;
	static boolean grading=false;
	static boolean parseTaxid=true;
	static boolean depthZeroProxy=true;
	static int globalTime=0;
	
	static {setK(4);}
	
}
