package bin;

import java.util.Arrays;

import shared.Parse;
import shared.Tools;

/** Mostly written by ChatGPT and modified by me */
public class SimilarityMeasures {
	
    public static void main(String[] args) {
        float[] sample1={0.1f, 0.2f, 0.3f, 0.4f};
        float[] sample2={0.1f, 0.2f, 0.4f, 0.3f};
        int[] sample1i={1, 2, 3, 4};
        int[] sample2i={1, 2, 4, 3};
        int[] sample3i={2, 4, 6, 8};
        int[] sample4i={8, 6, 4, 2};

        // Print the similarity vector
        System.out.println("Difference Vector Float12: "+Arrays.toString(calculateDifferenceVector(sample1, sample2)));
        System.out.println("Difference Vector Int12:   "+Arrays.toString(calculateDifferenceVector(sample1i, sample2i)));
        System.out.println("Difference Vector Int13:   "+Arrays.toString(calculateDifferenceVector(sample1i, sample3i)));
        System.out.println("Difference Vector Int14:   "+Arrays.toString(calculateDifferenceVector(sample1i, sample4i)));
    }
	
    public static boolean parse(String arg, String a, String b){
    	if(a.equals("null")){
    		//Do nothing
    	}else if(a.equals("cosine") || a.equals("cos")){
    		COSINE=Parse.parseBoolean(b);
    	}else if(a.equals("euclid") || a.equals("euc")){
    		EUCLID=Parse.parseBoolean(b);
    	}else if(a.equals("absolute") || a.equals("abs")){
    		ABSOLUTE=Parse.parseBoolean(b);
    	}else if(a.equals("jsd")){
    		JSD=Parse.parseBoolean(b);
    	}else if(a.equals("hellinger") || a.equals("hell") || a.equals("hel")){
    		HELLINGER=Parse.parseBoolean(b);
    	}else if(a.equals("ks") || a.equals("kst")){
    		KST=Parse.parseBoolean(b);
    	}else {
    		return false;
    	}
    	
    	return true;
    }

    public static float[] calculateDifferenceVector(float[] a, float[] b) {
//        float cosineSimilarity=cosineSimilarity(a, b);
        float cosineDifference=cosineDifference(a, b);
        float euclideanDistance=euclideanDistance(a, b);
        float absoluteDifference=absDif(a, b);
        float jensenShannonDivergence=jensenShannonDivergence(a, b);
        float hellingerDistance=hellingerDistance(a, b);
        float ksDifference=ksTest(a, b);

        return new float[] {
            cosineDifference,
            euclideanDistance,
            absoluteDifference,
            jensenShannonDivergence,
            hellingerDistance,
            ksDifference
        };
    }

    //For setting thresholds before neural net is implemented
    public static float calculateDifferenceAverage(int[] a, int[] b) {
    	float inva=1f/Tools.max(1, Tools.sum(a));
    	float invb=1f/Tools.max(1, Tools.sum(b));
        float cosineDifference=(COSINE ? cosineDifference(a, b, inva, invb) : 0);
        float euclideanDistance=(EUCLID ? euclideanDistance(a, b, inva, invb) : 0);
        float absoluteDifference=(ABSOLUTE ? absDif(a, b, inva, invb) : 0);
        float jensenShannonDivergence=(JSD ? jensenShannonDivergence(a, b, inva, invb) : 0);
        float hellingerDistance=(HELLINGER? hellingerDistance(a, b, inva, invb) : 0);
        float ksDifference=(KST ? ksTest(a, b, inva, invb) : 0);
        int div=(COSINE ? 1 : 0)+(EUCLID ? 1 : 0)+(ABSOLUTE ? 1 : 0)+(JSD ? 1 : 0)+(HELLINGER ? 1 : 0)+(KST ? 1 : 0);
        return (cosineDifference+euclideanDistance+absoluteDifference+
        		jensenShannonDivergence+hellingerDistance+ksDifference)/div;
    }

    public static float[] calculateDifferenceVector(int[] a, int[] b) {
    	float inva=1f/Tools.max(1, Tools.sum(a));
    	float invb=1f/Tools.max(1, Tools.sum(b));
//        float cosineSimilarity=cosineSimilarity(a, b, inva, invb);
        float cosineDifference=cosineDifference(a, b, inva, invb);
        float euclideanDistance=euclideanDistance(a, b, inva, invb);
        float absoluteDifference=absDif(a, b, inva, invb);
        float jensenShannonDivergence=jensenShannonDivergence(a, b, inva, invb);
        float hellingerDistance=hellingerDistance(a, b, inva, invb);
        float ksDifference=ksTest(a, b, inva, invb);

        return new float[] {
                cosineDifference,
            euclideanDistance,
            absoluteDifference,
            jensenShannonDivergence,
            hellingerDistance,
            ksDifference
        };
    }

    public static float cosineDifference(float[] a, float[] b) {
    	return 1-cosineSimilarity(a, b);
    }
    
    public static float cosineSimilarity(float[] a, float[] b) {
        float dotProduct=0f;
        float normVec1=0f;
        float normVec2=0f;

        for (int i=0; i<a.length; i++) {
        	float ai=a[i], bi=b[i];
            dotProduct+=ai*bi;
            normVec1+=ai*ai;
            normVec2+=bi*bi;
        }

        return (float)(dotProduct/(Math.sqrt(normVec1)*Math.sqrt(normVec2)));
    }

    public static float cosineDifference(int[] a, int[] b, float inva, float invb) {
    	return 1-cosineSimilarity(a, b, inva, invb);
    }
    
    public static float cosineSimilarity(int[] a, int[] b, float inva, float invb) {
        float dotProduct=0f;
        float normVec1=0f;
        float normVec2=0f;

        for (int i=0; i<a.length; i++) {
        	float ai=a[i]*inva, bi=b[i]*invb;
            dotProduct+=ai*bi;
            normVec1+=ai*ai;
            normVec2+=bi*bi;
        }

        return (float)(dotProduct/(Math.sqrt(normVec1)*Math.sqrt(normVec2)));
    }

    public static float euclideanDistance(float[] a, float[] b) {
        float sumSquaredDifferences=0f;

        for (int i=0; i<a.length; i++) {
        	float ai=a[i], bi=b[i];
        	float d=ai-bi;
            sumSquaredDifferences+=d*d;
        }

        return (float)Math.sqrt(sumSquaredDifferences);
    }

    public static float euclideanDistance(int[] a, int[] b, float inva, float invb) {
        float sumSquaredDifferences=0f;

        for (int i=0; i<a.length; i++) {
        	float ai=a[i]*inva, bi=b[i]*invb;
        	float d=ai-bi;
            sumSquaredDifferences+=d*d;
        }

        return (float)Math.sqrt(sumSquaredDifferences);
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

    public static float jensenShannonDivergence(float[] a, float[] b) {
        float kldSumA=0, kldSumB=0;
        for (int i=0; i<a.length; i++) {
        	float ai=a[i]+0.0005f, bi=b[i]+0.0005f;//Prevents zero values
        	float avgi=(ai+bi)*0.5f;
            kldSumA+=ai*Math.log(ai/avgi);
            kldSumA+=bi*Math.log(bi/avgi);
        }
        return (kldSumA+kldSumB)*invLog2*0.5f;
    }

    public static float jensenShannonDivergence(int[] a, int[] b, float inva, float invb) {
        float kldSumA=0, kldSumB=0;
        for (int i=0; i<a.length; i++) {
        	float ai=a[i]*inva+0.0005f, bi=b[i]*invb+0.0005f;//Prevents zero values
        	float avgi=(ai+bi)*0.5f;
            kldSumA+=ai*Math.log(ai/avgi);
            kldSumA+=bi*Math.log(bi/avgi);
        }
        return (kldSumA+kldSumB)*invLog2*0.5f;
    }

//    public static float jensenShannonDivergenceSlow(float[] a, float[] b) {
//        float[] avg=new float[a.length];
//        for (int i=0; i<a.length; i++) {
//        	float ai=a[i], bi=b[i];
//            avg[i]=(ai+bi)*0.5f;
//        }
//
//        return (kullbackLeiblerDivergence(a, avg)+kullbackLeiblerDivergence(b, avg))*0.5f;
//    }
//
//    public static float kullbackLeiblerDivergence(float[] p, float[] q) {
//        float sum=0f;
//        for (int i=0; i<p.length; i++) {
//        	float pi=p[i], qi=q[i];
//            if (p[i]!=0) {
//                sum+=p[i]*Math.log(pi/qi);
//            }
//        }
//        return sum*invLog2;
//    }

    public static float hellingerDistance(float[] a, float[] b) {
        float sum=0f;
        for (int i=0; i<a.length; i++) {
        	float ai=a[i], bi=b[i];
        	float d=(float)(Math.sqrt(ai)-Math.sqrt(bi));
            sum+=d*d;
        }
        return (float)Math.sqrt(sum)*invRoot2;
    }

    public static float hellingerDistance(int[] a, int[] b, float inva, float invb) {
        float sum=0f;
        for (int i=0; i<a.length; i++) {
        	float ai=a[i]*inva, bi=b[i]*invb;
        	float d=(float)(Math.sqrt(ai)-Math.sqrt(bi));
            sum+=d*d;
        }
        return (float)Math.sqrt(sum)*invRoot2;
    }
    
    /** This is a KS test for binned histograms, not raw values */
    public static float ksTest(float[] histogram1, float[] histogram2) {
        // Ensure both histograms have the same length
        if (histogram1.length!=histogram2.length) {
            throw new IllegalArgumentException("Histograms must have the same number of bins");
        }

        float cd1=0, cd2=0, dMax=0;

        // Compute the KS statistic (maximum absolute difference between the two CDFs)
        for (int i=0; i<histogram1.length; i++) {
        	cd1+=histogram1[i];
        	cd2+=histogram2[i];
            dMax=(float)Math.max(dMax, Math.abs(cd1-cd2));
        }

        return dMax;
    }
    
    /** This is a KS test for binned histograms, not raw values */
    public static float ksTest(int[] a, int[] b, float inva, float invb) {
        // Ensure both histograms have the same length
        if (a.length!=b.length) {
            throw new IllegalArgumentException("Histograms must have the same number of bins");
        }

        float cda=0, cdb=0, dMax=0;

        // Compute the KS statistic (maximum absolute difference between the two CDFs)
        for (int i=0; i<a.length; i++) {
        	float ai=a[i]*inva, bi=b[i]*invb;
        	cda+=ai;
        	cdb+=bi;
            dMax=(float)Math.max(dMax, Math.abs(cda-cdb));
        }

        return dMax;
    }

    private static final float root2=(float)Math.sqrt(2);
    private static final float log2=(float)Math.log(2);
    private static final float invRoot2=1/root2;
    private static final float invLog2=1/log2;

    //2531 kcps (times include contig loading)
    //26 clusters
    public static boolean COSINE=true;
    //2796 kcps
    //21 clusters
    public static boolean EUCLID=false;
    //2636 kcps
    //23 clusters at 4x threshold of cosine
    public static boolean ABSOLUTE=false;
    //183 kcps
    //22 clusters
    public static boolean JSD=false;
    //953 kcps
    //~22 at 2x threshold of cosine
    public static boolean HELLINGER=false;
    //1859 kcps
    //20 clusters 
    public static boolean KST=false;
    
}
