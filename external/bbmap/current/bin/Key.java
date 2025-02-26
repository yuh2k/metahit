package bin;

import shared.Tools;

class Key implements Cloneable {
	
	public Key(float gc, float cov, float cov2) {
		setValue(gc, cov, cov2);
	}
	
	public Key() {}

	public Key set(Bin a) {
		return setValue(a.gc(), a.depth(0), a.depth(1));
	}
	
	public Key setLevel(int gcLevel_, int covLevel_, int covLevel2_) {
		gcLevel=gcLevel_;
		covLevel=covLevel_;
		covLevel2=covLevel2_;
		assert(gcLevel>=0 && gcLevel<=(int)invGCLevel);
		assert(covLevel>=0 && covLevel<=maxDepthLevel) : covLevel_;
		assert(covLevel2>=0 && covLevel2<=maxDepthLevel) : covLevel_;
		return this;
	}
	
	public Key setValue(float gc, float cov, float cov2) {
		assert(gc>=0 && gc<=1) : gc;
		assert(cov>=0) : cov;
		assert(cov2>=0) : cov;
		return setLevel(quantizeGC(gc), quantizeDepth(cov), quantizeDepth(cov2));
	}
	
	@Override
	public boolean equals(Object other) {
		return equals((Key)other);
	}
	
	public boolean equals(Key b) {
		return gcLevel==b.gcLevel && covLevel==b.covLevel && covLevel2==b.covLevel2;
	}
	
	@Override
	public int hashCode() {
		return covLevel+(covLevel2<<10)+(gcLevel<<20);
	}
	
	@Override
	public Key clone() {
		try {
			return (Key)(super.clone());
		} catch (CloneNotSupportedException e) {
			throw new RuntimeException(e);
		}
	}

	public static int quantizeDepth(float depth) {
		float yf=depth*depth*16;
		long y=(long)yf;
		int zeros=Long.numberOfLeadingZeros(y);
		int level=Tools.min(maxDepthLevel, 64-zeros);
		return level;
	}

	public static int quantizeGC(float gc) {
		return (int)(Tools.mid(0,gc,1)*invGCLevel);
	}
	
	public String toString() {
		return "("+gcLevel+","+covLevel+","+covLevel2+")";
	}
	
	int gcLevel;
	int covLevel;
	int covLevel2;
	
	static final int maxDepthLevel=26;
	static final float gcLevelWidth=0.02f;
	static final float invGCLevel=1f/gcLevelWidth;
	
}
