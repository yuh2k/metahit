package bin;

public class BinStats {
	
	BinStats(){}
	BinStats(Bin b){
		id=b.id();
		taxid=b.taxid;
		size=b.size();
		contigs=b.numContigs();
		contam=b.contam;
		complt=b.completeness;
		gc=b.gc();
		depth=b.depth();
	}
	
	int id;
	int taxid;
	long size;
	int contigs;
	float contam;
	float complt;
	float gc;
	float depth;
	
}
