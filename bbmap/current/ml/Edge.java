package ml;

import java.util.Arrays;

public class Edge {

	public Edge(long id_, Cell c, Cell d, float w) {
		id=id_;
		source=c;
		dest=d;
		weight=w;
		assert(source!=null);
		assert(dest!=null);
		assert(source!=dest);
		assert(source.id()<dest.id());//feed-forward
	}
	
	public long id() {return id;}
	
	public String toString() {
//		return String.format("C"+source.id()+"->C"+dest.id()+",%.3f,w"+id()+"=%.4f,x=%.4f; ", source.value(), weight, weight2);
		return String.format("C"+source.id()+"->C"+dest.id()+",%.3f,w"+id()+"=%.4f,x=%.4f; ", source.value(), weight, delta);
	}
	
	void clearTemp() {
		assert(!CellNet.FAST);
//		assert(false);
//		weight2=0;
		delta=0;
	}
	
	void setFrom(Edge e, boolean setDelta) {
		assert(!CellNet.FAST);
		weight=e.weight;
//		if(setDelta) {weight2=e.weight2;}
		if(setDelta) {delta=e.delta;}
	}
	
//	public float weight() {
//		assert(CellNet.FAST || weight==dest.weights[source.lpos]) : weight+", "+
//				dest.weights[source.lpos]+", "+Arrays.toString(dest.weights);
//		return dest.weights[source.lpos];
////		return weight;
//	}
	
	public float weight() {
		assert(CellNet.FAST || weight==dest.weights[source.lpos]) : weight+", "+
				dest.weights[source.lpos]+", "+Arrays.toString(dest.weights);
//		return dest.weights[source.lpos];
		return CellNet.FAST ? dest.weights[source.lpos] : weight;
	}
	
//	float adjustWeight(float invSamples) {
////		CellNet.println(String.format("w"+id()+": %.5f -> %.5f * %.2f = %.5f",
////				weight, weight2, invSamples, weight2*invSamples));
//		weight=(float)(weight2*invSamples);
////		assert(false);
//		weight2=0;
//		return weight;
//	}
	
//	void incrementWeight2(double incr) {
//		assert(!CellNet.FAST);
////		CellNet.println(String.format("w"+id()+": weight2 %.5f -> %.5f", weight2, weight2+incr));
//		weight2+=incr;
//	}
	
	void incrementDelta(double incr) {
		assert(!CellNet.FAST);
//		CellNet.println(String.format("w"+id()+": weight2 %.5f -> %.5f", weight2, weight2+incr));
		delta+=incr;
	}
	
	void setWeight(float f) {
		assert(!CellNet.FAST);
//		CellNet.println(String.format("Set w"+id()+" %.5f -> %.5f", weight, weight+f));
		weight=f;
	}
	
	void accumulate(Edge e) {
		assert(!CellNet.FAST);
		assert(e.id==id);
		assert(e!=this);
//		weight2+=e.weight2;
		delta+=e.delta;
	}
	
	final Cell source;
	final Cell dest;
	private volatile float weight;
//	volatile double weight2;
	volatile double delta;
	private final long id;
	
}
