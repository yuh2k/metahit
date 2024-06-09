package ml;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Locale;
import java.util.Random;

import shared.Shared;
import shared.Tools;
import shared.Vector;
import structures.ByteBuilder;
import structures.FloatList;

public class CellNet implements Cloneable, Comparable<CellNet> {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	public CellNet() {//Dummy network
		seed=0;
		dims=null;
		layers=0;
		values=null;
		weights=weightsOut=null;
		eOverNet=null;
		list=null;
		net=null;
		finalLayer=null;
	}
	
	public CellNet(int[] dims_, long seed_, ArrayList<String> commands_) {
		seed=(seed_>=0 ? seed_ : new Random().nextLong()&Long.MAX_VALUE);
		commands=commands_;
		dims=dims_.clone();
		layers=dims.length;
		values=makeFloatMatrix(dims);
		eOverNet=makeFloatMatrix(dims);
		
		int cells=(int) shared.Vector.sum(dims);
		list=new ArrayList<Cell>(cells+1);
		net=makeNodes(dims, list, values, eOverNet);
		finalLayer=net[layers-1];
		
		weights=makeWeights();
		weightsOut=makeWeightsOut();
	}
	
	private float[][][] makeWeights(){
		float[][][] x=new float[dims.length][][];
		for(int layer=1; layer<x.length; layer++) {
			Cell[] current=net[layer];
			float[][] y=x[layer]=new float[current.length][];
			for(int i=0; i<y.length; i++) {
				y[i]=current[i].weights;
			}
		}
		return x;
	}
	
	private float[][][] makeWeightsOut(){
		float[][][] x=new float[dims.length][][];
		for(int layer=0; layer<x.length-1; layer++) {
			Cell[] current=net[layer];
			Cell[] next=net[layer+1];
			float[][] y=x[layer]=new float[current.length][];
			for(int i=0; i<y.length; i++) {
				y[i]=new float[next.length];
			}
		}
		return x;
	}
	
//	public void initFast() {
//		for(Cell c : list) {if(c!=null) {c.initFast();}}
//	}
//	
//	public void loadFast() {
//		for(Cell c : list) {if(c!=null) {c.loadFast();}}//TODO: Skip this for layer 0...
//	}
	
	public void randomize(float density) {
		{
			Random randy=Shared.threadLocalRandom(seed);
			if(FAST) {
				makeEdgesFast(net, list, randy, density);
			}else {
				makeEdgesSlow(net, list, randy, density);
			}
		}
		assert((Function.TYPE_RATES_CUM!=null)==(shared.Vector.sum(Function.TYPE_RATES)>0)) : Arrays.toString(Function.TYPE_RATES);
		if(Function.TYPE_RATES_CUM!=null || Cell.randomTypeRate>0) {
			Random randy=Shared.threadLocalRandom(seed);
			if(Function.TYPE_RATES_CUM!=null) {
				randomizeActivationA(net, list, randy, Function.TYPE_RATES_CUM);
			}else{
				randomizeActivationB(net, list, randy, Cell.randomTypeRate);
			}
		}
	}
	
	static void randomizeActivationA(Cell[][] net, ArrayList<Cell> list, Random randy, float[] cumRate) {
		assert(cumRate!=null && cumRate[cumRate.length-1]==1) : Arrays.toString(cumRate);
		
		for(int layerNum=1; layerNum<net.length-1; layerNum++) {//Only middle layers
			Cell[] layer=net[layerNum];
			for(Cell c : layer) {
//				c.type=randomType(randy, cumRate);
				c.function=Function.randomFunction(randy);
			}
		}
	}
	
	static void randomizeActivationB(Cell[][] net, ArrayList<Cell> list, Random randy, float rate) {
		assert(rate>0);
//		final int types=Cell.TYPES.length;
		final int types=Tools.min(Cell.MAX_TYPE+1, Function.TYPES.length);
//		assert(false) : types;
//		System.err.println("types="+types);
		for(int layerNum=1; layerNum<net.length-1; layerNum++) {//Only middle layers
			Cell[] layer=net[layerNum];
			for(Cell c : layer) {
				if(randy.nextFloat()<rate) {
//					int t=c.type;
//					while(t==c.type) {t=randy.nextInt(types);}
//					c.type=t;
					Function f=c.function;
					while(f==c.function) {f=Function.getFunction(randy.nextInt(types));}
					c.function=f;
				}
			}
		}
	}
	
//	void anneal(int iteration, float rate, Random randy) {
//		for(Cell c : list) {
//			{
//				float b=c.bias;
//				float r=randy.nextFloat();
//				float r2=r*r*(randy.nextBoolean() ? 1 : -1);
//				float b2=b+b*r2;
//				c.setBias(b2, true);
//			}
//		}
//	}
	
	private static float[][] makeFloatMatrix(int[] dims){
		//TODO: Note - these can be made 1 longer with a constant of 1 to include bias
		float[][] matrix=new float[dims.length][];
		for(int i=0; i<dims.length; i++) {
			matrix[i]=new float[dims[i]];
		}
		return matrix;
	}
	
	private static double[][] makeDoubleMatrix(int[] dims){
		//TODO: Note - these can be made 1 longer with a constant of 1 to include bias
		double[][] matrix=new double[dims.length][];
		for(int i=0; i<dims.length; i++) {
			matrix[i]=new double[dims[i]];
		}
		return matrix;
	}
	
	private static Cell[][] makeNodes(int[] dims, ArrayList<Cell> list, float[][] values, float[][] eOverNext){
		final Cell[][] net=new Cell[dims.length][];
		assert(list.isEmpty());
		list.add(null);
		int prevWidth=-1, width=dims[0], nextWidth;
		for(int layerNum=0; layerNum<dims.length; layerNum++) {
			final int type=(layerNum<dims.length-1) ? Cell.defaultActivationType : Cell.finalLayerType;
			nextWidth=(layerNum+1>=dims.length ? -1 : dims[layerNum+1]);
			Cell[] layer=net[layerNum]=new Cell[width];
			float[] lvals=values[layerNum];
			float[] eons=eOverNext[layerNum];
			assert(lvals.length==width) : layerNum+", "+lvals.length+", "+width;
			for(int i=0; i<width; i++){
				Cell c=new Cell(list.size(), type, i, layerNum, prevWidth, width, nextWidth, lvals, eons);
				layer[i]=c;
				list.add(c);
			}
			prevWidth=width;
			width=nextWidth;
		}
		return net;
	}
	
	//Fully-connected
	private static void makeEdgesSlow(Cell[][] net, ArrayList<Cell> list, Random randy, float density){
		long numEdges=0;
		for(int layerNum=1; layerNum<net.length; layerNum++) {
			Cell[] layer=net[layerNum], prev=net[layerNum-1];
			final boolean finalLayer=(layerNum==net.length-1);
//			int width=layer.length, pwidth=prev.length;
			for(Cell c : layer) {
				c.setBias(randomWeight(randy, 0.8f), true);
				final float[] weights=c.weights;
				for(int i=0; i<weights.length; i++) {
					Cell p=prev[i];
					float weight=randomWeight(randy, 0.5f);
					weight=(finalLayer || randy.nextFloat()<=density ? weight : 0);
					weights[i]=weight;
					numEdges++;
					Edge e=new Edge(numEdges, p, c, weight);
					c.inputs.add(e);
					p.outputs.add(e);
				}
//				assert(c.check());
			}
		}
	}
	
	//Fully-connected
	private static void makeEdgesFast(Cell[][] net, ArrayList<Cell> list, Random randy, float density){
		long numEdges=0;
		for(int layerNum=1; layerNum<net.length; layerNum++) {
			Cell[] layer=net[layerNum], prev=net[layerNum-1];
			final boolean finalLayer=(layerNum==net.length-1);
			for(Cell c : layer) {
				c.setBias(randomWeight(randy, 0.8f), true);//TODO: 0.8f should be a parameter
				final float[] weights=c.weights;
				if(weights!=null) {
					for(int i=0; i<weights.length; i++) {
//						Cell p=prev[i];
						float w=randomWeight(randy, 0.5f);
						weights[i]=(finalLayer || randy.nextFloat()<=density ? w : 0);
//						numEdges++;
//						Edge e=new Edge(numEdges, p, c, w);
//						c.inputs.add(e);
//						p.outputs.add(e);
					}
				}
//				assert(c.check());
			}
		}
	}
	
	//Distributed between -2 and 2.  probFlat=1 gives a uniform distribution from -1 to 1; 0 gives quadratic.
	private static float randomWeight(Random randy, float probFlat) {
		float weight=0;
//		while(weight<=0 || weight>=1){
//			weight=randy.nextFloat()/**randy.nextFloat()*/*(randy.nextBoolean() ? 1 : -1);
//			System.err.println(weight);
//		}
//		weight=randy.nextFloat()*randy.nextFloat()*2*(randy.nextBoolean() ? 1 : -1);
		if(randy.nextFloat()<=probFlat){
			weight=1-2*randy.nextFloat();
		}else{
			weight=randy.nextFloat()*randy.nextFloat()*2*(randy.nextBoolean() ? 1 : -1);
		}
		return weight;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/
	
	
	void anneal(float strength, Random randy) {
		//assert(check());
		if(CellNet.FAST) {
			for(Cell c : list) {
				if(c!=null) {c.annealFast(strength, randy);}
			}
		}else {
			for(Cell c : list) {
				if(c!=null) {c.annealSlow(strength, randy);}
			}
		}
		//assert(check());
	}
	
	public void processSample(Sample s, float alpha, boolean backProp, float weightMult) {
		//assert(check());
		
		applyInput(s.in);
		if(FAST) {
			feedForwardFast();
		}else {
			feedForwardSlow();
		}

//		s.errorValue=0;
		for(int i=0; i<finalLayer.length; i++) {
			s.result[i]=finalLayer[i].value();
//			s.errorValue+=(s.result[i]-s.goal[i]);
		}
		
		s.calcError(weightMult);
		
		if(!backProp) {return;}
		
//		double error=calcError(s.goal);
//		assert(error==s.errorMagnitude) : "\n"+s.result[0]+", "+s.goal[0]+", "+s.errorMagnitude+"\n"
//				+ finalLayer[0].value()+", "+s.goal[0]+", "+error+"\n";
		
//		s.errorMagnitude=(float)error;
//		s.weightedErrorMagnitude=(float)error*Cell.getWeightMultiplier(s.result[0], s.goal[0]);
		
//		errorSum+=s.errorMagnitude;
//		weightedErrorSum+=s.weightedErrorMagnitude;
		if(FAST) {
			backPropFast(s.goal, alpha, weightMult);
		}else {
			backPropSlow(s.goal, alpha, weightMult);
		}
		//assert(check());
	}
	
	public void applyInput(FloatList valuesIn) {
		//assert(check());
		assert(valuesIn.size==dims[0]) : valuesIn.size+"; "+Arrays.toString(dims);
		if(FAST) {
			Vector.copy(values[0], valuesIn.array);
			return;
		}
		for(int i=0; i<valuesIn.size; i++) {
			final Cell c=net[0][i];
//			//assert(c.check());
			float v=valuesIn.get(i);
//			assert(v>=0 && v<=1) : v;//Assume 0-1 input range
			assert(c.value()==values[0][i]);
//			assert(values[0]==c.values);
			net[0][i].setValue(v);
//			values[0][i]=v;
			assert(c.value()==values[0][i]);
//			//assert(c.check());
		}
		//assert(check());
	}
	
	public void applyInput(float[] valuesIn) {
		//assert(check());
		assert(valuesIn.length==dims[0]) : valuesIn.length+"; "+Arrays.toString(dims);
		if(FAST) {
			Vector.copy(values[0], valuesIn);
			return;
		}
		for(int i=0; i<valuesIn.length; i++) {
			final Cell c=net[0][i];
//			//assert(c.check());
			float v=valuesIn[i];
//			assert(v>=0 && v<=1) : v;//Assume 0-1 input range
//			assert(c.value()==values[0][i]) : v+", "+c.value()+", "+values[0][i];
//			assert(values[0]==c.values);
			net[0][i].setValue(v);
//			values[0][i]=v;
			assert(c.value()==values[0][i]);
//			//assert(c.check());
		}
		//assert(check());
	}
	
	public void feedForwardSlow(){
		//assert(check());
		for(int i=1; i<layers; i++){
			Cell[] layer=net[i];
			for(Cell c : layer) {
				double sum=c.summateSlow();
//				//assert(c.check());
				assert(values[i-1].length==c.inputs.size()) : i+", "+layer.length+", "+values[i-1].length+", "+c.inputs.size();
				c.setValue((float)Functions.sigmoid((float)(sum)));
//				//assert(c.check());
			}
		}
		//assert(check());
	}
	
	public float feedForwardFast(){
		//assert(check());
		for(int lnum=1; lnum<layers; lnum++){
			final float[] valuesIn=values[lnum-1];
			Cell[] layer=net[lnum];
			if(SPECIAL_FMA) {
				Vector.feedForward(layer, valuesIn);
			}else {
				for(int cnum=0; cnum<layer.length; cnum++) {
					Cell c=layer[cnum];
					//				//assert(c.check());
					c.summateFast(valuesIn);
					//				//assert(c.check());
				}
			}
		}
		return finalLayer[0].value();
	}
	
	public void backPropSlow(float[] truth, float alpha, float weightMult) {
		//assert(check());
		{
			for(int i=0; i<finalLayer.length; i++) {
				final Cell c=finalLayer[i];
				//assert(c.check());
				c.updateEdgesFinalLayerSlow(alpha, truth[i], weightMult);
				//assert(c.check());
			}
		}

		for(int lnum=layers-2; lnum>0; lnum--) {
			Cell[] layer=net[lnum];
			for(Cell c : layer){
				//assert(c.check());
				c.updateEdgesHiddenLayerSlow(alpha);
				//assert(c.check());
			}
		}
		//assert(check());
	}
	
	public void backPropFast(float[] truth, float alpha, float weightMult) {
		//assert(check());
		{
			final float[] valuesIn=values[values.length-2];
			for(int i=0; i<finalLayer.length; i++) {
				final Cell c=finalLayer[i];
				//assert(c.check());
				c.updateEdgesFinalLayerFast(alpha, truth[i], valuesIn, weightMult);
				//assert(c.check());
			}
		}
		
		if(TRANSPOSE_EACH_SAMPLE){
			transpose();//Allows you to ignore Edges entirely
		}
		
		for(int lnum=layers-2; lnum>0; lnum--) {
			final float[] valuesIn=values[lnum-1];
			final float[] eOverNetNext=eOverNet[lnum+1];
			Cell[] layer=net[lnum];
			float[][] weightsOutLnum=weightsOut[lnum];
			
			if(SPECIAL_FMA) {Vector.backPropFma(layer, eOverNetNext, weightsOutLnum);}
			for(int i=0; i<layer.length; i++){
				Cell c=layer[i];
				//assert(c.check());
				c.updateEdgesHiddenLayerFast(alpha, valuesIn, eOverNetNext, weightsOutLnum[i]);
				//assert(c.check());
			}
		}
		//assert(check());
	}
	
	void transpose() {
		for(int layer=0; layer<weights.length-1; layer++) {
			float[][] a=weights[layer+1], b=weightsOut[layer];
			for(int i=0; i<b.length; i++) {
				float[] d=b[i];
				for(int j=0; j<a.length; j++) {
					d[j]=a[j][i];
				}
			}
		}
	}
	
	public void applyChanges(int samples, float alpha) {
		//assert(check());
		float invSamples=1f/samples;
		assert(invSamples<=1) : samples+", "+invSamples;
		if(FAST) {
			for(Cell c : list) {
				if(c!=null) {
					//assert(c.check());
//					synchronized(c) {
						c.applyUpdatesFast(invSamples, alpha);
//					}
					//assert(c.check());
				}
			}
		}else {
			for(Cell c : list) {
				if(c!=null) {
					//assert(c.check());
//					synchronized(c) {
						c.applyUpdatesSlow(invSamples, alpha);
//					}
					//assert(c.check());
				}
			}
		}
		//assert(check());
		epochsTrained++;
		samplesTrained+=samples;
	}
	
	public float getOutput(int outnum){
		return net[layers-1][outnum].value();
	}
	
	public float[] getOutput(){
		Cell[] outLayer=net[layers-1];
		float[] out=new float[outLayer.length];
		for(int i=0; i<outLayer.length; i++) {
			out[i]=outLayer[i].value();
		}
		return out;
	}
	
	@Deprecated
	public double calcError(float[] truth){
		double error=0;
		for(int i=0; i<finalLayer.length; i++){
			Cell c=finalLayer[i];
			float t=truth[i];
			float e=c.calcError(t);
			error+=e;
		}
		return error;
	}

	
	/*--------------------------------------------------------------*/
	/*----------------            X            ----------------*/
	/*--------------------------------------------------------------*/
	
	public String toString() {
		StringBuilder sb=new StringBuilder();
		for(int layernum=0; layernum<layers; layernum++) {
			Cell[] layer=net[layernum];
			sb.append("\n* Layer "+layernum+", nodes="+layer.length+" *");
			for(Cell c : layer) {
				sb.append("\n"+c.toString());
			}
		}
		return sb.toString();
	}
	
	public ByteBuilder header() {
		ByteBuilder bb=new ByteBuilder();
		bb.append("##bbnet").nl();
		bb.append("#version ").append(version).nl();
		if(CONCISE) {bb.append("#concise").nl();}
		bb.append("#seed ").append(seed).nl();
		bb.append("#annealseed ").append(annealSeed).nl();
		bb.append("#layers ").append(layers).nl();
		if(epochsTrained>0) {bb.append("#epochs ").append(epochsTrained).nl();}
		if(samplesTrained>0) {bb.append("#samples ").append(samplesTrained).nl();}
		
		bb.append("#dims");
		for(int d : dims) {bb.space().append(d);}
		bb.nl();
		for(String s : commands) {bb.append(s).nl();}
		if(lastStats!=null) {bb.append("##stats ").append(lastStats).nl();}
		bb.append("##fpr ").append(fpRate, 6).nl();
		bb.append("##fnr ").append(fnRate, 6).nl();
		bb.append("##err ").append(errorRate, 6).nl();
		bb.append("##wer ").append(weightedErrorRate, 6).nl();
		bb.append("##ctf ").append(cutoff, 6).nl();
		return bb;
	}
	
	public ByteBuilder toBytes() {
		
		ByteBuilder bb=header();
		bb.append("#edges").nl();
		
		lastLinesWritten=6;
		long edgeCount=0;
		for(int lnum=1; lnum<layers; lnum++) {
			Cell[] layer=net[lnum];
			Cell[] prev=net[lnum-1];
			bb.append("##layer ").append(lnum).nl();
			for(Cell c : layer){
				if(CONCISE) {
					lastLinesWritten++;
					bb.append('C').append(c.id()).space().append(c.typeString());
					bb.space().append(c.bias(), 6, true);
					for(int i=0; i<c.weights.length; i++) {bb.space().append(c.weights[i], 6, true);}
					bb.nl();
					edgeCount+=c.weights.length;
				}else {
					lastLinesWritten+=2+c.weights.length;
					bb.append("C").append(c.id()).space().append(c.typeString()).nl();
					bb.append(0).space().append(c.id()).space().append(
							String.format(Locale.ROOT, "%.6f", c.bias())).space().append('b').append(c.id()).nl();
					if(FAST) {
						for(int i=0; i<c.weights.length; i++) {
							bb.append(prev[i].id()).space().append(c.id()).space().append(
									String.format(Locale.ROOT, "%.6f", c.weights[i])
									).space().append('w').append(edgeCount).nl();
							edgeCount++;
						}
					}else {
						for(Edge e : c.inputs) {
							bb.append(e.source.id()).space().append(c.id()).space().append(
									String.format(Locale.ROOT, "%.6f", e.weight())).space().append('w').append(e.id()).nl();
						}
					}
				}
			}
		}
		return bb;
	}
	
	public synchronized CellNet copy(boolean copyWeight2) {
		//assert(check());
		CellNet copy=new CellNet(dims, seed, commands);
		ArrayList<Cell> list2=copy.list;
		for(int i=1; i<list.size(); i++) {
			Cell c=list.get(i);
			//assert(c.check());
			Cell c2=list2.get(i);
			if(c.inputs!=null) {
				assert(!FAST);
				for(Edge e : c.inputs) {
					int sid=e.source.id();
					Cell source2=list2.get(sid);
					Edge e2=new Edge(e.id(), source2, c2, e.weight());
					c2.inputs.add(e2);
					source2.outputs.add(e2);
				}
			}
			c2.setFrom(c, copyWeight2);
		}
		copy.commands=commands;
		copy.errorRate=errorRate;
		copy.weightedErrorRate=weightedErrorRate;
		copy.fpRate=fpRate;
		copy.fnRate=fnRate;
		copy.tpRate=tpRate;
		copy.tnRate=tnRate;
		copy.cutoff=cutoff;
		copy.alpha=alpha;
		copy.annealStrength=annealStrength;
		copy.annealSeed=annealSeed;
		copy.epoch=epoch;
		copy.epochsTrained=epochsTrained;
		copy.lastStats=lastStats;
		copy.samplesTrained=samplesTrained;
		
//		assert(copy.check());
		return copy;
	}
	
	public CellNet setFrom(CellNet cn, boolean copyWeight2) {
		//assert(cn.check());
		//assert(check());
		ArrayList<Cell> list2=cn.list;
		for(int i=1; i<list.size(); i++) {
			Cell c=list.get(i);
			//assert(c.check());
			Cell c2=list2.get(i);
			c.setFrom(c2, copyWeight2);
//			if(c.inputs!=null) {
//				assert(c.inputs.size()==c.weights.length);
//				assert(c2.inputs.size()==c2.weights.length);
//				assert(c.inputs.size()==c2.inputs.size());
//				for(int j=0; j<c.inputs.size(); j++) {
//					Edge e=c.inputs.get(j);
//					Edge e2=c2.inputs.get(j);
//					e.setFrom(e2);
//					c.weights[j]=c2.weights[j];
//				}
//			}
		}
		//assert(cn.check());
		//assert(check());
		commands=cn.commands;
		errorRate=cn.errorRate;
		weightedErrorRate=cn.weightedErrorRate;
		fpRate=cn.fpRate;
		fnRate=cn.fnRate;
		tpRate=cn.tpRate;
		tnRate=cn.tnRate;
		cutoff=cn.cutoff;
		alpha=cn.alpha;
		annealStrength=cn.annealStrength;
		annealSeed=cn.annealSeed;
		epoch=cn.epoch;
		epochsTrained=cn.epochsTrained;
		lastStats=cn.lastStats;
		samplesTrained=cn.samplesTrained;
		return this;
	}
	
	public boolean check() {
		for(int i=1; i<list.size(); i++) {
			Cell c=list.get(i);
			boolean b=c.check();
			if(!b) {return false;}
		}
		return true;
	}
	
	public void addError(CellNet net2) {
//		errorSum+=net2.errorSum;
//		weightedErrorSum+=net2.weightedErrorSum;
		for(int i=1; i<list.size(); i++) {
			Cell c=list.get(i);
			Cell c2=net2.list.get(i);
			c.addError(c2);
		}
	}
	
	public void accumulate(CellNet net2) {
		synchronized(this) {
			synchronized(net2) {
				if(FAST) {accumulateFast(net2);}
				else {accumulateSlow(net2);}
			}
		}
	}
	
	private void accumulateSlow(CellNet net2) {
//		errorSum+=net2.errorSum;
//		weightedErrorSum+=net2.weightedErrorSum;
		for(int i=1; i<list.size(); i++) {
			Cell c=list.get(i);
			Cell c2=net2.list.get(i);
//			synchronized(c) {
				c.accumulateSlow(c2);
//			}
		}
	}
	
	private void accumulateFast(CellNet net2) {
//		errorSum+=net2.errorSum;
//		weightedErrorSum+=net2.weightedErrorSum;
		for(int i=1; i<list.size(); i++) {
			Cell c=list.get(i);
			Cell c2=net2.list.get(i);
			c.accumulateFast(c2);
		}
	}
	
	synchronized public void clear() {
//		errorSum=0;
//		weightedErrorSum=0;
		if(FAST) {
			for(Cell c : list) {
				if(c!=null) {c.clearTempFast();}
			}
		}else {
			for(Cell c : list) {
				if(c!=null) {c.clearTempSlow();}
			}
		}
	}
	
	public long countEdges() {
		long sum=0;
		for(float[][] y : weights) {
			if(y!=null) {
				for(float[] z : y) {
					if(z!=null) {
						sum+=z.length;
					}
				}
			}
		}
		return sum;
	}
	
	@Override
	public int compareTo(CellNet b) {
		if(b==null) {
			return 1;
		}else if(compareCode==compareWER && weightedErrorRate!=b.weightedErrorRate) {
			return weightedErrorRate<b.weightedErrorRate ? 1 : -1;
		}else if(compareCode==compareERR && errorRate!=b.errorRate) {
			return errorRate<b.errorRate ? 1 : -1;
		}else if(compareCode==compareFNR && fnRate!=b.fnRate) {
			return fnRate<b.fnRate ? 1 : -1;
		}else if(compareCode==compareFPR && fpRate!=b.fpRate) {
			return fpRate<b.fpRate ? 1 : -1;
		}
		
		if(weightedErrorRate!=b.weightedErrorRate) {
			return weightedErrorRate<b.weightedErrorRate ? 1 : -1;
		}else if(errorRate!=b.errorRate) {
			return errorRate<b.errorRate ? 1 : -1;
		}else if(fnRate!=b.fnRate) {
			return fnRate<b.fnRate ? 1 : -1;
		}else if(fpRate!=b.fpRate) {
			return fpRate<b.fpRate ? 1 : -1;
		}
		
		return 0;
	}
	
	final float pivot() {
		return compareCode==compareWER ? weightedErrorRate : 
			compareCode==compareERR ? errorRate : 
			compareCode==compareFNR ? fnRate :
			compareCode==compareFPR ? fpRate : weightedErrorRate;
	}
	
	/*--------------------------------------------------------------*/
	
	public int numLayers(){
		return dims.length;
	}
	
	public int numInputs(){
		return dims[0];
	}
	
	public int numOutputs(){
		return dims[dims.length-1];
	}
	
	/*--------------------------------------------------------------*/
	
	@Override
	public boolean equals(Object b) {
		return epoch==((CellNet)b).epoch;
	}

	@Override
	public int hashCode() {return epoch;}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	float errorRate=999, weightedErrorRate=999;
	float fpRate=999, fnRate=999, tpRate=-999, tnRate=-999;
	float alpha=-1, annealStrength=-1;
	public float cutoff=-1;
	int epoch=-1;
	int count=1; //for printing
	
	long epochsTrained=0;
	long samplesTrained=0;
	String lastStats=null;
	
	/*--------------------------------------------------------------*/
	
	final long seed;
	long annealSeed=-1;
	final int layers;
	final int[] dims; //Stores widths of layers
	final Cell[][] net;
	final Cell[] finalLayer;
	final ArrayList<Cell> list;
	
	final float[][] values;
	final float[][] eOverNet;
	final float[][][] weights;
	final float[][][] weightsOut;
	
	public ArrayList<String> commands;
	
	long lastLinesWritten=0;
	
	public static final boolean TRANSPOSE_EACH_SAMPLE=false;
	public static final boolean SPECIAL_FMA=true;//~20% faster when true
	public static boolean CONCISE=true;
	public static boolean FAST=true;
	public static boolean verbose=false;
	public static final int version=1;
	
	static int compareCode=0;
	final static int compareWER=0, compareERR=1, compareFNR=2, compareFPR=3;
	
	/*--------------------------------------------------------------*/
	
}
