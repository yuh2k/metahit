package ml;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Random;

import shared.Tools;

public class CellDouble extends Source {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	public CellDouble(int id_, int activationType, int lpos_, int layer_, int pWid, int wid, int nWid,
			float[] values_, double[] eOverNetArray_) {
		id=id_;
//		type=activationType;
		function=Function.getFunction(activationType);
		lpos=lpos_;
		layer=layer_;
		inputs=(pWid<1 || CellNet.FAST ? null : new ArrayList<Edge>(pWid));
		outputs=(nWid<1 || CellNet.FAST ? null : new ArrayList<Edge>(nWid));
		values=values_;
		eOverNetArray=eOverNetArray_;
		assert(values.length==wid);
		weights=(pWid<1 ? null : new float[pWid]);
//		weights2=(pWid<1 ? null : new double[pWid]);
		deltas=(pWid<1 ? null : new double[pWid]);
		assert(lpos>=0 & lpos<wid);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Non-Mutators         ----------------*/
	/*--------------------------------------------------------------*/
	
	public double summateSlow() {
		sum=bias;
//		CellNet.println("sum="+sum);
		for(Edge e : inputs) {
			sum+=e.source.value()*e.weight();
//			CellNet.println("sum="+sum);
		}
//		CellNet.println("sig="+activation(sum));
//		CellNet.println("sig -1="+activation(-1));
//		CellNet.println("sig 1="+activation(1));
//		assert(bias==0) : bias+", "+sum+", "+this;
		return sum;
	}
	
	public void summateFast(float[] valuesIn) {
		sum=bias;
//		CellNet.println("sum="+sum);
//		assert(valuesIn.length==inputs.size()) : valuesIn.length+", "+inputs.size();
		assert(valuesIn.length==weights.length) : valuesIn.length+", "+weights.length;
		for(int i=0; i<valuesIn.length; i++) {//TODO: Should this run from 1 to length?
//			final Edge e=inputs.get(i);//TODO: Remove
//			assert(valuesIn[i]==e.source.value());
//			assert(weights[i]==e.weight()) : weights[i]+", "+e.weight();
			sum+=valuesIn[i]*weights[i];
		}
		final float v=(float)activation(sum);
		setValue(v);
	}
	
	public float calcError(float ideal) {
		float e=ideal-value();
		return 0.5f*e*e;
	}
	
	synchronized public boolean check() {
		assert(false);
		if(value!=values[lpos] && (value!=-1 && values[lpos]!=0)) {
			assert(false) : id+", "+layer+", "+lpos+": "+value+", "+values[lpos]+", "+Arrays.toString(values);
			return false;
		}
		if(eOverNet!=eOverNetArray[lpos]) {
			assert(false) : id+", "+layer+", "+lpos+": "+eOverNet+", "+eOverNetArray[lpos]+", "+Arrays.toString(values);
			return false;
		}
		if(layer>0) {
			assert(weights!=null);
			assert(inputs!=null);
			assert(weights.length==inputs.size()) : weights.length+", "+inputs.size();
//			assert(weights2.length==weights.length);
			assert(deltas.length==weights.length);
			for(int i=0; i<weights.length; i++) {
				Edge e=inputs.get(i);
				if(e.weight()!=weights[i]) {
					assert(false) : id+", "+layer+", "+lpos+": "+e.weight()+", "+weights[i];
					return false;
				}
//				if(e.weight2!=weights2[i]) {
//					assert(false) : id+", "+layer+", "+lpos+": "+e.weight2+", "+weights2[i];
//					return false;
//				}
				if(e.delta!=deltas[i]) {
					assert(false) : id+", "+layer+", "+lpos+": "+e.delta+", "+deltas[i];
					return false;
				}
			}
//			for(int i=0; i<weights2.length; i++) {
//				Edge e=outputs.get(i);
//				if(e.weight2!=weights2[i]) {
//					assert(false) : id+", "+layer+", "+lpos+": "+e.weight2+", "+weights2[i];
//					return false;
//				}
//			}
		}else {
//			assert(inputs==null || outputs.isEmpty())  : id+", "+layer+", "+lpos+": "+weights2+": "+outputs+"\n"+toString()+"\n";
			assert(inputs==null || outputs.isEmpty())  : id+", "+layer+", "+lpos+": "+deltas+": "+outputs+"\n"+toString()+"\n";
		}
		return true;
	}

	/*--------------------------------------------------------------*/
	/*----------------           Mutators           ----------------*/
	/*--------------------------------------------------------------*/
	
	void applyUpdatesSlow(float invSamples, float alpha) {
		if(layer<1) {return;}
		adjustBias(invSamples, alpha);
		for(int i=0; i<weights.length; i++) {
			Edge e=inputs.get(i);
			assert(e.weight()==weights[i]);
//			assert(e.weight2==weights2[i]);
			assert(e.delta==deltas[i]);
//			weights[i]=e.adjustWeight(invSamples);
//			final float w3=(float)(weights2[i]*invSamples);
			final float d3=(float)(deltas[i]*invSamples*alpha);
			final float w4=weights[i]+d3;
//			weights2[i]=e.weight2=0;
			deltas[i]=e.delta=0;
			weights[i]=w4;
			e.setWeight(w4);
		}
//		for(Edge e : inputs) {
//			e.adjustWeight(invSamples);
//		}
//		loadFast();
	}
	
	void applyUpdatesFast(float invSamples, float alpha) {
		if(layer<1) {return;}
		adjustBias(invSamples, alpha);
		for(int i=0; i<weights.length; i++) {
			final float w=weights[i];
			
//			final float w3=(float)(weights2[i]*invSamples);
//			weights2[i]=0;
			
			final float d3=(float)(deltas[i]*invSamples*alpha);
			final float w4=w+d3;
//			assert(Tools.absdif(w3, w4)<0.00001) : "w3="+w3+", w4="+w4+", d3="+d3+", w="+weights[i]+", d"+deltas[i]+", w2="+weights2[i];
			deltas[i]=0;

			weights[i]=w4;
		}
	}
	
	public void addError(float e) {
		assert(error>=0);
		error+=e;
	}
	
	public void clearError() {
		error=0;
	}
	
	public void clearTempSlow() {
//		bias2=0;
		biasDelta=0;
		error=0;
		if(layer<1) {return;}
		for(int i=0; i<weights.length; i++) {
			Edge e=inputs.get(i);
			e.clearTemp();
		}
//		Arrays.fill(weights2, 0);
		Arrays.fill(deltas, 0);
	}
	
	public void clearTempFast() {
//		bias2=0;
		biasDelta=0;
		error=0;
		if(layer<1) {return;}
//		Arrays.fill(weights2, 0);
		Arrays.fill(deltas, 0);
	}
	
	public void setBias(float b2, boolean ignoreAssertion) {
		assert(layer>0);
//		float dif=Tools.absdif(bias, b2);
//		assert(ignoreAssertion || dif<0.1f || dif<0.2f*bias) : dif+", "+bias+", "+b2; //just checking for bugs.
//		if(!(ignoreAssertion || dif<0.1f || dif<0.2f*bias)) {System.out.print("*");}
//		if(!ignoreAssertion) {System.err.println(bias);}
		bias=b2;
	}
	
	private void adjustBias(float invSamples, float alpha) {
		assert(layer>0);
		assert(alpha>0) : alpha;
//		System.err.println(String.format("b"+id()+": %.5f -> %.5f * %.2f = %.5f",
//				bias, bias2, invSamples, bias2*invSamples));
//		float b=((float)bias2)*invSamples;
		float bFromDelta=(float)(bias+biasDelta*invSamples*alpha*biasAlphaMult);
//		assert(Tools.absdif(b, bFromDelta)<0.00001) : "b="+b+", bFD="+bFromDelta+", bias="+bias+", bias2="+bias2+", bD="+biasDelta+", invS="+invSamples+", a="+alpha+", bAM="+biasAlphaMult
//			+", layer="+layer+" cid="+id;
		setBias(bFromDelta, false);
//		bias2=0;
		biasDelta=0;
	}
	
	public void addError(CellDouble c2) {
		error+=c2.error;
	}
	
	@Override
	void setValue(float v) {
		assert(layer==0 || v==(float)activation(sum)) : v+", "+activation(sum)+", "+sum;
//		assert(value==values[lpos]);
		values[lpos]=value=v;
	}
	
	/*--------------------------------------------------------------*/
	
	//TODO: Replace with array of function subclasses
	final double activation(double x) {
		return function.activate(x);
//		double y=(type==Function.SIGMOID ? Functions.sigmoid(x) : 
//			type==Function.TANH ? Functions.tanh(x) :
//				type==Function.RSLOG ? Functions.rslog(x) :
//					type==Function.MSIG ? Functions.mSig(x) :
//						type==Function.SWISH ? Functions.swish(x) :
//							Double.NaN);
//		assert(!Double.isNaN(y)) : x+", "+y+", "+type;
//		return y;
	}
	
//	private final double derivativeX(double x) {
//		double d=(type==SIGMOID ? Functions.sigmoidDerivativeX(x) : 
//			type==TANH ? Functions.tanhDerivativeX(x) :
//				type==RSLOG ? Functions.rslogDerivativeX(x) :
//				Functions.swishDerivativeX(x));
//		assert(!Double.isNaN(d)) : x+", "+d;
//		return d;
//	}
	
	private final double derivativeXFX(double x, double fx) {
		double d=function.derivativeXFX(x, fx);
//		double d=(type==Function.SIGMOID ? Functions.sigmoidDerivativeFX(fx) : 
//			type==Function.TANH ? Functions.tanhDerivativeFX(fx) :
//				type==Function.RSLOG ? Functions.rslogDerivativeX(x) :
//					type==Function.MSIG ? Functions.mSigDerivativeXFX(x, fx) :
//						type==Function.SWISH ? Functions.swishDerivativeX(x) :
//							Double.NaN);
//		assert(Tools.absdif(d, derivativeX(x))<0.0001) : "\nx="+x+", sum="+sum+"\n"+
//				"d="+d+", dx="+derivativeX(x)+"\n"+
//				"fx="+fx+", value="+value+", sig(x)="+Functions.sigmoid(x)+"\n"+
//				"sigdX(x)="+Functions.sigmoidDerivativeX(x)+
//				", sigdFX(fx)="+Functions.sigmoidDerivativeFX(fx)+
//				"\nswishdX(x)="+Functions.swishDerivativeX(x)+
//				", swishdFX(x, fx)="+Functions.swishDerivativeXFX(x, fx)+
//				"\n";
		assert(!Double.isNaN(d)) : x+", "+fx+", "+d;
		return d;
	}
	
	/*--------------------------------------------------------------*/
	
	void updateEdgesFinalLayerSlow(float alpha, float target, float weightMult) {
		assert(alpha>0) : alpha;
		//assert(check());
		final float v=value();
		eTotalOverOut=calcETotalOverOut(v, target, weightMult);
		outOverNet=(float)derivativeXFX(sum, v);
		eOverNet=eOverNetArray[lpos]=eTotalOverOut*outOverNet;
		
		for(int i=0; i<weights.length; i++) {
			final Edge e=inputs.get(i);
//			assert(e.source.check());
			float netOverWeight=e.source.value();
			double eTotalOverWeight=eOverNet*netOverWeight;
			
//			double incr=e.weight()-alpha*eTotalOverWeight;
//			e.incrementWeight2(incr);
//			weights2[i]+=incr;
//			assert(weights2[i]==e.weight2);
			
			double dincr=-eTotalOverWeight;
			e.incrementDelta(dincr);
			deltas[i]+=dincr;
			assert(deltas[i]==e.delta);
		}
		{
//			float product=(float) (eTotalOverOut*outOverNet);
//			bias2+=bias-alpha*eOverNet*biasAlphaMult;
			biasDelta-=eOverNet;
//			CellNet.println("C"+id()+": bias="+bias+", bias2="+bias2);
		}
		//assert(check());
	}
	
	void updateEdgesFinalLayerFast(float alpha, float target, float[] valuesIn, float weightMult) {
		assert(alpha>0) : alpha;
		//assert(check());
		final float v=value();
		eTotalOverOut=calcETotalOverOut(v, target, weightMult);
		outOverNet=(float)derivativeXFX(sum, v);
//		final double eTotalOverOut_X_outOverNet=eTotalOverOut*outOverNet;
		eOverNet=eOverNetArray[lpos]=eTotalOverOut*outOverNet;
		
		for(int i=0; i<weights.length; i++) {
//			assert(e.source.check());
			final float netOverWeight=valuesIn[i];//e.source.value();
//			assert(valuesIn[i]==e.source.value) : valuesIn[i]+", "+e.source.value+", "+valuesIn.length+", "+inputs.size();
			final double eTotalOverWeight=eOverNet*netOverWeight;
			
//			final double w=weights[i];
//			final double incr=w-alpha*eTotalOverWeight;
//			weights2[i]+=incr;

			deltas[i]-=eTotalOverWeight;
			//TODO: Alpha is not really needed here;
			//the result could be accumulated and multiplied by alpha at the end
		}
		{
			assert(layer>0);
//			bias2+=bias-alpha*eOverNet*biasAlphaMult;
			biasDelta-=eOverNet;
		}
		//assert(check());
	}
	
	/*--------------------------------------------------------------*/
	
	void updateEdgesHiddenLayerSlow(float alpha) {
		assert(alpha>0) : alpha;
		//assert(check());
		eTotalOverOut=0;
		outOverNet=value()*(1-value());
		
		for(int i=0; i<outputs.size(); i++) {
			final Edge o=outputs.get(i);
//			assert(o.dest.check());
//			double eOverNetDest=o.dest.eTotalOverOut*o.dest.outOverNet;
//			assert(eOverNetDest==o.dest.eOverNet);
			double eOverNetDest=o.dest.eOverNet;
			float netOverOut=o.weight();
			double eOverOut=eOverNetDest*netOverOut;
			eTotalOverOut+=eOverOut;
		}
		eOverNet=eOverNetArray[lpos]=eTotalOverOut*outOverNet;
		
		for(int i=0; i<inputs.size(); i++) {
			final Edge e=inputs.get(i);
//			assert(e.source.check());
			float netOverWeight=e.source.value();
			double eTotalOverWeight=(eTotalOverOut*outOverNet*netOverWeight);
//			e.incrementWeight2(e.weight()-alpha*eTotalOverWeight);
			
//			final double w=e.weight();
//			final double incr=w-alpha*eTotalOverWeight;
//			weights2[i]+=incr;
//			e.incrementWeight2(incr);
//			assert(weights2[i]==e.weight2);
			
			e.incrementDelta(-eTotalOverWeight);
			assert(deltas[i]==e.delta);
		}
		if(layer>0){
			float product=(float) (eTotalOverOut*outOverNet);
//			bias2+=bias-alpha*product*biasAlphaMult;//Bias adjusts slower than edges
			biasDelta-=product;
		}
		//assert(check());
	}
	
	void updateEdgesHiddenLayerFast(float alpha, float[] valuesPrev, double[] eOverNetNext, float[] weightsOut) {
		assert(alpha>0) : alpha;
		//assert(check());
		eTotalOverOut=0;
		final float v=value();
		assert(v==values[lpos]) : v+", "+values[lpos]+//Also fires on NaN, but that shouldn't happen...
			"\n"+Arrays.toString(values)+
			"\n"+Arrays.toString(weights)+"\n";
		outOverNet=(float)derivativeXFX(sum, v);
		
		for(int dest=0; dest<eOverNetNext.length; dest++){
//			Edge o=outputs.get(dest);
//			final float netOverOut=o.weight();//Random access
			
			final float netOverOut=weightsOut[dest];
			
//			assert(o.dest.eTotalOverOut==eTotalOverOutNext[dest]);
//			assert(o.dest.outOverNet==outOverNetNext[dest]);
//			final double eOverNetDest=eTotalOverOutNext[dest]*outOverNetNext[dest];//TODO: Should be already calculated by dest
//			assert(eOverNetDest==o.dest.eOverNet);
			final double eOverNetDest=eOverNetNext[dest];
//			assert(eOverNetDest==o.dest.eOverNet);
			final double eOverOut=eOverNetDest*netOverOut;
			eTotalOverOut+=eOverOut;
		}
		eOverNet=eOverNetArray[lpos]=eTotalOverOut*outOverNet;
		
//		final double eTotalOverOut_X_outOverNet=eTotalOverOut*outOverNet;
		for(int source=0; source<valuesPrev.length; source++) {
			final float netOverWeight=valuesPrev[source];
			final double eTotalOverWeight=eOverNet*netOverWeight;
			
//			final double incr=weights[source]-alpha*eTotalOverWeight;
//			weights2[source]+=incr;

			deltas[source]-=eTotalOverWeight;
		}
//		if(layer>0){
		{
			assert(layer>0) : layer;
//			bias2+=bias-alpha*eOverNet*biasAlphaMult;//Bias adjusts slower than edges
			biasDelta-=eOverNet;
		}
		//assert(check());
	}
	
	/*--------------------------------------------------------------*/
	
	public void accumulateSlow(CellDouble c2) {
		//assert(check());
		//assert(c2.check());
//		System.err.println(">>>"+c2);
//		System.err.println(">>>"+this);
		error+=c2.error;
//		bias2+=c2.bias2;
		biasDelta+=c2.biasDelta;
//		if(layer>0) {
//			final double[] c2w2=c2.weights2;
//			for(int i=0; i<weights2.length; i++) {
//				Edge e=inputs.get(i);
//				Edge e2=c2.inputs.get(i);
//				assert(e.weight2==weights2[i]);
//				assert(e.weight()==weights[i]);
//				assert(e2.weight2==c2.weights2[i]) : e2.weight2+", "+c2.weights2[i];
//				assert(e2.weight()==c2.weights[i]);
//				e.accumulate(e2);
//				weights2[i]+=c2w2[i];
//			}
//		}
		if(layer>0) {
			final double[] c2w2=c2.deltas;
			for(int i=0; i<deltas.length; i++) {
				Edge e=inputs.get(i);
				Edge e2=c2.inputs.get(i);
				assert(e.delta==deltas[i]);
				assert(e.weight()==weights[i]);
				assert(e2.delta==c2.deltas[i]) : e2.delta+", "+c2.deltas[i];
				assert(e2.weight()==c2.weights[i]);
				e.accumulate(e2);
				deltas[i]+=c2w2[i];
			}
		}
//		System.err.println(">>>"+this);
	}
	
	public void accumulateFast(CellDouble c2) {
		//assert(check());
		//assert(c2.check());
		error+=c2.error;
//		bias2+=c2.bias2;
		biasDelta+=c2.biasDelta;
//		if(layer>0) {
//			final double[] c2w2=c2.weights2;
//			for(int i=0; i<weights2.length; i++) {
//				weights2[i]+=c2w2[i];
//			}
//		}
		if(layer>0) {
			final double[] c2d=c2.deltas;
			for(int i=0; i<deltas.length; i++) {
				deltas[i]+=c2d[i];
			}
		}
	}
	
	/*--------------------------------------------------------------*/
	
	public void annealSlow(float strength, Random randy) {
		if(layer<1){return;}
		if(annealBias) {
			final float abs=Math.abs(bias);
			float xb=strength*(randy.nextFloat()-randy.nextFloat())*biasAnnealMult;
			if(abs<0.2f) {xb=xb*(Tools.max(abs*5, 0.2f));} //Weaker anneal for weaker bias //TODO: Make lower limit even weaker
			setBias(bias+xb, true);
		}
		for(int i=0; i<weights.length; i++) {
			final Edge e=inputs.get(i);
			final float w=weights[i];
			assert(w==e.weight());
			final float abs=Math.abs(w);
			float xe=strength*(randy.nextFloat()-randy.nextFloat());//*Tools.max(e.weight(), 0.01f);
//			if(abs<0.2f) {xe=xe*(Tools.max(abs*5, 0.2f));} //Weaker anneal for weaker weight
			if(abs<lowWeightAnnealCutoff) {xe=xe*(Tools.max(abs*lowWeightAnnealMult, lowWeightAnnealCutoff));} //Weaker anneal for weaker weight
			final float w2=w+xe;
			e.setWeight(w2);
			weights[i]=w2;
			assert(w2==e.weight());
		}
	}
	
	public void annealFast(float strength, Random randy) {
		if(layer<1){return;}
		if(annealBias) {
			final float abs=Math.abs(bias);
			float xb=strength*(randy.nextFloat()-randy.nextFloat())*biasAnnealMult;
			if(abs<0.2f) {xb=xb*(Tools.max(abs*5, 0.2f));} //Weaker anneal for weaker bias //TODO: Make lower limit even weaker
			setBias(bias+xb, true);
		}
		for(int i=0; i<weights.length; i++) {
//			final Edge e=inputs.get(i);
			final float w=weights[i];
//			assert(w==e.weight());
			final float abs=Math.abs(w);
			float xe=strength*(randy.nextFloat()-randy.nextFloat());//*Tools.max(e.weight(), 0.01f);
//			if(abs<0.2f) {xe=xe*(Tools.max(abs*5, 0.2f));} //Weaker anneal for weaker weight
			if(abs<lowWeightAnnealCutoff) {xe=xe*(Tools.max(abs*lowWeightAnnealMult, lowWeightAnnealCutoff));} //Weaker anneal for weaker weight
			final float w2=w+xe;
//			e.setWeight(w2);
			weights[i]=w2;
//			assert(w2==e.weight());
		}
	}
	
	public void setFrom(CellDouble c, boolean copyDelta) {
		eTotalOverOut=c.eTotalOverOut;
		outOverNet=c.outOverNet;
		bias=c.bias;
//		type=c.type;
		function=c.function;
		value=c.value;
		sum=c.sum;
		error=c.error;
		if(copyDelta) {
//			bias2=c.bias2;
			biasDelta=c.biasDelta;
		}
		if(layer==0){return;}
		if(CellNet.FAST) {
			for(int i=0; i<weights.length; i++) {
				weights[i]=c.weights[i];
			}
			if(copyDelta) {
//				for(int i=0; i<weights2.length; i++) {
//					weights2[i]=c.weights2[i];
//				}
				for(int i=0; i<deltas.length; i++) {
					deltas[i]=c.deltas[i];
				}
			}
		}else {
			if(copyDelta) {
//				for(int i=0; i<weights2.length; i++) {
//					weights2[i]=c.weights2[i];
//				}
				for(int i=0; i<deltas.length; i++) {
					deltas[i]=c.deltas[i];
				}
			}
			for(int i=0; i<weights.length; i++) {
				weights[i]=c.weights[i];
				inputs.get(i).setFrom(c.inputs.get(i), copyDelta);
			}
		}
		//assert(check());
		//assert(c.check());
		assert(CellNet.FAST || inputs.size()==0 || inputs.size()==c.inputs.size()) : "\n"+inputs+"\n"+c.inputs;
		assert(c.weights.length==weights.length) : "\n"+inputs+"\n"+c.inputs;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------          Overrides           ----------------*/
	/*--------------------------------------------------------------*/
	
	@Override
	public String toString() {
		StringBuilder sb=new StringBuilder();
//		sb.append("C"+id()+": v="+String.format("%.4f, b=%.4f, b2=%.4f, e=%.5f", value(), bias, bias2, error));
		sb.append("C"+id()+": v="+String.format("%.4f, b=%.4f, b2=%.4f, e=%.5f", value(), bias, biasDelta, error));
		if(!terminal()) {
			sb.append(", Edges: {");
			if(CellNet.FAST) {
				int prevBase=id-lpos-weights.length;
				for(int i=0; i<weights.length; i++) {
					//TODO: Could use source array or source value array here
//					return String.format("C"+(prevBase+i)+"->C"+id+",%.3f,w?=%.4f,w2=%.4f,d=%.4f; ", -1, weights[i], weights2[i], deltas[i]);
					return String.format("C"+(prevBase+i)+"->C"+id+",%.3f,w?=%.4f,d=%.4f; ", -1, weights[i], deltas[i]);
				}
			}else {
				for(Edge e : inputs) {
					sb.append(e.toString());
				}
			}
			sb.setLength(sb.length()-2);
			sb.append("}");
		}
		return sb.toString();
	}

	@Override
	public boolean terminal() {
//		return inputs==null || inputs.isEmpty();
		return layer==0;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------           Getters            ----------------*/
	/*--------------------------------------------------------------*/

	public int id() {return id;}
	public float bias() {return bias;}
	
	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	public static float toWeightedError(double rawError, float v, float target, float weightMult) {
		if(weightMult==0) {return (float)rawError;}
//		assert(rawError>=0) : rawError;
		float mult=toErrorMult(v, target, weightMult);
		double incr=toErrorIncr(rawError, v, target);
		assert(incr==0 || incr>=0 == rawError>=0) : incr+", "+rawError;
		return (float)((rawError+incr)*mult);
	}
	
	public static double toErrorIncr(double rawError, float v, float target) {
		double incr=0;
		if(target>cutoffForTraining) {
			if(v<=cutoffForTraining+spread) {incr=fnErrorIncr;}
		}else if(target<cutoffForTraining) {
			if(v>=cutoffForTraining-spread) {incr=fpErrorIncr;}
		}
		double ret=incr*(rawError>=0 ? 1 : -1);
		return ret;
	}
	
//	public static float toErrorMult(float v, float target) {
//		if(true) {return toErrorMult(v, target, 1);}
//		final float mult;
////		assert(cutoff-spread==lowThresh);
////		assert(cutoff+spread==highThresh);
////		if(v>target) {
////			mult=positiveErrorMult;
////			if(v>lowThresh && target<lowThresh) {
////				mult*=falsePositiveErrorMult;
////			}
////		}else{
////			mult=negativeErrorMult;
////			if(v<highThresh && target>=highThresh) {
////				mult*=falseNegativeErrorMult;
////			}
////		}
//		if(v>target) {
//			if(v>cutoff-spread && target<cutoff+spread) {
//				mult=falsePositiveErrorMult;
//			}else{
//				mult=positiveErrorMult;
//			}
//		}else{
//			if(v<cutoff+spread && target>=cutoff-spread) {
//				mult=falseNegativeErrorMult;
//			}else{
//				mult=negativeErrorMult;
//			}
//		}
//		return mult;
//	}
	
	public static float toErrorMult(float v, float target, float multFraction) {
		final float mult;
		if(v>target) {
			if(v>cutoffForTraining-spread && target<cutoffForTraining+spread) {
				mult=falsePositiveErrorMult;
			}else{
				mult=positiveErrorMult;
			}
		}else{
			if(v<cutoffForTraining+spread && target>=cutoffForTraining-spread) {
				mult=falseNegativeErrorMult;
			}else{
				mult=negativeErrorMult;
			}
		}
		return ((mult-1)*multFraction)+1; //multFraction=0.5, for example, returns halfway between 1.0 and mult
	}
	
	public static float calcETotalOverOut(float v, float target, float weightMult) {
		float eTotalOverOut=v-target; //Aka out-target
		final float ret=toWeightedError(eTotalOverOut, v, target, weightMult);
		return ret;
//		final float mult=toErrorMult(v, target);
//		final float ret2=(float)(eTotalOverOut*mult);
//		
//		if(v>target) {
//			eTotalOverOut*=positiveErrorMult;
//			if(v>lowThresh && target<lowThresh) {
//				eTotalOverOut*=falsePositiveErrorMult;
//			}
//		}else{
//			eTotalOverOut*=negativeErrorMult;
//			if(v<highThresh && target>=highThresh) {
//				eTotalOverOut*=falseNegativeErrorMult;
//			}
//		}
//		assert(Tools.absdif(ret, ret2)<0.00001) : ret+", "+ret2;
//		assert(Tools.absdif(eTotalOverOut, ret2)<0.00001) : eTotalOverOut+", "+ret2;
//		
//		return (float)eTotalOverOut;
	}
	
	public static void setLowWeightAnnealCutoff(float c) {
		assert(c>=0 && c<=1);
		lowWeightAnnealCutoff=c;
		lowWeightAnnealMult=Tools.max(1f, 1f/Tools.max(lowWeightAnnealCutoff, 0.0001f));
	}
	
	/*--------------------------------------------------------------*/
	
	public final String typeString() {
		return function.name();
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	private double eTotalOverOut;
	private float outOverNet;
	private double eOverNet;
	private double sum=0;
	
	private float bias;
//	private double bias2;
	private double biasDelta;
	public double error;
	private final int id;
	final int lpos;//position within layer
	final int layer;
	int nextWeight=0; //For use when loading networks
	
//	int type;//TODO: Change to Function.
	Function function;
	
	/*--------------------------------------------------------------*/
	
	//Edges are not used in fast mode.
	public final ArrayList<Edge> inputs;
	public final ArrayList<Edge> outputs;
	
	/*--------------------------------------------------------------*/
	
	public final float[] weights;
//	private final double[] weights2;
	private final double[] deltas;
	private final float[] values; 
	private final double[] eOverNetArray;
	
	/*--------------------------------------------------------------*/
	/*----------------           Statics            ----------------*/
	/*--------------------------------------------------------------*/
	
	public static int MAX_TYPE=Function.TANH;
	public static int defaultActivationType=Function.SIG;
	public static int finalLayerType=Function.RSLOG;
	public static float randomTypeRate=0.0f;

	static float biasAlphaMult=1f;
	static float biasAnnealMult=0.5f;
	static boolean annealBias=true;
	
	private static float lowWeightAnnealCutoff=0.2f;
	private static float lowWeightAnnealMult=1f/lowWeightAnnealCutoff;
	
	static float cutoffForTraining=0.5f;
	static boolean setCutoffForTraining=false;

	static float positiveErrorMult=1.65f;
	static float falsePositiveErrorMult=12.5f;

	static float negativeErrorMult=0.825f;
	static float falseNegativeErrorMult=2.7f;
	
	static float fnErrorIncr=0.01f;
	static float fpErrorIncr=0.00f;
	static float spread=0.025f;
	
}
