package barcode;

import java.util.HashMap;

import align2.BandedAligner;
import dna.AminoAcid;
import shared.Tools;
import structures.ByteBuilder;

public class Barcode implements Comparable<Barcode> {
	
	/*--------------------------------------------------------------*/
	/*----------------         Constructors         ----------------*/
	/*--------------------------------------------------------------*/

	public Barcode(String s){this(s, 0);}
	public Barcode(String s, long c){
		name=s;
		count=c;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Methods           ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Gets one of two barcodes for dual indexes */
	public Barcode getBarcodeForPairnum(int pairnum, char delimiter) {
		return new Barcode(getStringForPairnum(pairnum, delimiter), count);
	}
	public String getStringForPairnum(int pairnum, char delimiter) {
		assert(pairnum==0 || pairnum==1);
		final int pos, pos2;
		if(delimiter>0){//If there is a delimiter
			pos=name.indexOf(delimiter);
			pos2=pos+1;
			assert(pos>=0) : pos+", "+Character.toString(delimiter)+", "+name;
		}else{
			boolean even=((name.length()&1)==0);
			pos=(name.length()/2)+1;
			pos2=(even ? pos : pos+1);
		}
		String s=(pairnum==0 ? name.substring(0, pos) : name.substring(pos2, name.length()));
		return s;
	}

	public Barcode left(char delimiter){return getBarcodeForPairnum(0, delimiter);}
	public Barcode right(char delimiter){return getBarcodeForPairnum(1, delimiter);}
	public String leftString(char delimiter){return getStringForPairnum(0, delimiter);}
	public String rightString(char delimiter){return getStringForPairnum(1, delimiter);}
	
	/*--------------------------------------------------------------*/
	
	public int countUndefined() {
		int sum=0;
		for(int i=0; i<length(); i++){
			char c=charAt(i);
			sum+=(Tools.isLetter(c) && !AminoAcid.isFullyDefined(c)) ? 1 : 0;
		}
		return sum;
	}
	
	/** Returns -1 if false, and the homopolymer character 2-bit encoding if true */
	public byte checkHomopolymer() {
		int match=0;
		int nonletter=0;
		final byte mer=AminoAcid.baseToNumber[charAt(0)];
		if(mer<0) {return mer;}
		for(int i=1; i<length(); i++){
			char c=charAt(i);
			byte n=AminoAcid.baseToNumber[c];
			if(n==mer) {
				match++;
			}else if(Tools.isLetter(c)) {
				return -1;
			}else {
				nonletter++;
			}
		}
		return (match>0 && nonletter<=1 ? mer : 0);
	}
	
	public int hdist(final Barcode b){return hdist(b.name);}
	public int hdist(final String b){
		assert(length()==b.length());
		final int min=Tools.min(length(), b.length());
		int subs=0;
		for(int i=0; i<min; i++){
			final char ca=charAt(i), cb=b.charAt(i);
			subs+=(ca==cb ? 0 : 1);
		}
		return subs;//+Tools.absdif(length(), b.length());
	}

	public int edist(final Barcode b, BandedAligner bandy){return edist(b.name, bandy);}
	public int edist(final String b, BandedAligner bandy){
		int dist=hdist(b);
		if(dist>1){dist=bandy.alignForward(getBytes(), b.getBytes(), 0, 0, length(), true);}
		return dist;
	}
	//Bandy is, e.g., new BandedAlignerConcrete(21);
	
	/*--------------------------------------------------------------*/

	public long count(){return count;}
	public void setCount(long x){count=x;}
	public int length() {return name.length();}
	public char charAt(int i) {return name.charAt(i);}
	public byte[] getBytes() {return name.getBytes();}
	
	/*--------------------------------------------------------------*/
	/*----------------           Mutators           ----------------*/
	/*--------------------------------------------------------------*/
	
	public void increment() {count++;}
	public void increment(Barcode b) {count+=b.count;}
	public void increment(long x) {count+=x;}
	
	/*--------------------------------------------------------------*/
	/*----------------           Overrides          ----------------*/
	/*--------------------------------------------------------------*/

	@Override
	public int hashCode(){
		return name.hashCode();
	}
	
	@Override
	public int compareTo(Barcode b) {
		if(count!=b.count){return count>b.count ? -1 : 1;}
		return name.compareTo(b.name);
	}
	
	@Override
	public boolean equals(Object b) {
//		return(b!=null && getClass()==b.getClass() && equals((Barcode)b));//Proper way
		return(equals((Barcode)b));//Faster way
	}
	public boolean equals(Barcode b) {return(b!=null && name.equals(b.name));}//Ignores count
	
	@Override
	public String toString(){
		return name+"\t"+count;
	}
	
	public ByteBuilder appendTo(ByteBuilder bb){
		return bb.append(name).tab().append(count);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	public final String name;
	private long count=0;
	
}
