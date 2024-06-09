package shared;

import java.io.File;
import java.util.ArrayList;

import fileIO.ByteFile;

/** Similar speed, but less powerful.
 * Main advantage is having a bounded memory footprint for very long lines.
 * */
public final class LineParser2 implements LineParser {
	
	/*--------------------------------------------------------------*/
	/*----------------             Main             ----------------*/
	/*--------------------------------------------------------------*/
	
	//For testing
	//Syntax: LineParser fname/literal delimiter 
	public static void main(String[] args) {
		assert(args.length==2);
		String fname=args[0];
		String dstring=args[1];
		assert(dstring.length()==1);
		
		final ArrayList<byte[]> lines;
		if(new File(fname).exists()){
			lines=ByteFile.toLines(fname);
		}else{
			lines=new ArrayList<byte[]>(1);
			lines.add(fname.getBytes());
		}
		
		LineParser2 lp=new LineParser2(dstring.charAt(0));
		for(byte[] line : lines) {
			lp.set(line);
			System.out.println(lp);
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Constructors         ----------------*/
	/*--------------------------------------------------------------*/

	public LineParser2(byte delimiter_) {delimiter=delimiter_;}

	public LineParser2(int delimiter_) {
		assert(delimiter_>=0 && delimiter_<=127);
		delimiter=(byte)delimiter_;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	@Override
	public LineParser2 set(byte[] line_) {
		reset();
		line=line_;
		return this;
	}
	
	@Override
	public LineParser2 set(byte[] line_, int maxTerm) {
		return set(line_);
	}
	
	@Override
	public LineParser2 clear() {
		line=null;
		a=b=currentTerm=-1;
		return this;
	}
	
	@Override
	public LineParser2 reset() {
		a=b=currentTerm=-1;
		return this;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Parse Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	public int parseInt() {
		advance();
		return Parse.parseInt(line, a, b);
	}
	
	public long parseLong() {
		advance();
		return Parse.parseLong(line, a, b);
	}
	
	public float parseFloat() {
		advance();
		return Parse.parseFloat(line, a, b);
	}
	
	public double parseDouble() {
		advance();
		return Parse.parseDouble(line, a, b);
	}
	
	public byte parseByte(int offset) {
		advance();
		int index=a+offset;
		assert(index<b);
		return line[index];
	}
	
	public String parseString() {
		int len=advance();
		return new String(line, a, len);
	}
	
	/*--------------------------------------------------------------*/
	
	@Override
	public int parseInt(int term) {
		advanceTo(term);
		return Parse.parseInt(line, a, b);
	}

	@Override
	public long parseLong(int term) {
		advanceTo(term);
		return Parse.parseLong(line, a, b);
	}

	@Override
	public float parseFloat(int term) {
		advanceTo(term);
		return Parse.parseFloat(line, a, b);
	}

	@Override
	public double parseDouble(int term) {
		advanceTo(term);
		return Parse.parseDouble(line, a, b);
	}

	@Override
	public byte parseByte(int term, int offset) {
		advanceTo(term);
		int index=a+offset;
		assert(index<b);
		return line[index];
	}

	@Override
	public String parseString(int term) {
		int len=advanceTo(term);
		return new String(line, a, len);
	}
	
	@Override
	public boolean startsWith(String s) {
		return Tools.startsWith(line, s);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Advance Methods       ----------------*/
	/*--------------------------------------------------------------*/
	
	public final int advance() {
		currentTerm++;
		b++;
		a=b;
		while(b<line.length && line[b]!=delimiter){b++;}
		return b-a;
	}
	
	public void advanceBy(int terms) {
		for(; terms>0; terms--) {
			advance();
		}
	}
	
	//Advances to term before toTerm
	public void advanceToBefore(int toTerm) {
		assert(toTerm>=currentTerm) : "Can't advance backwards: "+currentTerm+">"+toTerm;
		for(toTerm--; currentTerm<toTerm;) {
			advance();
		}
	}
	
	//Advances to actual term
	private int advanceTo(int toTerm) {
		assert(toTerm>=currentTerm) : "Can't advance backwards: "+currentTerm+">"+toTerm;
		for(toTerm--; currentTerm<=toTerm;) {
			advance();
		}
		return b-a;
	}
	
	public void incrementA(int incr) {
		a+=incr;
	}
	
	public void incrementB(int incr) {
		b+=incr;
	}
	
	public void setBounds(int a_, int b_) {
		a=a_;
		b=b_;
	}
	
	public boolean hasMore() {
		return b<line.length;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Methods           ----------------*/
	/*--------------------------------------------------------------*/
	
	@Override
	public String toString() {
		return /*toList().toString()+"\n"*/"a="+a+", b="+b+", line.length="+line.length;
	}
	
	public ArrayList<String> toList(){
		ArrayList<String> list=new ArrayList<String>();
		do{
			list.add(parseString());
		}while(b<line.length);
		return list;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	private int a=-1;
	private int b=-1;
	private int currentTerm=-1;
	private byte[] line;
	
	public final byte delimiter;
	
}
