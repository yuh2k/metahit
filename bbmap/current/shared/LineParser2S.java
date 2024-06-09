package shared;

import java.io.File;
import java.util.ArrayList;

import fileIO.ByteFile;
import fileIO.TextFile;

/** Similar speed, but less powerful.
 * Main advantage is having a bounded memory footprint for very long lines.
 * */
public final class LineParser2S implements LineParser {
	
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
		
		final String[] lines;
		if(new File(fname).exists()){
			lines=TextFile.toStringLines(fname);
		}else{
			lines=new String[] {fname};
		}
		
		LineParser2S lp=new LineParser2S(dstring.charAt(0));
		for(String line : lines) {
			lp.set(line);
			System.out.println(lp);
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Constructors         ----------------*/
	/*--------------------------------------------------------------*/

	public LineParser2S(char delimiter_) {delimiter=delimiter_;}

	public LineParser2S(int delimiter_) {
		assert(delimiter_>=0 && delimiter_<=Character.MAX_VALUE);
		delimiter=(char)delimiter_;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/

	@Override
	public LineParser2S set(byte[] line_) {
		assert(false) : "Use byte version.";
		line=new String(line_);
		return this;
	}

	@Override
	public LineParser2S set(byte[] line_, int maxTerm) {
		return set(line_);
	}
	
	public LineParser2S set(String line_) {
		reset();
		line=line_;
		return this;
	}
	
	public LineParser2S set(String line_, int maxTerm) {
		return set(line_);
	}
	
	public LineParser2S clear() {
		line=null;
		a=b=currentTerm=-1;
		return this;
	}
	
	public LineParser2S reset() {
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
		return (byte)line.charAt(index);
	}
	
	public String parseString() {
		int len=advance();
		return line.substring(a, b);
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
		return (byte)line.charAt(index);
	}

	@Override
	public String parseString(int term) {
		int len=advanceTo(term);
		return line.substring(a, b);
	}
	
	@Override
	public boolean startsWith(String s) {
		return line.startsWith(s);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Advance Methods       ----------------*/
	/*--------------------------------------------------------------*/
	
	public final int advance() {
		currentTerm++;
		b++;
		a=b;
		while(b<line.length() && line.charAt(b)!=delimiter){b++;}
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
	
	/*--------------------------------------------------------------*/
	/*----------------            Methods           ----------------*/
	/*--------------------------------------------------------------*/
	
	@Override
	public String toString() {
		return toList().toString();
	}
	
	public ArrayList<String> toList(){
		ArrayList<String> list=new ArrayList<String>();
		do{
			list.add(parseString());
		}while(b<line.length());
		return list;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	private int a=-1;
	private int b=-1;
	private int currentTerm=-1;
	private String line;
	
	public final char delimiter;
	
}
