package shared;

import java.io.File;
import java.util.ArrayList;

import fileIO.ByteFile;
import structures.IntList;

/**
 * Finds delimiters of a text line efficiently, to allow for parsing.
 * For example:<br>
 * Integer.parseInt("a b c 22 jan".split(" ")[3])<br>
 * could be redone as:<br>
 * LineParser lp=new LineParser(' ')<br>
 * lp.set("a b c 22 jan".toBytes()).parseInt(3)<br>
 * 
 * Uses memory proportional to 4*(# delimiters per line); for constant memory, use LineParser2.
 * 
 * @author BBushnell
 *
 */
public final class LineParser1 implements LineParser {
	
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
		
		LineParser lp=new LineParser1(dstring.charAt(0));
		for(byte[] line : lines) {
			lp.set(line);
			System.out.println(lp);
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Constructors         ----------------*/
	/*--------------------------------------------------------------*/

	public LineParser1(byte delimiter_) {delimiter=delimiter_;}

	public LineParser1(int delimiter_) {
		assert(delimiter_>=0 && delimiter_<=127);
		delimiter=(byte)delimiter_;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	@Override
	public LineParser set(byte[] line_) {
		clear();
		line=line_;
		for(int len=advance(); b<line.length; len=advance()) {
			bounds.add(b);
		}
		bounds.add(b);
		return this;
	}
	
	@Override
	public LineParser set(byte[] line_, int maxTerm) {
		clear();
		line=line_;
		for(int term=0; term<=maxTerm; term++) {
			int len=advance();
			bounds.add(b);
		}
		return this;
	}
	
	@Override
	public LineParser clear() {
		line=null;
		a=b=-1;
		bounds.clear();
		return this;
	}
	
	@Override
	public LineParser reset() {
		//Does nothing for this class
		return this;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Parse Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	public int terms() {return bounds.size();}
	
	@Override
	public int parseInt(int term) {
		setBounds(term);
		return Parse.parseInt(line, a, b);
	}
	
	@Override
	public long parseLong(int term) {
		setBounds(term);
		return Parse.parseLong(line, a, b);
	}
	
	@Override
	public float parseFloat(int term) {
		setBounds(term);
		return Parse.parseFloat(line, a, b);
	}
	
	@Override
	public double parseDouble(int term) {
		setBounds(term);
		return Parse.parseDouble(line, a, b);
	}
	
	@Override
	public byte parseByte(int term, int offset) {
		setBounds(term);
		final int index=a+offset;
		assert(index<b);
		return line[index];
	}
	
	@Override
	public String parseString(int term) {
		final int len=setBounds(term);
		return new String(line, a, len);
	}
	
	@Override
	public boolean startsWith(String s) {
		return Tools.startsWith(line, s);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Private Methods       ----------------*/
	/*--------------------------------------------------------------*/
	
	private int setBounds(int term){
		a=(term==0 ? 0 : bounds.get(term-1)+1);
		b=bounds.get(term);
		return b-a;
	}
	
	private int advance() {
		b++;
		a=b;
		while(b<line.length && line[b]!=delimiter){b++;}
		return b-a;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Methods           ----------------*/
	/*--------------------------------------------------------------*/
	
	@Override
	public String toString() {
		return toList().toString();
	}
	
	@Override
	public ArrayList<String> toList(){
		ArrayList<String> list=new ArrayList<String>(bounds.size);
		for(int i=0; i<bounds.size; i++){
			list.add(parseString(i));
		}
		return list;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	private final IntList bounds=new IntList();
	
	private int a=-1;
	private int b=-1;
	private byte[] line;
	
	public final byte delimiter;
	
}
