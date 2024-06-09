package ml;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;

import fileIO.ByteFile;
import shared.LineParser2;
import shared.Parse;
import shared.Tools;
import structures.IntList;

public class CellNetParser {
	
	static CellNet parse(String fname, boolean nullOnFailure) {
		if(nullOnFailure && !new File(fname).exists()) {return null;}
		CellNetParser cnp=new CellNetParser(fname);
		return cnp.net;
	}

	public static CellNet load(String fname) {return load(fname, true);}
	public static CellNet load(String fname, boolean nullOnFailure) {
		return parse(fname, nullOnFailure);
	}
	
	private CellNetParser(String fname_){
		this(ByteFile.toLines(fname_));
		fname=fname_;
	}
	
	private CellNetParser(ArrayList<byte[]> lines_){
		lines=lines_;
		
		parseHeader();
		
		net=new CellNet(dims, seed, commands);
		net.epochsTrained=epochs;
		net.samplesTrained=samples;
		net.annealSeed=annealSeed;
		net.cutoff=cutoff;
		assert(layers==net.layers);
		posFirstEdge=pos;
		parseEdges();
	}
	
	public void parseHeader() {
		while(pos<lines.size()) {
			byte[] line=lines.get(pos);
			
			if(line.length<1){
				//ignore
			}else if(Tools.startsWith(line, "#")){//header
				if(Tools.startsWith(line, "##ctf") || Tools.startsWith(line, "#ctf")){
					cutoff=parseFloat(line);
				}else if(Tools.startsWith(line, "##")){
					//Comment; ignore
				}else if(Tools.startsWith(line, "#version")){
					version=parseInt(line);
				}else if(Tools.startsWith(line, "#layers")){
					layers=parseInt(line);
				}else if(Tools.startsWith(line, "#seed")){
					seed=parseLong(line);
				}else if(Tools.startsWith(line, "#annealseed")){
					annealSeed=parseLong(line);
				}else if(Tools.startsWith(line, "#epochs")){
					epochs=parseLong(line);
				}else if(Tools.startsWith(line, "#samples")){
					samples=parseLong(line);
				}else if(Tools.startsWith(line, "#concise")){
					concise=true;
				}else if(Tools.startsWith(line, "#dims")){
					dims=parseIntArray(line, delimiter, true);
					assert(layers==dims.length) : layers+", "+Arrays.toString(dims);
				}else if(Tools.startsWith(line, "#edges")){
					//ignore
				}else if(Tools.startsWith(line, "#CL")){
					commands.add(new String(line));
				}else if(Tools.startsWith(line, "#")){
					assert(false) : "\nUnexpected header line: '"+new String(line)+"'"
							+ "\nComments should start with ##\n";
				}
			}else{
				break; //A cell or edge
			}
			pos++;
		}
//		assert(false) : pos+", "+new String(lines.get(pos));
	}
	
	private void parseEdges() {
		pos=posFirstEdge;
		long numEdges=0;
		
		int numCells=(int) shared.Vector.sum(dims);
		LineParser2 lp=new LineParser2(delimiter);
		
		while(pos<lines.size()) {
			byte[] line=lines.get(pos);
			
			if(line.length<1){
				//ignore
			}else if(Tools.startsWith(line, "##")){
				//ignore
			}else if(Tools.startsWith(line, 'C')){
				//parse cell
//				int a=1, b=1;
//				while(b<line.length && line[b]!=delimiter){b++;}
//				int cid=Parse.parseInt(line, a, b);
//				b++;
//				a=b;
//
//				while(b<line.length && line[b]!=delimiter){b++;}
//				String s=new String(line, a, b-a);
//				int type=Tools.find(s, Function.TYPES);
//				assert(type>=0) : type+", "+s+"\n'"+new String(line)+"'"+", "+a+", "+b;
//				Cell c=net.list.get(cid);
////				c.type=type;
//				c.function=Function.getFunction(type);
//				//assert(false) : cid;
				
				lp.set(line);
				lp.setBounds(0, 0);
				int cid=lp.parseInt();
//				assert(false) : cid;
				String s=lp.parseString();
				int type=Tools.find(s, Function.TYPES);
				assert(type>=0) : type+", "+s+"\n'"+new String(line)+"'";
				Cell c=net.list.get(cid);
				c.function=Function.getFunction(type);
				assert(c.function.type()==type);
				
				if(concise) {
					c.setBias(lp.parseFloat(), true);
					for(int i=0; i<c.weights.length; i++) {c.weights[i]=lp.parseFloat();}
					assert(!lp.hasMore()) : lp.toString();
				}
			}else{
				//parse edge
//				int a=0, b=0;
//
//				while(b<line.length && line[b]!=delimiter){b++;}
//				assert(b>a) : "Missing source : '"+new String(line)+"'";
//				int sid=Parse.parseInt(line, a, b);
//				b++;
//				a=b;
//
//				while(b<line.length && line[b]!=delimiter){b++;}
//				assert(b>a) : "Missing dest : '"+new String(line)+"'";
//				int cid=Parse.parseInt(line, a, b);
//				b++;
//				a=b;
//
//				while(b<line.length && line[b]!=delimiter){b++;}
//				assert(b>a) : "Missing weight : '"+new String(line)+"'";
//				float weight=Parse.parseFloat(line, a, b);
//				b++;
//				a=b;
				
				lp.set(line);
				final int sid=lp.parseInt();
				final int cid=lp.parseInt();
				final float weight=lp.parseFloat();

				Cell c=net.list.get(cid);
				if(sid==cid || sid<1) {
					c.setBias(weight, true);
//					System.err.println("Set bias to "+weight);
				}else {
					c.weights[c.nextWeight]=weight;
					c.nextWeight++;
					if(!CellNet.FAST){
						Cell s=net.list.get(sid);

						numEdges++;
						Edge e=new Edge(numEdges, s, c, weight);
						c.inputs.add(e);
						s.outputs.add(e);
						assert(c.nextWeight==c.inputs.size());
//						System.err.println("Added "+e);
					}
				}
			}
			pos++;
		}
		assert(CellNet.FAST || net.check());
	}
	
	boolean hasMore(){
		return pos<lines.size();
	}
	
	byte[] nextLine(){
		if(pos>=lines.size()){return null;}
		byte[] line=lines.get(pos);
		pos++;
		return line;
	}
	
	String fname;
	final ArrayList<byte[]> lines;
	private final CellNet net;
	long seed;
	long annealSeed=-1;
	long epochs=0;
	long samples=0;
	int version;
	int layers=-1;
	boolean concise=false;
	int[] dims;
	int pos=0;
	float cutoff=0.5f;
	final int posFirstEdge;
	ArrayList<String> commands=new ArrayList<String>();
	
	public static final byte delimiter=' ';
	
	/*--------------------------------------------------------------*/
	/*----------------           Parsing            ----------------*/
	/*--------------------------------------------------------------*/
	
	private static String parseString(byte[] line){
		int idx=Tools.indexOf(line, delimiter);
		String s=new String(line, idx+1, line.length-idx-1);
		return s;
	}
	
	private static float parseFloat(byte[] line){
		int idx=Tools.indexOf(line, delimiter);
		return Parse.parseFloat(line, idx+1, line.length);
	}
	
	private static int parseInt(byte[] line){
		int idx=Tools.indexOf(line, delimiter);
		return Parse.parseInt(line, idx+1, line.length);
	}
	
	private static long parseLong(byte[] line){
		int idx=Tools.indexOf(line, delimiter);
		return Parse.parseLong(line, idx+1, line.length);
	}
	
	public static int[] parseIntArray(final byte[] line, final byte delimiter, boolean parseTitle){
		int a=0, b=0;
		IntList list=new IntList(3);
		
		if(parseTitle) {
			while(b<line.length && line[b]!=delimiter){b++;}
			assert(b>a) : "Missing Title: "+new String(line);
			b++;
			a=b;
		}
		
		while(a<line.length) {
			while(b<line.length && line[b]!=delimiter){b++;}
			assert(b>a) : "Missing element "+list.size+": '"+new String(line)+"'";
			int x=Parse.parseInt(line, a, b);
//			assert(x>0) : new String(line);
			list.add(x);
			b++;
			a=b;
		}
		return list.toArray();
	}
	
	
}
