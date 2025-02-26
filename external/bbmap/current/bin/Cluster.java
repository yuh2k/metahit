package bin;

import java.util.ArrayList;
import java.util.Iterator;

import shared.Tools;
import sketch.Sketch;
import sketch.SketchMakerMini;
import stream.Read;
import structures.IntHashMap;
import structures.IntHashSet;

public class Cluster extends Bin {

	public Cluster(int id_) {id=id_;}
	
	public Cluster(Contig c) {
		id=c.id();
		add(c);
	}
	
	@Override
	public boolean isCluster() {return true;}
	
	@Override
	public Cluster toCluster() {return this;}
	
//	@Override
//	public int clusterID() {return id();}
	
	@Override
	public Cluster cluster() {return this;}
	
	@Override
	public final int id() {return id;}
	
	@Override
	@Deprecated
	public void setID(int id_) {
		throw new RuntimeException("Not supported.");
	}
	
	public Cluster add(Bin b) {
		if(b.getClass()==Cluster.class) {return add((Cluster)b);}
		else {return add((Contig)b);}
	}
	
	public Cluster add(Cluster clust) {
		assert(clust!=this) : id();
		for(Contig c : clust.contigs) {add(c);}
		clust.contigs=null;
//		assert(isValid());
		return this;
	}
	
	boolean sameCluster(Bin b) {
		if(contigSet.contains(b.id())) {return true;}
		return b.isCluster() && ((Cluster)b).contigSet.contains(id());
	}
	
	public Cluster add(Contig c) {
		timestamp=globalTime;
//		verbose=(id()==52 || id()==3 || id()==3 || c.id()==52 || c.id()==3 || c.id()==3);
//		if(verbose) {
//			System.err.println("Adding "+c.id()+" to "+this.id()+(c.cluster==null ? "" : " from "+c.cluster.id));
//			assert(c.cluster==null || c.cluster.id()!=3);
//		}
//		assert(c.isValid());
		assert(c.cluster!=this);
		assert(!contigSet.contains(c.id()));
		c.cluster=this;
		if(contigs.size()==0) {
			topHit=c.topHit;
			secondHit=c.secondHit;
			taxid=c.taxid;
			genusTaxid=c.genusTaxid;
			labelTaxid=c.labelTaxid;
			sketchedSize=c.sketchedSize;
			entropy=c.entropy;
		}else {
//			assert(isValid());
			entropy=(entropy*size()+c.entropy*c.size())/(size+c.size());
		}
		contigs.add(c);
		contigSet.add(c.id());
		kmers+=c.kmers;
		invKmers=1f/kmers;
		if(counts==null) {counts=c.counts.clone();}
		else {Tools.add(counts, c.counts);}
		
		if(pairMap==null && c.pairMap!=null) {
			pairMap=new IntHashMap(5);
			pairMap.incrementAll(c.pairMap);//Remember, the targets are contig IDs.
		}
		
		if(numDepths()==0) {//Empty cluster
			assert(size==0 || grading);
			for(int i=0, max=c.numDepths(); i<max; i++) {appendDepth(c.depth(i));}
		}else {//Nonempty cluster
			assert(c.numDepths()==numDepths()) : c.numDepths()+", "+numDepths();
			float mult=1f/(size+c.size());
			float cmult=mult*c.size();
			float bmult=mult*size;
			for(int i=0, max=c.numDepths(); i<max; i++) {
				float cdepth=c.depth(i)*cmult;
				float bdepth=depth(i)*bmult;
				setDepth(cdepth+bdepth, i);
			}
		}
		if(numDepths()>1) {
			fillNormDepth();
		}
		avgDepthValid=false;
		size+=c.size();
		gcSum+=c.gcSum;
		assert(gcSum>0 || grading);
//		assert(isValid());
//		assert(c.isValid());
		return this;
	}
	
	@Override
	public long size() {return size;}

	@Override
	public Sketch toSketch(SketchMakerMini smm, Read r) {
		String name=Long.toString(id());
		if(r==null) {r=new Read(null, null, name, id());}
		r.id=name;
		r.numericID=id();
		for(Contig c : contigs) {
			r.bases=c.bases;
			smm.processReadNucleotide(r);
		}
		return smm.toSketch(0);
	}
	
	@Override
	public boolean isValid() {
		if(id()<0) {
			assert(false) : id()+", "+contigSet;
			return false;
		}
		if(numDepths()<1) {
			assert(false) : id()+", "+contigSet;
			return false;
		}
		if(counts==null) {
			assert(false) : id()+", "+contigSet;
			return false;
		}
		if(contigs.isEmpty()) {
			assert(false) : id()+", "+contigSet;
			return false;
		}
		for(Contig c : contigs) {
			if(c.cluster!=this) {
				assert(c.cluster!=null);
				assert(false) : "Cluster "+id()+" contains a contig "+c.id()+" that points to "+c.cluster().id()+"\n"
						+ " set="+contigSet+", c.set="+c.cluster.contigSet+"\n"+contigs+"\n"+c.cluster.contigs+"\n";
				return false;
			}
			if(!contigSet.contains(c.id())) {
				assert(false) : id()+", "+contigSet;
				return false;
			}
		}
		return true;
	}
	
	@Override
	public int numContigs() {return contigs==null ? 0 : contigs.size();}
	
	@Override
	public Iterator<Contig> iterator() {
		return contigs.iterator();
	}
	
	final int id;
	public long size;
	public ArrayList<Contig> contigs=new ArrayList<Contig>(8);
	public IntHashSet contigSet=new IntHashSet(5);
	private int timestamp=-1;
	
}
