package bin;

import java.util.ArrayList;
import java.util.Iterator;

import dna.AminoAcid;
import shared.Shared;
import shared.Tools;
import sketch.Sketch;
import sketch.SketchMakerMini;
import stream.Read;
import structures.ByteBuilder;

public class Contig extends Bin {

	public Contig(String name_, byte[] bases_, int id_) {
		name=name_;
		shortName=ContigRenamer.toShortName(name);
		bases=bases_;
		id=id_;
	}

//	public Contig(String name_, byte[] bases_, int id_) {
//		name=name_;
//		shortName=ContigRenamer.toShortName(name);
//		bases=bases_;
//		setID(id_);
//	}
	
	@Override
	public boolean isCluster() {return false;}
	
	@Override
	public Cluster toCluster() {
		assert(cluster==null);
		return new Cluster(this);
	}
	
	boolean sameCluster(Bin b) {
		return cluster!=null && cluster.contigSet.contains(b.id());
	}
	
//	@Override
//	public int clusterID() {return clusterID;}
	
	@Override
	public Cluster cluster() {return cluster;}
	
	@Override
	public void setID(int id_) {
//		assert(id==-1 && id_>-1);
		id=id_;
	}
	
	@Override
	public int id() {
//		assert(id>=0);
		return id;
	}
	
	@Override
	public boolean isValid() {
		if(numDepths()<1) {return false;}
		if(counts==null) {return false;}
		if(cluster!=null && !cluster.contigSet.contains(id())) {
			assert(false) : id()+", "+cluster.id+", "+cluster.contigSet;
			return false;
		}
		return true;
	}
	
	public void loadCounts() {
		assert(kmers==0);
		counts=new int[canonicalKmers];
		kmers=countKmers(bases, counts);
		invKmers=1f/Tools.max(1, kmers);
		for(byte b : bases) {
			int x=AminoAcid.baseToNumber[b];
			gcSum+=(x==1 || x==2) ? 1 : 0;
		}
		//This assertion can fail; I saw an all A/T contig
//		assert(gcSum>0) : gcSum+", "+kmers+", "+new String(bases);
	}
	
	/** In fasta format */
	public void appendTo(ByteBuilder bb, int cluster) {
		bb.append('>').append(name);
		if(cluster>=0) {bb.tab().append("cluster_").append(cluster);}
		bb.nl();
		final int wrap=Shared.FASTA_WRAP;
		for(int i=0; i<bases.length; i+=wrap) {
			//Now with modified append I can just append(bases, wrap)
			bb.append(bases, i, wrap).nl();
		}
	}
	
	@Override
	public long size() {return bases.length;}

	@Override
	public Sketch toSketch(SketchMakerMini smm, Read r) {
		String name=Long.toString(id());
		if(r==null) {r=new Read(null, null, name, id());}
		r.id=name;
		r.numericID=id();
		r.bases=bases;
		smm.processReadNucleotide(r);
		return smm.toSketch(0);
	}
	
	public final ByteBuilder toCov(ByteBuilder bb) {
		if(bb==null) {bb=new ByteBuilder();}
		bb.append(shortName);
		bb.tab().append(id());
		bb.tab().append(size());
		for(int i=0, max=numDepths(); i<max; i++) {
			bb.tab().append(depth(i), 2);
		}
		ArrayList<KeyValue> list=KeyValue.toList(pairMap);
		if(list!=null) {
			for(KeyValue ip : list) {
				bb.tab().append(ip.key).tab().append(ip.value);
			}
		}
		return bb;
	}
	
	@Override
	public int numContigs() {return 1;}
	
	@Override
	public Iterator<Contig> iterator() {
		return new ContigIterator();
	}
	
	private class ContigIterator implements Iterator<Contig> {

		@Override
		public boolean hasNext() {
			return hasMore;
		}

		@Override
		public Contig next() {
			if(!hasMore) {return null;}
			hasMore=false;
			return Contig.this;
		}
		
		boolean hasMore=true;
		
	}
	
	private int id=-1;
	public Cluster cluster=null;
	public final String name;
	public final String shortName;
	public final byte[] bases;
	
}
