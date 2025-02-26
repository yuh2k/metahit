package bin;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.concurrent.locks.ReadWriteLock;

import fileIO.FileFormat;
import shared.Shared;
import shared.Tools;
import stream.SamLine;
import stream.SamLineStreamer;
import structures.IntHashMap;
import structures.ListNum;
import template.Accumulator;
import template.ThreadWaiter;

public class SamLoader implements Accumulator<SamLoader.LoadThread> {
	
	public SamLoader(PrintStream outstream_) {
		outstream=outstream_;
	}
	
	@Deprecated
	public void load(ArrayList<String> fnames, HashMap<String, Contig> contigMap, IntHashMap[] graph) {
		//Contig list should already be sorted and numbered.
		ArrayList<Contig> list=new ArrayList<Contig>(contigMap.values());
		Collections.sort(list);
		for(int i=0; i<list.size(); i++) {list.get(i).setID(i);}
		load(fnames, contigMap, list, graph);
	}
	
	/** Spawn process threads */
	public void load(ArrayList<String> fnames, HashMap<String, Contig> contigMap, 
			ArrayList<Contig> contigs, IntHashMap[] graph){
		
		//Do anything necessary prior to processing
		
		//Determine how many threads may be used
		final int threads=fnames.size();
		
		//Fill a list with LoadThreads
		ArrayList<LoadThread> alpt=new ArrayList<LoadThread>(threads);
		for(int i=0; i<threads; i++){
			String fname=fnames.get(i);
			LoadThread lt=new LoadThread(fname, i, contigMap, contigs, graph);
			alpt.add(lt);
		}
		
		//Start the threads and wait for them to finish
		boolean success=ThreadWaiter.startAndWait(alpt, this);
		errorState&=!success;
		
		//Do anything necessary after processing
		
	}
	
	@Override
	public synchronized void accumulate(LoadThread t) {
		synchronized(t) {
			readsIn+=t.readsInT;
			basesIn+=t.basesInT;
			bytesIn+=t.bytesInT;
			errorState|=(t.success);
		}
	}

	@Override
	public ReadWriteLock rwlock() {return null;}

	@Override
	public synchronized boolean success() {
		return errorState;
	}
	
	class LoadThread extends Thread {
		
		LoadThread(final String fname_, final int sample_, HashMap<String, 
				Contig> contigMap_, ArrayList<Contig> contigs_, IntHashMap[] graph_) {
			fname=fname_;
			sample=sample_;
			contigMap=contigMap_;
			contigs=contigs_;
			graph=graph_;
		}
		
		@Override
		public void run() {
			synchronized(this) {runInner();}
		}
		
		private void runInner() {
			long[] depthArray=new long[contigs.size()];
			SamLineStreamer ss=null;
			outstream.println("Loading "+fname);
			final int streamerThreads=Tools.min(3, Shared.threads());
			FileFormat ff=FileFormat.testInput(fname, FileFormat.SAM, null, true, false);
			ss=new SamLineStreamer(ff, streamerThreads, false, -1);
			ss.start();
			processSam_Thread(ss, depthArray);
			
			postprocess(depthArray);
			depthArray=null;
			success=true;
		}
		
		private void postprocess(long[] depthArray) {
			for(int cnum=0; cnum<depthArray.length; cnum++) {
				Contig c=contigs.get(cnum);
				float depth=depthArray[cnum]*1f/Tools.max(1, c.size());
				synchronized(c) {c.setDepth(depth, sample);}
			}
		}
		
		void processSam_Thread(SamLineStreamer ss, long[] depthArray) {
			ListNum<SamLine> ln=ss.nextLines();
			ArrayList<SamLine> reads=(ln==null ? null : ln.list);

			while(ln!=null && reads!=null && reads.size()>0){

				for(int idx=0; idx<reads.size(); idx++){
					SamLine sl=reads.get(idx);
					if(sl.mapped()) {
						addSamLine(sl, depthArray);
						readsInT++;
						basesInT+=(sl.seq==null ? 0 : sl.length());
						bytesInT+=(sl.countBytes());
					}
				}
				ln=ss.nextLines();
				reads=(ln==null ? null : ln.list);
			}
		}
		
		private boolean addSamLine(SamLine sl, long[] depthArray) {
			final String rname=ContigRenamer.toShortName(sl.rname());
			final Contig c1=contigMap.get(rname);
			if(c1==null) {return false;}//Contig not found; possibly too short
			assert(c1!=null) : "Can't find contig for rname "+rname;
			final int cid=c1.id();
			depthArray[cid]+=sl.mappedNonClippedBases();
			
			if(graph==null || sl.ambiguous() || !sl.hasMate() || 
					!sl.nextMapped() || sl.pairedOnSameChrom() || sl.mapq<minMapq) {return true;}
			if(minMateq>0) {
				int mateq=sl.mateq();
				if(mateq>=0 && mateq<minMateq) {return true;}
			}
			final String rnext=ContigRenamer.toShortName(sl.rnext());
			assert(rnext!=null && !"*".equals(rnext) && !"=".equals(rnext));
			
			final Contig c2=contigMap.get(rnext);
			if(c2==null) {return false;}//Contig not found
			assert(c2!=null) : "Can't find contig for rnext "+rnext;
			
			IntHashMap destMap=graph[cid];
			if(destMap==null) {
				synchronized(graph) {
					if(graph[cid]==null) {graph[cid]=new IntHashMap(5);}
					destMap=graph[cid];
				}
			}
			synchronized(destMap) {
				destMap.increment(c2.id());
			}
			return true;
		}
		
		final String fname;
		final int sample;
		final HashMap<String, Contig> contigMap;
		final ArrayList<Contig> contigs;
		final IntHashMap[] graph;
		long readsInT=0;
		long basesInT=0;
		long bytesInT=0;
		boolean success=false;
	}
	
	public PrintStream outstream=System.err;
	public long readsIn=0;
	public long basesIn=0;
	public long bytesIn=0;
	public int minMapq=4;
	public int minMateq=4;
	
	public boolean errorState=false;
	
}
