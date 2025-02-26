package bin;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.concurrent.locks.ReadWriteLock;

import shared.LineParserS1;
import shared.LineParserS4;
import shared.Shared;
import shared.Tools;
import structures.IntLongHashMap;
import template.Accumulator;
import template.ThreadWaiter;
import tracker.EntropyTracker;

public class SpectraCounter implements Accumulator<SpectraCounter.LoadThread> {
	
	public SpectraCounter(PrintStream outstream_, boolean parseTax_, boolean parseDepth_) {
		outstream=outstream_;
		parseTax=parseTax_;
		parseDepth=parseDepth_;
	}
	
	/** Spawn process threads */
	public void makeSpectra(ArrayList<Contig> contigs){
		
		//Do anything necessary prior to processing
		sizeMap=(parseTax ? new IntLongHashMap(1021) : null);
		
		//Determine how many threads may be used
		int threads=Tools.mid(1, contigs.size()/4, Shared.threads());
		//Fill a list with LoadThreads
		ArrayList<LoadThread> alpt=new ArrayList<LoadThread>(threads);
		for(int i=0; i<threads; i++){
			LoadThread lt=new LoadThread(contigs, i, threads);
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
			errorState|=(t.success);
		}
	}

	@Override
	public ReadWriteLock rwlock() {return null;}

	@Override
	public synchronized boolean success() {return errorState;}
	
	class LoadThread extends Thread {
		
		LoadThread(ArrayList<Contig> contigs_, int tid_, int threads_) {
			contigs=contigs_;
			tid=tid_;
			threads=threads_;
			sizeMapT=(parseTax ? new IntLongHashMap(1021) : null);
		}
		
		@Override
		public void run() {
			synchronized(this) {
				runInner();
				if(parseTax) {
					synchronized(sizeMap) {
						sizeMap.incrementAll(sizeMapT);
					}
				}
			}
		}
		
		private void runInner() {
			LineParserS1 lps=new LineParserS1('_');
			LineParserS4 lpt=new LineParserS4(",,=,");
			for(int i=tid; i<contigs.size(); i+=threads) {
				Contig c=contigs.get(i);
				synchronized(c) {
					contigsProcessedT++;
					basesProcessedT+=c.size();
					c.loadCounts();
					c.fillNormDepth();
					c.entropy=(calcEntropy ? et.averageEntropy(c.bases, false) : 1);
					if(parseTax) {
						int tid=DataLoader.parseTaxID(c.name);
						if(tid>0) {
							c.labelTaxid=tid;
							if(BinObject.validation) {
								sizeMapT.increment(tid, c.size());
							}else {
								if(!BinObject.validation) {c.taxid=tid;}
							}
						}
					}
					if(parseDepth) {
						boolean b=DataLoader.parseAndSetDepth(c, lps, lpt);
						assert(b) : "Could not parse depth from "+c.name;
					}
					
					assert(c.counts!=null && c.kmers>0);
				}
			}
			success=true;
		}
		
		final IntLongHashMap sizeMapT;
		final int tid;
		final int threads;
		final ArrayList<Contig> contigs;
		final EntropyTracker et=new EntropyTracker(5, 50, false);
		boolean success=false;
		int contigsProcessedT=0;
		long basesProcessedT=0;
	}
	
	public PrintStream outstream=System.err;

	public final boolean parseTax;
	public final boolean parseDepth;
	public IntLongHashMap sizeMap;
	public boolean errorState=false;
	public static boolean calcEntropy=false;
	
}
