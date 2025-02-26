package bin;

import java.util.ArrayList;
import java.util.Collections;

import fileIO.ByteStreamWriter;
import shared.Tools;

public class ChartMaker {
	
	static void makeChartFromBins(String fname, ArrayList<Bin> list) {
		Collections.sort(list, new BinComparator());
		ByteStreamWriter bsw=new ByteStreamWriter(fname, true, false, false);
		bsw.start();
		bsw.print("#Bin\tSize\tClean\tDirty\n");
		double size=0, clean=0, dirty=0;
		int i=0;
		for(Bin b : list) {
			double contam=b.size()*Tools.max(0, b.contam);
			size+=b.size();
			clean+=(b.size()-contam);
			dirty+=contam;
			bsw.print(i).tab().print((long)size).tab().print((long)clean).tab().print((long)dirty).println();
			i++;
		}
		bsw.poisonAndWait();
	}
	
	static void makeChartFromBinStats(String fname, ArrayList<BinStats> list) {
		Collections.sort(list, new BinStatsComparator());
		ByteStreamWriter bsw=new ByteStreamWriter(fname, true, false, false);
		bsw.start();
		bsw.print("#Bin\tSize\tClean\tDirty\n");
		double size=0, clean=0, dirty=0;
		int i=0;
		for(BinStats b : list) {
			double contam=b.size*Tools.max(0, b.contam);
			size+=b.size;
			clean+=(b.size-contam);
			dirty+=contam;
			bsw.print(i).tab().print((long)size).tab().print((long)clean).tab().print((long)dirty).println();
			i++;
		}
		bsw.poisonAndWait();
	}
	
}
