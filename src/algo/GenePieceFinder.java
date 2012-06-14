package algo;

import java.util.ArrayList;
import java.util.List;

import model.CGene;
import model.Helix;
import utils.Utils;

public class GenePieceFinder implements GeneFinder {
	
	public static final int HEARTBEAT = 500000;
	public static final int BLOCK_SIZE = 1000;

	@Override
	public List<CGene> findGenes(String gbk) {
		List<CGene> all = new ArrayList<CGene>();
		for (int i = 0; i < gbk.length(); i+=BLOCK_SIZE/2) {
			if (i % HEARTBEAT == 0) {
				Utils.log(String.format("heartbeat %d/%d", i, gbk.length()));
			}
			int[] codes = Utils.toInt(gbk.substring(i,Math.min(i+BLOCK_SIZE,gbk.length())));
			List<Helix> helixes = HelixFinder.findHelixesUnefficient(codes, 4, 4, 1000000);
			Helix.addShift(helixes, i);
			DynamicSolver solver = new DynamicSolver(helixes);
			List<CGene> res = solver.solve();
			all.addAll(res);
		}
		return all;
	}
}
