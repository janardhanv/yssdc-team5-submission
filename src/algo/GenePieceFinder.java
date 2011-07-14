package algo;

import java.util.ArrayList;
import java.util.List;

import model.CGene;
import model.Helix;
import utils.Utils;

public class GenePieceFinder implements GeneFinder {

	@Override
	public List<CGene> findGenes(String gbk) {
		List<CGene> all = new ArrayList<CGene>();
		for (int i = 0; i < gbk.length(); i+=200) {
			if (i % 1000000 == 0) {
				Utils.log(String.format("heartbeat %d/%d", i, gbk.length()));
			}
			int[] codes = Utils.toInt(gbk.substring(i,Math.min(i+400,gbk.length())));
			List<Helix> helixes = HelixFinder.findHelixesUnefficient(codes, 4, 4, 1000000);
			Helix.addShift(helixes, i);
			DynamicSolver solver = new DynamicSolver(helixes);
			List<CGene> res = solver.solve();
			all.addAll(res);
		}
		return all;
	}
}
