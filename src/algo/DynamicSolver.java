package algo;
import io.GbkReader;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

import model.Constraints;
import model.GlobalContext;
import model.Helix;
import utils.Utils;

public class DynamicSolver {
	public static class State {
		public Helix h;
		public int answer;
		public int cumulativeAnswer;
		public State(Helix h) {
			super();
			this.h = h;
			answer = -1;
		}
		public double getScore() {
			return answer*2.0/h.getSpan();
		}
		@Override
		public String toString() {
			return "Score "+getScore()+" span "+h.getSpan()+" top "+h;
		}
	}
	
	public State[] a;

	public DynamicSolver(List<Helix> h) {
		super();
		a = new State[h.size()];
		for (int i = 0; i < h.size(); ++i)
			a[i] = new State(h.get(i));
	}
	
	public List<State> solve() {
		Arrays.sort(a,new Comparator<State>() {
			@Override
			public int compare(State x, State y) {
				if (x.h.start == y.h.start && x.h.end == y.h.end)
					return 0;
				return x.h.end < y.h.end || x.h.end == y.h.end && x.h.start > y.h.start ? -1 : 1;
			}
		});
		State dummy = new State(new Helix(-1,-1,0));
		dummy.cumulativeAnswer = 0;
		List<State> res = new ArrayList<State>();
		for (int i = 0; i < a.length; i++) {
			ArrayList<State> inside = new ArrayList<State>();
			inside.add(dummy);
			for (int j = 0; j < i; j++) {
				if (!(a[i].h.start < a[j].h.getLeft() && a[j].h.getRight() < a[i].h.end)
						|| a[i].answer == 0)
					continue;
				int p = inside.size()-1;
				while (inside.get(p).h.getRight() >= a[j].h.getLeft())
					--p;
				a[j].cumulativeAnswer = inside.get(p).cumulativeAnswer + a[j].answer;
				a[j].cumulativeAnswer = Math.max(a[j].cumulativeAnswer, inside.get(inside.size()-1).cumulativeAnswer);
				inside.add(a[j]);
			}
			a[i].answer = inside.get(inside.size()-1).cumulativeAnswer;
			// there is something inside or we have a valid hairpin
			if (a[i].answer > 0 || a[i].h.isHairpin())
				a[i].answer += a[i].h.len;
			if (a[i].h.getSpan() >= Constraints.GENE_MIN_LEN && a[i].getScore() >= Constraints.PAIRED_MIN_RATIO)
				if (!a[i].h.isFAT())
					res.add(a[i]);
		}
		return res;
	}

	public static void main(String[] args) throws IOException {
		String gbk = GbkReader.read("c:\\yandex\\Tests\\gbk_for_students\\ref_chr7_00.gbk");
		GlobalContext.init(gbk);
		int total = 0;
		List<State> all = new ArrayList<DynamicSolver.State>();
		for (int i = 0; i < gbk.length(); i+=200) {
			//if (i%10000 == 0) System.err.println(i);
			int[] codes = Utils.toInt(gbk.substring(i,Math.min(i+400,gbk.length())));
			List<Helix> helixes = HelixFinder.findHelixesUnefficient(codes, 4, 4, 1000000);
			Helix.addShift(helixes, i);
			//System.out.println(helixes.size());
			DynamicSolver solver = new DynamicSolver(helixes);
			List<State> res = solver.solve();
			total += res.size();
			all.addAll(res);
		}
		Collections.sort(all, new Comparator<State>() {
			@Override
			public int compare(State arg0, State arg1) {
				return arg0.h.start - arg1.h.start; 
				//return Double.compare(arg1.getScore(), arg0.getScore());
			}
		});
		System.out.println(total);
		for (int i = 0; i < all.size() && i < 100; ++i)
			System.out.println(all.get(i) + " fat " + all.get(i).h.getFAT());
	}
}
