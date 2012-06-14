package algo;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

import utils.LowerBound;

import model.CGene;
import model.FamilyCriteria;

public class GeneSelection {
	public static class State implements LowerBound.IntKey {
		public CGene g;
		public State best;
		public State prev;
		public int answer;
		public State(CGene g) {
			super();
			this.g = g;
			this.answer = 0;
		}
		@Override
		public int getKey() {
			return g.maxPosition;
		} 
	}
	public static List<CGene> selectOptimal(List<CGene> genes, FamilyCriteria criteria) {
		List<State> a = new ArrayList<GeneSelection.State>();
		a.add(new State(CGene.dummy));
		for (CGene g : genes) {
			a.add(new State(g));
		}
		Collections.sort(a, new Comparator<State>() {
			@Override
			public int compare(State o1, State o2) {
				return o1.g.maxPosition - o2.g.maxPosition;
			}
		});
		for (int i = 1; i < a.size(); i++) {
			State s = a.get(i);
			State x = LowerBound.lastBelow(a.subList(0, i), a.get(i).g.minPosition);
			s.answer = x.answer + criteria.getWeight(s.g);
			if (a.get(i-1).answer > s.answer) {
				s.answer = a.get(i-1).answer;
				s.best = a.get(i-1).best;
			} else {
				s.best = s;
				s.prev = x.best;
			}
		}
		List<CGene> res = new ArrayList<CGene>();
		for (State s = a.get(a.size()-1).best; s != null; s = s.prev) {
			res.add(s.g);
		}
		return res;
	}
}
