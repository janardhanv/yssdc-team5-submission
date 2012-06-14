package algo;
import io.GbkReader;
import io.SolutionWriter;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashSet;
import java.util.List;

import org.apache.commons.lang.time.StopWatch;

import model.CGene;
import model.Constraints;
import model.FamilyCriteria;
import model.GlobalContext;
import model.Helix;
import utils.LowerBound;
import utils.Utils;

public class DynamicSolver {
	public static final int GENE_LEN_MAX = 1000;
	
	public static class State implements LowerBound.IntKey {
		public Helix h;
		public int answer;
		public int cumulativeAnswer;
		public List<State> inner;
		// last included item
		public State best;
		// if included, previous one
		public State prev;
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
		public void collectHelixes(List<Helix> helixes) {
			helixes.add(h);
			for (State s: inner) {
				s.collectHelixes(helixes);
			}
		}
		public void writeDebug(String prefix) {
			Collections.sort(inner,new Comparator<State>() {
				@Override
				public int compare(State x, State y) {
					return x.h.getRight() < y.h.getRight() ? -1 : 1;
				}
			});
			System.out.println(prefix+h+" [");
			for (State state : inner) {
				state.writeDebug(prefix+"  ");
			}
			System.out.println(prefix+" ]");
		}
		@Override
		public int getKey() {
			return h.getRight();
		}
	}
	
	public State[] a;
	private int nStates;
	State dummy = new State(new Helix(-1,-1,0));
	{
		dummy.cumulativeAnswer = 0;
	}

	public DynamicSolver(List<Helix> h) {
		super();
		a = new State[h.size()];
		for (int i = 0; i < h.size(); ++i)
			a[i] = new State(h.get(i));
	}
	
	private State[] inside;
	private int nInside;
	
	private void solveInternal(State beginning, int startIndex, int startCoord, int endCoord) {
		nInside = 0;
		inside[nInside++] = beginning;
		for (int j = startIndex; j < nStates; j++) {
			if (a[j].h.getRight() > endCoord)
				break;
			if (startCoord > a[j].h.getLeft()/* || a[j].answer == 0*/)
				continue;
			int p = LowerBound.lastBelow(inside, 0, nInside-1, a[j].h.getLeft());
			a[j].cumulativeAnswer = inside[p].cumulativeAnswer + a[j].answer;
			if (inside[nInside-1].cumulativeAnswer > a[j].cumulativeAnswer) {
				a[j].cumulativeAnswer = inside[nInside-1].cumulativeAnswer;
				a[j].best = inside[nInside-1].best;
			} else {
				a[j].prev = inside[p].best;
				a[j].best = a[j];
			}
			inside[nInside++] = a[j];
		}
	}
	
	private List<State> getSolutionFromInside() {
		List<State> states = new ArrayList<DynamicSolver.State>();
		for (State s = inside[nInside-1].best; s != null; s = s.prev) {
			states.add(s);
		}
		return states;
	}
	
	public List<CGene> solve() {
		Arrays.sort(a,new Comparator<State>() {
			@Override
			public int compare(State x, State y) {
				return x.h.getRight() - y.h.getRight();
			}
		});
		List<CGene> res = new ArrayList<CGene>();
		inside = new State[a.length+1];
		nStates = 0;
		for (int i = 0; i < a.length; i++) {
			int startPos = LowerBound.lastBelow(a, 0, nStates-1, a[i].h.start+1)+1;
			solveInternal(dummy, startPos, a[i].h.start+1, a[i].h.end-1);
			a[i].inner = getSolutionFromInside();
			a[i].answer = inside[nInside-1].cumulativeAnswer;
			// there is something inside or we have a valid hairpin
			if (a[i].answer > 0 || a[i].h.isHairpin()) {
				a[i].answer += a[i].h.len;
				a[nStates++] = a[i];
			}
			/*if (a[i].h.getSpan() >= Constraints.GENE_MIN_LEN && a[i].getScore() >= Constraints.PAIRED_MIN_RATIO)
				if (!a[i].h.isFAT()) {
					List<Helix> helixes = new ArrayList<Helix>();
					a[i].collectHelixes(helixes);
					res.add(new CGene(helixes));
					//if (Math.abs(new CGene(helixes).getScore() - a[i].getScore()) > 0.000001) System.err.println("BAD!");
				}*/
		}
		HashSet<Integer> startPoses = new HashSet<Integer>();
		for (int i = 0; i < nStates; i++) {
			int left = a[i].h.getLeft();
			int searchStart = a[i].h.getRight()+1;
			//if (!startPoses.add(left)) continue;
			a[i].best = a[i];
			a[i].prev = null;
			a[i].cumulativeAnswer = a[i].answer;
			solveInternal(a[i], i, searchStart, searchStart+GENE_LEN_MAX);
			for (int j = 0; j < nInside; j++) {
				if (inside[j].best != inside[j]) continue;
				int span = inside[j].best.h.getRight() - left + 1;
				double score = 2.0 * inside[j].best.cumulativeAnswer / (double) span; 
				if (span >= Constraints.GENE_MIN_LEN && score >= Constraints.PAIRED_MIN_RATIO)
					if (!Utils.isFAT(GlobalContext.codes, left, left+span-1)) {
						List<Helix> helixes = new ArrayList<Helix>();
						for (State s = inside[j].best; s != null; s = s.prev) {
							s.collectHelixes(helixes);
						}
						res.add(new CGene(helixes));
//						if (Math.abs(new CGene(helixes).getScore() - score) > 0.000001)
						if (new CGene(helixes).getScore() > 1)System.err.println("BAD!");
					}
			}
		}
		return res;
	}

	public static void main(String[] args) throws IOException {
		StopWatch tm = new StopWatch();
		tm.start();
		//String gbk = GbkReader.read("c:\\yandex\\Tests\\gbk_for_students\\ref_chr7_00.gbk");
		String gbk = GbkReader.read("../ref_chr7_00.gbk");
		//String gbk = GbkReader.read("c:\\yandex\\test.txt");
		gbk = gbk.substring(235700, 236500+10000);
		GlobalContext.init(gbk);
		int total = 0;
		List<CGene> all = new ArrayList<CGene>();
		for (int i = 0; i < gbk.length(); i+=20000) {
			//if (i%10000 == 0) System.err.println(i);
			int[] codes = Utils.toInt(gbk.substring(i,Math.min(i+1500,gbk.length())));
			List<Helix> helixes = HelixFinder.findHelixesUnefficient(codes, 4, 4, 1000000);
			Helix.addShift(helixes, i);
			//System.out.println(helixes.size());
			DynamicSolver solver = new DynamicSolver(helixes);
			List<CGene> res = solver.solve();
			total += res.size();
			all.addAll(res);
		}
		Collections.sort(all, new Comparator<CGene>() {
			@Override
			public int compare(CGene arg0, CGene arg1) {
				//return arg0.minPosition - arg1.minPosition; 
				return Double.compare(arg1.getScore(), arg0.getScore());
			}
		});
		System.out.println(total);
		for (int i = 0; i < all.size() && i < 100; ++i)
			System.out.println(all.get(i));
		System.out.println("Elapsed "+tm);
		tm.reset();
		tm.start();

		System.out.println("Candidates " + all.size());
		List<CGene> filtered = GeneSelection.selectOptimal(all, new FamilyCriteria.Weight());
		System.out.println("After "+filtered.size());
		SolutionWriter.write(filtered, "../solution1.txt");
		
		System.out.println("Elapsed "+tm);
	}
}
