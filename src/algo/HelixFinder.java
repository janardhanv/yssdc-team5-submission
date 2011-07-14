package algo;
import java.util.ArrayList;
import java.util.List;

import model.Helix;
import utils.Utils;


public class HelixFinder {
	public static List<Helix> findHelixesClosestOnly(int[] a, int len, int mindist, int maxdist) {
		List<Helix> res = new ArrayList<Helix>();
		int[] fwd = Utils.fwdCodes(a, len);
		int[] back = Utils.backInvCodes(a, len);
		int pw = 1;
		for (int i = 0; i < len; ++i)
			pw *= 4;
		Integer[] prev = new Integer[pw];  
		for (int i = 0; i < a.length-len+1; ++i) {
			if (i-mindist-len >= 0 && fwd[i-mindist-len] != -1)
				prev[fwd[i-mindist-len]] = i-mindist-len;
			if (i-maxdist-1-len >= 0 && fwd[i-maxdist-1-len] != -1 &&
				prev[fwd[i-maxdist-1-len]].equals(i-maxdist-1-len))
				prev[fwd[i-maxdist-1-len]] = null;
			if (back[i] != -1 && prev[back[i]] != null) {
				Helix h = new Helix(prev[back[i]]+len-1, i, len);
				if (h.getDist() >= mindist && (a[h.start+1]^1) == a[h.end-1])
					continue;
				while (h.len <= h.start && h.end+h.len < a.length && (a[h.start-h.len]^1) == a[h.end+h.len])
					++h.len;
				res.add(h);
			}
		}
		return res;
	}
	
	public static List<Helix> findHelixesUnefficient(int[] a, int len, int mindist, int maxdist) {
		List<Helix> res = new ArrayList<Helix>();
		int[] fwd = Utils.fwdCodes(a, len);
		int[] back = Utils.backInvCodes(a, len);
		for (int i = 0; i < a.length-len+1; ++i) {
			for (int j = i+len+mindist; j < a.length-len+1; ++j) {
				if (fwd[i] == back[j] && fwd[i]!=-1) {
					Helix h = new Helix(i+len-1, j, len);
					if (h.getDist() > maxdist) continue;
					if (h.getDist() >= mindist && (a[h.start+1]^1) == a[h.end-1]) continue;
					while (h.len <= h.start && h.end+h.len < a.length && (a[h.start-h.len]^1) == a[h.end+h.len])
						++h.len;
					res.add(h);
				}
			}
		}
		return res;
	}	
}
