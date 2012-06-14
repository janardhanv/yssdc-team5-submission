package model;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.StringTokenizer;

import utils.Utils;

public class Helix {
	// (start,end) - index of innermost pair in the helix
	public int start;
	public int end;
	public int len;
	public Helix(int start, int end, int len) {
		super();
		this.start = start;
		this.end = end;
		this.len = len;
	}
	public int getDist() {
		return end-start-1;
	}
	public int getLeft() {
		return start-len+1;
	}
	public int getRight() {
		return end+len-1;
	}
	public int getSpan() {
		return end-start-1+2*len;
	}
	@Override
	public String toString() {
		return GlobalContext.gbk.substring(start-len+1,start+1)+"@"+start+"..."+
			GlobalContext.gbk.substring(end,end+len)+"@"+end+" len="+len+" dist="+getDist();
	}
	public void addShift(int shift) {
		start += shift;
		end += shift;
	}
	public static void addShift(Collection<Helix> a, int shift) {
		for (Helix helix : a) {
			helix.addShift(shift);
		}
	}
	
	public static String serialize(Helix helix) {
		return helix.start + " " + helix.end + " " + helix.len;
	}
	
	public static String serialize(List<Helix> helixes) {
		String res = "";
		for (int i = 0; i < helixes.size(); ++i) {
			if (i > 0) res += "#";
			res += serialize(helixes.get(i));
		}
		return res;
	}
	
	public static Helix deserialize(String s) {
		//Utils.log("desering helix: " + s);
		String[] ss = s.split(" ");
		int start = Integer.parseInt(ss[0]);
		int end = Integer.parseInt(ss[1]);
		int len = Integer.parseInt(ss[2]);
		return new Helix(start, end, len);
	}
	
	public static List<Helix> deserializeList(String s) {
		StringTokenizer st = new StringTokenizer(s, "#");
		List<Helix> res = new ArrayList<Helix>();
		while (st.hasMoreTokens()) {
			Helix helix = deserialize(st.nextToken());
			//Utils.log("after deser:" + helix);
			res.add(helix);
		}
		return res;
	}
	
	
	public boolean isHairpin() {
		return Constraints.HAIRPIN_MIN_LEN <= getDist() && getDist() <= Constraints.HAIRPIN_MAX_LEN; 
	}
	// Uses global context!!!
	public double getFAT() {
		return Utils.getFAT(GlobalContext.codes, getLeft(), getRight());
	}
	public boolean isFAT() {
		return Utils.isFAT(GlobalContext.codes, getLeft(), getRight());
	}
}
