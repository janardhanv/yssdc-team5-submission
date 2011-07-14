package model;

import java.util.Collection;

import utils.Utils;

public class Helix {
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
