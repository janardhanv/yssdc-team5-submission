package model;

import java.util.List;

public class CGene {

	public List<Helix> helixes;
	public int minPosition;
	public int maxPosition;
	public int pairs;
	
	private CGene(int fix) {
		minPosition = maxPosition = -1;
	}
	
	public static CGene dummy = new CGene(-1);
	
	public CGene(List<Helix> helixes) {
		super();
		this.helixes = helixes;
		minPosition = helixes.get(0).getLeft();
		maxPosition = helixes.get(0).getRight();
		pairs = 0;
		for (Helix helix : helixes) {
			minPosition = Math.min(minPosition, helixes.get(0).getLeft());
			maxPosition = Math.max(maxPosition, helixes.get(0).getRight());
			pairs += helix.len;
		}
	}

	@Override
	public String toString() {
		return "Gene "+getScore()+" "+minPosition+"-"+maxPosition+" pairs "+pairs;
	}

	public double getScore() {
		return pairs*2.0/(maxPosition - minPosition + 1);
	}

}
