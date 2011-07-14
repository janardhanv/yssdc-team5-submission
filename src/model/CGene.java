package model;

import java.util.Collections;
import java.util.Comparator;
import java.util.List;

import algo.DynamicSolver.State;

public class CGene {

	public List<Helix> helixes;
	int minPosition;
	int maxPosition;
	int pairs;
	
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
		return "Gene "+getScore()+" "+minPosition+"-"+maxPosition;
	}

	private double getScore() {
		return pairs*2.0/(maxPosition - minPosition);
	}

}
