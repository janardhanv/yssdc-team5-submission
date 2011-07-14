package model;

import java.util.Collections;
import java.util.Comparator;
import java.util.List;

import algo.DynamicSolver.State;

public class CGene {

	public List<Helix> helixes;
	int minPosition;
	int maxPosition;
	int score;
	
	
	
	public CGene(List<Helix> helixes) {
		super();
	}

	@Override
	public String toString() {
		Collections.sort(helixes, new Comparator<Helix>() {
			@Override
			public int compare(Helix arg0, Helix arg1) {
				return arg0.start - arg1.start; 
			}
		});

	}

	public static void main(String[] args) {
		// TODO Auto-generated method stub

	}

}
