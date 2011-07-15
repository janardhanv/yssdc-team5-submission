package model;

import java.util.Collection;
import java.util.List;

import utils.Utils;

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
			minPosition = Math.min(minPosition, helix.getLeft());
			maxPosition = Math.max(maxPosition, helix.getRight());
			pairs += helix.len;
		}
	}
	
	public void addShift(int shift) {
		Helix.addShift(helixes, shift);
		minPosition += shift;
		maxPosition += shift;
	}
	
	public static void addShift(Collection<CGene> a, int shift) {
		for (CGene cGene : a) {
			cGene.addShift(shift);
		}
	}
	
	public static String serialize(CGene gene) {
		return Helix.serialize(gene.helixes);
	}
		
	public static CGene deserialize(String s) {
		//Utils.log("deserializing " + s);
		List<Helix> helixes = Helix.deserializeList(s);
		CGene gene = new CGene(helixes);
		//Utils.log("after deser" +  gene);
		return gene;
	}
		
	
	@Override
	public String toString() {
		return "Gene "+getScore()+" "+minPosition+"-"+maxPosition+" pairs "+pairs;
	}

	public double getScore() {
		return pairs*2.0/(maxPosition - minPosition + 1);
	}

	public static int totalPairs(List<CGene> a) {
		int total = 0;
		for (CGene cGene : a) {
			total += cGene.pairs;
		}
		return total;
	}

}
