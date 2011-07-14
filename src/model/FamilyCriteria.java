package model;

public interface FamilyCriteria {
	public int getWeight(CGene gene);
	
	public static class Length implements FamilyCriteria {
		@Override
		public int getWeight(CGene gene) {
			return 1;
		}
	}
	
	public static class Weight implements FamilyCriteria {
		@Override
		public int getWeight(CGene gene) {
			return gene.pairs;
		}
	}
}
