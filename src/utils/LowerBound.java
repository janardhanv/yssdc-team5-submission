package utils;

import java.util.List;

public class LowerBound {
	public interface IntKey {
		public int getKey();
	}
	
	public static <T extends IntKey> T lastBelow(List<T> a, int key) {
		int left = -1, right = a.size()-1;
		while (left < right) {
			int mid = (left + right + 1) / 2;
			if (a.get(mid).getKey() < key)
				left = mid;
			else
				right = mid-1;
		}
		if (left == -1)
			return null;
		return a.get(left);
	}
	
	public static <T extends IntKey> int lastBelow(T[] a, int left, int right, int key) {
		while (left < right) {
			int mid = (left + right + 1) / 2;
			if (a[mid].getKey() < key)
				left = mid;
			else
				right = mid-1;
		}
		return left;
	}	
}
