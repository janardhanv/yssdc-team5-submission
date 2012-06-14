package utils;

import io.GbkReader;

import java.io.IOException;

import model.Constraints;
import model.GlobalContext;
import algo.HelixFinder;

public class StatsCalculator {
	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		//args[0] = "c:\\yandex\\Tests\\gbk_for_students\\ref_chr11_01.gbk";
		String gbk = GbkReader.read(args[0]);
		//gbk = gbk.substring(0,5000000);
		//gbk = gbk.substring(2636000,2636200+1000);
		gbk = gbk.substring(2638200,2638200+1000);
		System.out.println(gbk);
		GlobalContext.gbk = gbk;
		
		int[] codes = Utils.toInt(gbk);
		int[] freq = new int[5];
		for (int i = 0; i < codes.length; i++) {
			int c = codes[i];
			if (c == -1) c = 4;
			++freq[c];
		}
		for (int i = 0; i <= 4; i++) {
			System.out.println(String.format("%c %d", (Utils.ALPHABET+"N").charAt(i), freq[i]));
		}
		
		System.out.println("Found helixes "+HelixFinder.findHelixesClosestOnly(codes, 4, 4, 20).size());
		/*for (Helix helix : h) {
			System.out.println(helix);
		}*/
		for (int len = 4; len <= 20; ++len) {
			System.out.println("Of size "+len+" "+HelixFinder.findHelixesClosestOnly(codes, Constraints.HELIX_MIN_LEN, len, len).size());
		}
		
		int mx = 0;
		int shift = 0;
		for (int b = 0; b+1000 <= gbk.length(); b+=500) {
			String piece = gbk.substring(b,b+1000);
			int[] pieceCodes = Utils.toInt(piece);
			int x = HelixFinder.findHelixesUnefficient(pieceCodes, 4, 4, 10000).size();
			mx = Math.max(x, mx);
			if (mx == x)
				shift = b;
			if (x > 4000) System.out.println("@ "+b+" "+x+" fat="+Utils.getFAT(pieceCodes)+" fat="+Utils.getFAT(pieceCodes));
		}
		System.out.println("Max per block: "+mx+" @ "+shift);
	}

}
