package io;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;

public class GbkReader {
	public static String read(String filename) throws IOException {
		BufferedReader r = new BufferedReader(new FileReader(filename));
		int len = Integer.parseInt(r.readLine());
		String s = r.readLine().toUpperCase().replace(" ", "");
		if (s.length() != len)
			throw new IOException("Bad GBK file!");
		return s;
	}
}
