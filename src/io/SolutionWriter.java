package io;

import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.List;

import model.Helix;

public class SolutionWriter {
	
	public static void write(List<Helix> a, String filename) throws IOException {
		PrintWriter out = new PrintWriter(new FileWriter(filename));
		out.close();
	}

	public static void main(String[] args) {
	}

}
