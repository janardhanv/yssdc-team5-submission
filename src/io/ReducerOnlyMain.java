package io;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import model.CGene;
import model.FamilyCriteria;
import utils.Utils;
import algo.GeneFinder;
import algo.GenePieceFinder;
import algo.GeneSelection;

public class ReducerOnlyMain {
	
	static GeneFinder geneFinder = new GenePieceFinder();
	
	public static void main(String[] args) throws IOException {
		if (args.length < 1 || args.length > 2) {
			System.err.println("Usage: java -jar rna.jar <gbk> [<output>]");
			return;
		}
		String inputFile = args[0];
		String outputFile = args.length >= 2 ? args[1] : inputFile.substring(0, inputFile.lastIndexOf('.')) + ".ans";
		
		BufferedReader in = new BufferedReader(new FileReader(inputFile));
		List<CGene> genes = new ArrayList<CGene>();
		for (String s; (s = in.readLine()) != null; ) {
			genes.add(CGene.deserialize(s));
		}
		Utils.log("data loaded, "+genes.size()+" c-genes");
		
		List<CGene> filtered = GeneSelection.selectOptimal(genes, new FamilyCriteria.Weight());
		Utils.log("optimal blocks selected, "+filtered.size()+" in total, "+CGene.totalPairs(filtered)+" pairs");

		SolutionWriter.write(filtered, outputFile);
		
		Utils.log("output written");
	}

}
