package io;

import java.io.IOException;
import java.util.List;

import model.CGene;
import model.FamilyCriteria;
import model.GlobalContext;
import utils.Utils;
import algo.GeneFinder;
import algo.GenePieceFinder;
import algo.GeneSelection;

public class Main {
	
	static GeneFinder geneFinder = new GenePieceFinder();
	
	public static void main(String[] args) throws IOException {
		if (args.length < 1 || args.length > 2) {
			System.err.println("Usage: java -jar rna.jar <gbk> [<output>]");
			return;
		}
		String inputFile = args[0];
		String outputFile = args.length >= 2 ? args[1] : inputFile.substring(0, inputFile.lastIndexOf('.')) + ".ans";
		
		String gbk = GbkReader.read(inputFile);
		GlobalContext.init(gbk);
		Utils.log("data loaded");
		
		List<CGene> all = geneFinder.findGenes(gbk);
		Utils.log("c-genes search done, "+all.size()+" c-genes found");

		List<CGene> filtered = GeneSelection.selectOptimal(all, new FamilyCriteria.Weight());
		Utils.log("optimal blocks selected, "+filtered.size()+" in total");

		SolutionWriter.write(filtered, outputFile);
		
		Utils.log("output written");	
	}

}
