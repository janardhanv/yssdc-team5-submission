package io;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;
import java.util.StringTokenizer;

import model.CGene;
import model.FamilyCriteria;

import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.FileSystem;
import org.apache.hadoop.fs.Path;
import org.apache.hadoop.fs.permission.FsAction;
import org.apache.hadoop.io.Text;
import org.apache.hadoop.util.GenericOptionsParser;

import utils.Utils;

import algo.GeneSelection;

public class Hadoop {


	private static void analyzeResult(Path file, String outputFile) throws IOException {
		FileSystem fs = FileSystem.get(new Configuration());
		if (!fs.exists(file))
			return;

		BufferedReader in = new BufferedReader(new InputStreamReader(fs
				.open(file)));
		
		List<CGene> all = new ArrayList<CGene>();
		int totalScore = 0;
		while (in != null && in.ready()) {
			String s = in.readLine();
			//Utils.log("analazing string " +  s);
			StringTokenizer st = new StringTokenizer(s, "\t");
			st.nextToken(); // read cgene
			String value = st.nextToken();
			CGene gene = CGene.deserialize(value);
			totalScore += gene.pairs;
			Utils.log("found " + gene);
			all.add(gene);
		}

		Utils.log("optimal blocks selected, " + all.size() + " in total");
		Utils.log("total score: " + totalScore);
		SolutionWriter.write(all, outputFile);
	}
	
	public static void reduce(Path inputFile, Path outputFile)  throws IOException {
		FileSystem fs = FileSystem.get(new Configuration());
		if (!fs.exists(inputFile))
			return;
		
		BufferedReader in = new BufferedReader(new InputStreamReader(fs
				.open(inputFile)));
		
		List<CGene> all = new ArrayList<CGene>();
		while (in != null && in.ready()) {
			String s = in.readLine();
			//Utils.log("analazing string " +  s);
			StringTokenizer st = new StringTokenizer(s, "\t");
			st.nextToken(); // read cgene
			String value = st.nextToken();
			//Utils.log("value: " + value);
			CGene gene = CGene.deserialize(value.toString());
			//Utils.log("Cgene: " + gene);
			all.add(gene);
			
		}
		BufferedWriter out = new BufferedWriter(new OutputStreamWriter(fs.create(outputFile, true)));
		List<CGene> filtered = GeneSelection.selectOptimal(all, new FamilyCriteria.Weight());
		for (CGene cGene : filtered) {
			String outValue = CGene.serialize(cGene);
			System.out.println("writing to file " + cGene);
			out.write("cgene" + "\t" + outValue + "\n");
		}
		out.flush();
		out.close();
	}
	
	// program to reduce input of mappers and analyze result
	// genes a stored in a form
	// cgene \t helix # helix ...
	public static void main(String[] args) throws Exception{
		Configuration conf = new Configuration();
	    String[] otherArgs = new GenericOptionsParser(conf, args).getRemainingArgs();
	    if (otherArgs.length != 3) {
	      System.err.println("Usage: [this] <in> <out> <outputFile>");
	      System.exit(2);
	    }
	    
	    Path inputFile = new Path(otherArgs[0]);
	    Path outputFile = new Path(otherArgs[1]);
	    reduce(inputFile, outputFile);
	    analyzeResult(outputFile, otherArgs[2]);
	}

}
