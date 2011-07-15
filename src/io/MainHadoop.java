package io;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.StringTokenizer;

import model.CGene;
import model.FamilyCriteria;

import org.apache.hadoop.io.LongWritable;
import org.apache.hadoop.io.Text;
import org.apache.hadoop.mapreduce.Mapper;
import org.apache.hadoop.mapreduce.Reducer;
import org.apache.hadoop.mapreduce.Mapper.Context;

import algo.GeneFinder;
import algo.GenePieceFinder;
import algo.GeneSelection;

public class MainHadoop {
	
	public static class MainHadoopMapper extends
			Mapper<Integer, Text, Text, Text> {
		private final static Text OUT_KEY = new Text(new String("cgene"));
		static GeneFinder geneFinder = new GenePieceFinder();
		
		// key - offset
		// value - substring of rna
		public void map(Integer key, Text value, Context context)
				throws IOException, InterruptedException {
			List<CGene> genes = geneFinder.findGenes(value.toString());
			for (CGene gene : genes) {
				gene.addShift(key);
				Text outValue = new Text(CGene.serialize(gene));
				context.write(OUT_KEY, outValue);
			}
		}
	}

	public static class MainHadoopReducer extends
			Reducer<Text, LongWritable, Text, LongWritable> {
		
		public void reduce(Text key, Iterable<Text> values, Context context)
				throws IOException, InterruptedException {
				
			List<CGene> all = new ArrayList<CGene>();
			for (Text value : values) {
				CGene gene = CGene.deserialize(value.toString());
				all.add(gene);
				
			}
			List<CGene> filtered = GeneSelection.selectOptimal(all, new FamilyCriteria.Weight());
		}
	}	
	
	
}
