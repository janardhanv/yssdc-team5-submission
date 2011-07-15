package io;


import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.StringTokenizer;

import model.CGene;
import model.FamilyCriteria;

import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.FileSystem;
import org.apache.hadoop.fs.Path;
import org.apache.hadoop.io.LongWritable;
import org.apache.hadoop.io.Text;
import org.apache.hadoop.mapreduce.Job;
import org.apache.hadoop.mapreduce.Mapper;
import org.apache.hadoop.mapreduce.Reducer;
import org.apache.hadoop.mapreduce.Mapper.Context;
import org.apache.hadoop.mapreduce.lib.input.FileInputFormat;
import org.apache.hadoop.mapreduce.lib.output.FileOutputFormat;
import org.apache.hadoop.util.GenericOptionsParser;

import algo.GeneFinder;
import algo.GenePieceFinder;
import algo.GeneSelection;

public class MainHadoop {
	private final static Text OUT_KEY = new Text(new String("cgene"));
	
	public static class MainHadoopMapper extends
			Mapper<Integer, Text, Text, Text> {
		

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
			Reducer<Text, Text, Text, Text> {
		
		public void reduce(Text key, Iterable<Text> values, Context context)
				throws IOException, InterruptedException {
				
			List<CGene> all = new ArrayList<CGene>();
			for (Text value : values) {
				CGene gene = CGene.deserialize(value.toString());
				all.add(gene);
				
			}
			List<CGene> filtered = GeneSelection.selectOptimal(all, new FamilyCriteria.Weight());
			for (CGene cGene : filtered) {
				Text outValue = new Text(CGene.serialize(cGene));
				context.write(OUT_KEY, outValue);
			}
		}
	}	
	
	
	public static void main(String[] args) throws Exception {
	    Configuration conf = new Configuration();
	    String[] otherArgs = new GenericOptionsParser(conf, args).getRemainingArgs();
	    if (otherArgs.length != 2) {
	      System.err.println("Usage: [this] <in> <out>");
	      System.exit(2);
	    }

	    Job job = new Job(conf, "find genes");
	    job.setJarByClass(MainHadoop.class);
	    job.setMapperClass(MainHadoopMapper.class);
	    job.setReducerClass(MainHadoopReducer.class);
	    //job.setCombinerClass(MainHadoopMapper.class);
	    job.setOutputKeyClass(Text.class);
	    job.setOutputValueClass(LongWritable.class);
	    FileInputFormat.addInputPath(job, new Path(otherArgs[0]));
	    Path outputpath = new Path(otherArgs[1]);
	    FileSystem fs = FileSystem.get(new Configuration());
	    if(fs.exists(outputpath))
	      fs.delete(outputpath, true);
	    FileOutputFormat.setOutputPath(job, outputpath);
	    boolean result = job.waitForCompletion(true);
	    
	    System.exit(result ? 0 : 1);
	  }
	
}
