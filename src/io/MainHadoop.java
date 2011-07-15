package io;


import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import javax.rmi.CORBA.Util;

import model.CGene;
import model.FamilyCriteria;
import model.GlobalContext;

import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.FileSystem;
import org.apache.hadoop.fs.Path;
import org.apache.hadoop.io.LongWritable;
import org.apache.hadoop.io.Text;
import org.apache.hadoop.mapreduce.Job;
import org.apache.hadoop.mapreduce.Mapper;
import org.apache.hadoop.mapreduce.Reducer;
import org.apache.hadoop.mapreduce.lib.input.FileInputFormat;
import org.apache.hadoop.mapreduce.lib.output.FileOutputFormat;
import org.apache.hadoop.util.GenericOptionsParser;

import utils.Utils;
import algo.GeneFinder;
import algo.GenePieceFinder;
import algo.GeneSelection;

public class MainHadoop {
	private final static Text OUT_KEY = new Text(new String("cgene"));
	
	public static class MainHadoopMapper extends
			Mapper<LongWritable, Text, Text, Text> {
		
		static GeneFinder geneFinder = new GenePieceFinder();
		
		// key - offset
		// value - substring of rna
		public void map(LongWritable key, Text value, Context context)
				throws IOException, InterruptedException {
			
			//Utils.log("get " + key );
			String gbk = value.toString().toUpperCase();
			GlobalContext.init(gbk);
			
			List<CGene> genes = geneFinder.findGenes(gbk);
			//Utils.log("c-genes search done, "+genes.size()+" c-genes found");
			for (CGene gene : genes) {
				//Utils.log("shifting " + gene);
				gene.addShift((int) key.get());
				//Utils.log("after shifting: " + gene);
				Text outValue = new Text(CGene.serialize(gene));
				context.write(OUT_KEY, outValue);
			}
		}
	}

	public static class MainHadoopReducer extends
			Reducer<Text, Text, Text, Text> {
		
		public void reduce(Text key, Iterable<Text> values, Context context)
				throws IOException, InterruptedException {
			//Utils.log("reduce: " + key);
			//Utils.log("value: " + values);
			List<CGene> all = new ArrayList<CGene>();
			for (Text value : values) {
				//Utils.log("value: " + value);
				CGene gene = CGene.deserialize(value.toString());
				//Utils.log("Cgene: " + gene);
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
	    job.setOutputValueClass(Text.class);
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
