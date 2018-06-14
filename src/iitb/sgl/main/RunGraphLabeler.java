package iitb.sgl.main;

import iitb.shared.gm.inference.ScalableMessagePassing.Method;
import iitb.shared.ArrayUtils;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.Properties;
import java.util.Random;
import java.util.StringTokenizer;

import iitb.utils.Stopwatch;

import iitb.sgl.data.SocialGraph;
import iitb.sgl.data.PokecGraph;
import iitb.sgl.data.TwitterGraph;
import iitb.sgl.inference.LabelingByInference;
import iitb.sgl.learning.Learner;
import iitb.sgl.learning.JL;
import iitb.sgl.learning.NCL;
import sun.reflect.generics.reflectiveObjects.NotImplementedException;

/**
 * MAIN CLASS OF THIS PROJECT.
 * Run a learner followed by an inference through this class.
 */

public class RunGraphLabeler {

	enum LearningAlgo { 
		JL, 	// Joint Likelihood training
		NCL 	// Node Conditional Likelihood training
	};
	
	String inputPath;		// Input path of Social Graph
	String outputPath;		// Directory to put output files
	String learningAlgo;	// Choice of Learning algorithm
	int percentObserved;	// Percent of Nodes observed in Social Graph (e.g. 2%, 5%, 10%)
	int rSeed;				// Random seed used to choose observed nodes in social graph
	SocialGraph snGraph;	// Social Graph (can be TwitterGraph or PokecGraph)
	
	static boolean isTwitter = false;
	boolean shouldLearn = true;		// bool indicating should optimize weights or not
	double[] initWeights = null;	// Initialized weights
	String initWeightsFile = "/mnt/a99/d0/gaurish/config/ObsNbr_Weights.log";
	
	public static void main(String[] args) throws Exception {
		
		String inputPath = isTwitter 
				? "/mnt/bag/wwt/twitter/models/gModel_Twitter/subGraph/" 
				: "/mnt/bag/wwt/twitter/models/gModel_Pokec/subGraph/";
		
		String outputPath = "/mnt/bag/wwt/twitter/outputs/";	
	
		int percentObserved = 2;
		int seed = 2;
		double C = 1E-8;

		String algo = "NCL";
		
		new RunGraphLabeler(inputPath, outputPath, percentObserved, seed, algo, C);
	}
		
	public RunGraphLabeler(String inputPath, String outputPath, int percentObserved, int rSeed, String learningAlgo, double C) throws Exception{
		
		this.inputPath = inputPath;
		this.outputPath = outputPath;
		this.percentObserved = percentObserved;
		this.rSeed = rSeed;
		this.learningAlgo = learningAlgo;
		
		System.out.println("Running - " + learningAlgo + ", % Observed = " + percentObserved + ", Seed = " + rSeed);
		
		System.out.println("Loading Social Graph...");
		if(isTwitter)
			this.snGraph = new TwitterGraph(inputPath, percentObserved, rSeed);
		else
			this.snGraph = new PokecGraph(inputPath, percentObserved, rSeed);
		
		double[] weights = null;
		LearningAlgo algo = LearningAlgo.valueOf(learningAlgo);
		Stopwatch sw = new Stopwatch();

		if(shouldLearn){
			System.out.println("Running Learner...");
			initializeWeights();
			Learner learner = null;
			
			if(algo == LearningAlgo.JL){
				learner = new JL(snGraph, new Properties(), new Random(1), 1, initWeights, C);
			}
			else if(algo == LearningAlgo.NCL){
				learner = new NCL(snGraph, new Properties(), new Random(1), 1, false, initWeights, C);
			}
			else{
				System.out.println("Invalid Learning Algorithm! " + learningAlgo);
				throw new NotImplementedException();
			}
			learner.learnWeights();
			weights = learner.getWeights();
		}
		else{
			weights = loadWeightsFromFile("/mnt/bag/wwt/twitter/RunGraphLabeler_" + (isTwitter? "Twitter" : "Pokec") + "/" + learningAlgo + "_" + percentObserved + "_" + rSeed + ".log");
		}
		
		if(weights != null) {
			System.out.println("Optimal Weights = " + ArrayUtils.toStrArr(weights));
		}
		System.out.println("Learning Execution Time " + sw.elapsedTimes());
		
		sw.restart();
		System.out.println("Running Inference...");
		
		if(algo == LearningAlgo.JL || algo == LearningAlgo.NCL){
			String preRclDataFile = outputPath + learningAlgo + "_" + percentObserved+"_r" + rSeed;
			LabelingByInference inference = new LabelingByInference(snGraph, percentObserved, preRclDataFile);
			inference.runInference(weights, false, Method.BP);	
		}

		System.out.println("Total Execution Time " + sw.elapsedTimes());
		
		// TOOD: Run ConfidenceModel for NCL to add confidence measurement to each label prediction.
	}
	
	public void initializeWeights() throws Exception{
		if(!new File(initWeightsFile).exists()) {
			return;		// if initWeights = null, Learner will initialize weights to some default value
		}
		BufferedReader br = new BufferedReader(new FileReader(initWeightsFile));
		String line;
		while((line=br.readLine()) != null){
			if(line.contains("Optimal Weights")){
				String arrStr = line.substring(line.indexOf("[")+1, line.indexOf("]"));
				StringTokenizer st = new StringTokenizer(arrStr, ",");
				initWeights = new double[st.countTokens()];
				for(int i=0; st.hasMoreTokens(); ++i)
					initWeights[i] = Double.parseDouble(st.nextToken());
			}
		}
		br.close();
	}
	
	public static double[] loadWeightsFromFile(String logFile) throws Exception{
		double[] weights = null;
		BufferedReader br = new BufferedReader(new FileReader(logFile));
		String line;
		while((line = br.readLine()) != null){
			if(line.contains("Optimal Weights")){
				line = line.substring(line.indexOf("[")+1, line.indexOf("]"));
				StringTokenizer st = new StringTokenizer(line,",");
				weights = new double[st.countTokens()];
				for(int i = 0; st.hasMoreTokens(); ++i)
					weights[i] = Double.parseDouble(st.nextToken());
				break;
			}
		}
		br.close();
		return weights;
	}
}

