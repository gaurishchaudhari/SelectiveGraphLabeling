package iitb.sgl.learning;

import java.util.Arrays;
import java.util.BitSet;
import java.util.Properties;
import java.util.Random;

import iitb.sgl.data.NodeFeatureGenerator;
import iitb.sgl.data.EdgeFeatureGenerator;
import iitb.sgl.data.SocialGraph;
import iitb.sgl.data.NodeFeatured;
import iitb.shared.optimization.ConvexOptimization;
import iitb.shared.optimization.GradientCompute;
import iitb.shared.optimization.LBFGSTrainer;
import sun.reflect.generics.reflectiveObjects.NotImplementedException;

public class LearnerImpl implements GradientCompute, Learner {
	
	protected SocialGraph snGraph;
	protected Properties options;
	protected Random random;
	protected int debugLvl;
	
	protected NodeFeatureGenerator nodeFeatureGenerator;
	protected EdgeFeatureGenerator edgeFeatureGenerator;
	protected int numNodeFeatures;
	protected int numEdgeFeatures;
	protected int numFeatures;
	protected int[] nodeFeatureIds;
	protected int[] nodeFeatureLabels;
	
	protected BitSet candidateLabels;
	int[] labeling;
	
	protected double[] weights;					// parameters to be learnt
	protected ConvexOptimization trainer;		// optimizer
	
	public double INIT_W = 0.1;					// Initial value for all W's 
	public double C_VALUE = 30;					// Regularization parameter, C
	public static double GRAD_DELTA = 10e-5;	// Delta value to terminate optimization loop
	public static int MAX_ITERS = 100;			// Maximum number of optimization loop iterations
	
	public LearnerImpl(SocialGraph snGraph, Properties options, Random random, int debugLvl, double[] initWeights, double C) throws Exception{
		this(snGraph, options, random, debugLvl);
		
		weights = new double[this.numFeatures];
		if(initWeights != null){
			assert(initWeights.length == this.numFeatures);
			System.arraycopy(initWeights, 0, weights, 0, initWeights.length);
		}
		else
			Arrays.fill(weights, INIT_W);
		this.C_VALUE = C;
	}
	
	private LearnerImpl(SocialGraph snGraph, Properties options, Random random, int debugLvl) throws Exception {
		this.snGraph = snGraph;
		this.options = options;
		this.random = random; 
		this.debugLvl = debugLvl;
		this.trainer = new LBFGSTrainer(debugLvl,
				Integer.parseInt(options.getProperty("mForHessian","7")),
				Integer.parseInt(options.getProperty("maxIters","1000")),
				Double.parseDouble(options.getProperty("epsForConvergence","0.001")), true);		
		this.nodeFeatureGenerator = null;
		this.numNodeFeatures = 0;
		if(snGraph instanceof NodeFeatured){
			this.nodeFeatureGenerator = new NodeFeatureGenerator(snGraph);
			this.numNodeFeatures = nodeFeatureGenerator.numFeatures();
			this.nodeFeatureIds = new int[numNodeFeatures];
			this.nodeFeatureLabels = new int[nodeFeatureIds.length];
		}
		this.edgeFeatureGenerator = new EdgeFeatureGenerator(snGraph);
		this.numEdgeFeatures = edgeFeatureGenerator.numFeatures();
		this.numFeatures = numNodeFeatures + numEdgeFeatures;
		this.candidateLabels = new BitSet(snGraph.numLabels);
		System.out.println("# Features: " + numNodeFeatures+" + " + numEdgeFeatures + " = " + numFeatures);
	}
	
	public void initLabeling(){
		labeling = new int[snGraph.graph.getNumNodes()];
		for(int i=0; i<labeling.length; ++i){
			if(snGraph.isNodeObserved(i))
				labeling[i] = snGraph.nodeLabels.get(i);
			else
				labeling[i] = -1;
		}
	}
	
	@Override
	public void learnWeights() throws Exception {
		throw new NotImplementedException();
	}

	@Override
	public double computeFunctionGradient(double[] lambda, double[] grad) throws Exception {
		throw new NotImplementedException();
	}
	
	@Override
	public double[] getWeights(){
		return this.weights;
	}
}
