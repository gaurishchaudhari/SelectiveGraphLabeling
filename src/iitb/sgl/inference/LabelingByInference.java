package iitb.sgl.inference;

import iitb.shared.gm.inference.ScalableMessagePassing.Method;
import iitb.shared.gm.inference.SparseScalableMP;
import iitb.shared.ArrayUtils;
import iitb.shared.gm.PottsPotential;
import iitb.shared.gm.UDGM;

import java.io.FileWriter;
import java.util.Random;

import iitb.sgl.data.EdgeFeatureGenerator;
import iitb.sgl.data.NodeFeatureGenerator;
import iitb.sgl.data.NodeFeatured;
import iitb.sgl.data.PokecGraph;
import iitb.sgl.data.SocialGraph;

/**
 * Run an Inference Algorithm for predicting locations + Estimating confidence
 * using marginals.
 */

public class LabelingByInference {
	String precisionRecallDataFile = "";

	SocialGraph snGraph;
	UDGM model;
	EdgeFeatureGenerator edgeFeatureGenerator;
	NodeFeatureGenerator nodeFeatureGenerator;
	int numNodeFeatures;
	int numEdgeFeatures;
	int numFeatures;

	public LabelingByInference(SocialGraph snGraph, int perKnown, String precisionRecallDataFile) {
		this.snGraph = snGraph;
		this.nodeFeatureGenerator = null;
		this.numNodeFeatures = 0;
		if (snGraph instanceof NodeFeatured) {
			this.nodeFeatureGenerator = new NodeFeatureGenerator(snGraph);
			this.numNodeFeatures = nodeFeatureGenerator.numFeatures();
		}
		this.edgeFeatureGenerator = new EdgeFeatureGenerator(snGraph);
		this.numEdgeFeatures = edgeFeatureGenerator.numFeatures();
		this.numFeatures = numNodeFeatures + numEdgeFeatures;
		
		this.precisionRecallDataFile = precisionRecallDataFile;
		
		System.out.println(numNodeFeatures+" "+numEdgeFeatures+" "+numFeatures );
	}

	public void runInference(double[] weights, boolean isUniformPotts, Method method) throws Exception {
		
		this.model = new UDGM(snGraph.graph, snGraph.numLabels, new PottsPotential(0.5, 0), 1);
		if(isUniformPotts)
			// Set only node potentials
			SGModelUtils.setNodePotentials(snGraph, model, weights, nodeFeatureGenerator, edgeFeatureGenerator);
		else{
			System.out.println("Using Weights " + ArrayUtils.toStrArr(weights));
			SGModelUtils.setPotentials(snGraph, model, weights, nodeFeatureGenerator, edgeFeatureGenerator);
		}
		
		SparseScalableMP solution = new SparseScalableMP(model, 1, method, false, snGraph.observedNodes);
		solution.getTopK(1);

		FileWriter fw = new FileWriter(precisionRecallDataFile);
		int total = 0;
		int correct = 0;
		int recall = 0;
		double softCount = 0;
		Random r = new Random(5);
		for (int user = snGraph.graph.getNumNodes() - 1; user >= 0; --user) {
			int trueIndex = snGraph.nodeLabels.get(user);
			int predictedIndex = solution.getSolution(user);
			if (snGraph.isNodeObserved(user)) {
				assert (trueIndex == predictedIndex);
				continue;
			}
			double confidence = solution.getConfidence(user);
			assert (confidence > 0);
			if (predictedIndex == -1)
				predictedIndex = r.nextInt(snGraph.numLabels);
			else
				++recall;
			++total;
			
			double corr = 1 - Math.abs(trueIndex - predictedIndex) / (snGraph.numLabels * 1.0);
			softCount += corr;
			if (trueIndex == predictedIndex)
				++correct;
			fw.write(confidence + " " + (trueIndex == predictedIndex ? 1 : 0) + " " + user+" "+trueIndex+" "+predictedIndex + "\n");
		}
		fw.close();
		System.out.println("Test Set Accuracy : " + correct + "/" + total + " = " + (correct * 100.0 / total));
		System.out.println("Test Set SOFT-Accuracy : " + softCount + "/" + total + " = " + (softCount * 100.0 / total));
		double precision = correct*100.0/recall;
		System.out.println("Recall = " + recall + " Precision = "+precision+" F-Measure = "+(2*precision*recall/(precision+recall)));
	}

	public static void main(String[] args) throws Exception {
		String inputPath = "/mnt/bag/wwt/twitter/models/gModel_Pokec/";
		String outputPath = "/mnt/bag/wwt/twitter/pubPlots/temp/";
		int perKnown = 2;
		int rSeed = 2;
		String method= "Cvx_JointL_EM";
		LabelingByInference inference = new LabelingByInference(new PokecGraph(inputPath, perKnown, rSeed), perKnown, outputPath + method + "_" + perKnown+"_r" + rSeed);
		
		double[] weights = new double[]{0.0034, 0.0024, 0.0016, 0.0013, 0.0042, 0.0016, 0.0017, 0.0046, 0.0099};
		inference.runInference(weights, method.equals("Uniform")? true : false, method.equals("SoftICA") ? Method.MeanField : Method.BP);
		System.out.println("Done !");
	}
}