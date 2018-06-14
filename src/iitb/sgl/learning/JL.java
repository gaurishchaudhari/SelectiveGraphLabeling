package iitb.sgl.learning;

import java.util.Arrays;
import java.util.Properties;
import java.util.Random;

import iitb.sgl.data.SocialGraph;
import iitb.sgl.inference.SGInference;
import iitb.sgl.inference.SGModelUtils;
import iitb.shared.ArrayUtils;
import iitb.shared.gm.NodePotentialSparse;
import iitb.shared.gm.PottsPotential;
import iitb.shared.gm.UDGM;
import iitb.shared.gm.inference.ScalableMessagePassing.Method;

public class JL extends LearnerImpl {
	
	private UDGM model;
	private SGInference solution;
	private double[] marginal;
	private float[] expectedFeatures;
	
	private LambdaInference lmbdInference;
	
	private UDGM mStepModel;
	private SGInference mStepSolution;
	
	private UDGM tmpModel;
	private double[] tmpMarginal;
	
	public JL(SocialGraph snGraph, Properties options, Random random, int debugLvl, double[] initWeights, double C)
			throws Exception {
		super(snGraph, options, random, debugLvl, initWeights, C);
	}
	
	@Override
	public void learnWeights() throws Exception {
		System.out.println("Initial Weights: " + ArrayUtils.toStrArr(weights));	
		model = new UDGM(snGraph.graph, snGraph.numLabels, new PottsPotential(0.5, 0), 2);
		mStepModel = new UDGM(snGraph.graph, snGraph.numLabels, new PottsPotential(0.5, 0), 2);
		tmpModel = new UDGM(snGraph.graph, snGraph.numLabels, new PottsPotential(0.5, 0), 2);
		expectedFeatures = new float[this.numFeatures];
		marginal = new double[snGraph.numLabels];
		tmpMarginal = new double[snGraph.numLabels];
		lmbdInference = new LambdaInference(snGraph.graph, snGraph.numLabels, nodeFeatureGenerator, edgeFeatureGenerator);	
		MAX_ITERS = 100;
		double oldNorm = 0, newNorm = 0;
		int iter = 0;
		do{
			SGModelUtils.setPotentials(snGraph, this.model, this.weights, nodeFeatureGenerator, edgeFeatureGenerator);
			solution = new SGInference(model, 1, Method.BP, false, snGraph.observedNodes);
			solution.setIsDegreeInvGamma(true);
			solution.getTopK(1);
			Arrays.fill(expectedFeatures, 0);
			for(int node = snGraph.graph.getNumNodes() - 1; node >= 0; --node){
				if(nodeFeatureGenerator != null && !snGraph.observedNodes.get(node)){
					solution.getNodeSumMarginal(node, marginal, null);
					nodeFeatureGenerator.getFeatureValues(node, nodeFeatureIds, nodeFeatureLabels);
					for (int k = 0; k < nodeFeatureIds.length && nodeFeatureIds[k] >= 0; k++) {
						int fNum = nodeFeatureIds[k];
						int lbl = nodeFeatureLabels[k];
						expectedFeatures[fNum] += nodeFeatureGenerator.nFVal * marginal[lbl];
					}
				}
				for(int j=snGraph.graph.getNumNeighbours(node) - 1; j >= 0; --j){
					int nbr = snGraph.graph.getNeighbour(node, j);
					if(node < nbr){
						double mu = solution.getEdgeSameLabelProbability(node, nbr);
						for(int fNum = 0; fNum < numEdgeFeatures; ++fNum)
							expectedFeatures[numNodeFeatures + fNum] += mu * edgeFeatureGenerator.getFeatureValue(node, nbr, fNum);
					}
				}
			}
			oldNorm = ArrayUtils.norm2(weights);
			maximizationStep();
			newNorm = ArrayUtils.norm2(weights);
			if(debugLvl>0) System.out.println("E-Step " + iter+" Weights: "+ArrayUtils.toStrArr(weights)+" || "+newNorm);
			++iter;		
		}while(iter < MAX_ITERS && Math.abs(newNorm - oldNorm) > GRAD_DELTA);
	}
	
	public void maximizationStep() throws Exception{
		double oldNorm = 0, newNorm = 0;
		int iter = 0;
		do{					
			resetPotentials(mStepModel, weights);
			mStepSolution = new SGInference(mStepModel, 0, Method.BP, false, null);
			mStepSolution.setIsDegreeInvGamma(true);
			mStepSolution.setMaxIters(1);
			mStepSolution.getTopK(1);	
			lmbdInference.storeLambdaValues(mStepSolution);
			oldNorm = ArrayUtils.norm2(weights);
			trainer.optimize(weights, this);
			newNorm = ArrayUtils.norm2(weights);
			if(debugLvl>0) System.out.println("M-Step "+iter+" Weights: "+ArrayUtils.toStrArr(weights)+" || "+newNorm);
			++iter;		
		}while(iter < MAX_ITERS && Math.abs(newNorm - oldNorm) > GRAD_DELTA);
	}
	
	@Override
	public double computeFunctionGradient(double[] lambda, double[] grad) throws Exception {
		resetPotentials(tmpModel, lambda);
		double logZ = lmbdInference.getConvexifiedLogZ(tmpModel);
		double dataLL = ArrayUtils.dotProduct(lambda, expectedFeatures) - logZ;
		for(int i = 0; i < grad.length; ++i)
			grad[i] = -expectedFeatures[i];
		for(int node = snGraph.graph.getNumNodes()-1; node>=0; --node){
			if(nodeFeatureGenerator != null && !snGraph.observedNodes.get(node)){	// nothing is observed. so no && condition
				lmbdInference.getLambdaNodeMarginal(tmpModel, node, tmpMarginal, null);
				nodeFeatureGenerator.getFeatureValues(node, nodeFeatureIds, nodeFeatureLabels);
				for (int k = 0; k < nodeFeatureIds.length && nodeFeatureIds[k] != -1; k++) {
					int fNum = nodeFeatureIds[k];
					int lbl = nodeFeatureLabels[k];
					grad[fNum] += nodeFeatureGenerator.nFVal * tmpMarginal[lbl];
				}
			}
			for(int j = snGraph.graph.getNumNeighbours(node) - 1; j >= 0; --j){
				int nbr = snGraph.graph.getNeighbour(node, j);
				if(node < nbr) {
					double mu = 0;
					if(nodeFeatureGenerator != null)
						mu = lmbdInference.getLambdaEdgeSameLabelProbability(tmpModel, node, nbr);
					else
						mu = getEqualLabelEdgeProbability(tmpModel.getEdgePotential(node, 0, nbr, 0));
					for(int fNum = 0; fNum < numEdgeFeatures; ++fNum)
						grad[numNodeFeatures + fNum] += mu * edgeFeatureGenerator.getFeatureValue(node, nbr, fNum);
				}
			}
		}
		// Regularization stuff
		dataLL += -C_VALUE * ArrayUtils.norm2Square(lambda);
		for(int k=0; k<grad.length; ++k){ 
			grad[k] += 2 * C_VALUE * lambda[k];
		}
		return -dataLL;
	}
	
	protected void resetPotentials(UDGM model, double[] weights){
		NodePotentialSparse nonActiveNodePot = new NodePotentialSparse(0);
		for(int node = snGraph.graph.getNumNodes()-1; node >= 0; --node){
			if(nodeFeatureGenerator != null){
				NodePotentialSparse nPotS = new NodePotentialSparse(0);
				nodeFeatureGenerator.getFeatureValues(node, nodeFeatureIds, nodeFeatureLabels);
				for (int k = 0; k < nodeFeatureIds.length && nodeFeatureIds[k]!= -1; k++) {
					int fNum = nodeFeatureIds[k];
					int lbl = nodeFeatureLabels[k];
					nPotS.addPotential(lbl, weights[fNum] * nodeFeatureGenerator.nFVal);
				}
				model.setNodePotentialTable(node, nPotS);
			} else
				model.setNodePotentialTable(node, nonActiveNodePot);
			
			for(int j=snGraph.graph.getNumNeighbours(node)-1; j>=0; --j){
				int nbr = snGraph.graph.getNeighbour(node, j);
				if(node < nbr){
					double potts = 0;
					for(int fNum = 0; fNum < this.numEdgeFeatures; ++ fNum)
						potts += weights[numNodeFeatures + fNum] * edgeFeatureGenerator.getFeatureValue(node, nbr, fNum);
					model.setEdgePotentialTable(node, nbr, new PottsPotential(potts, 0));
				}
			}
		}	
	}
	
	public double getEqualLabelEdgeProbability(double edgePotts){
		assert(nodeFeatureGenerator == null);
		return Math.exp(edgePotts) / ( Math.exp(edgePotts) + snGraph.numLabels - 1 );
	}
}
