package iitb.sgl.learning;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Properties;
import java.util.Random;

import gnu.trove.iterator.TIntIntIterator;
import iitb.sgl.data.SocialGraph;
import iitb.sgl.inference.SGInference;
import iitb.sgl.inference.SGModelUtils;
import iitb.shared.ArrayUtils;
import iitb.shared.gm.NodeMarginalSparse;
import iitb.shared.gm.PottsPotential;
import iitb.shared.gm.UDGM;
import iitb.shared.gm.inference.ScalableMessagePassing.Method;

/**
 * Arbitrary Graph Learning
 * Node-level Conditional Likelihood (NCL) training
 */

public class NCL extends LearnerImpl {
	
	MiniGraphModelFactory miniGMFactory;
	SGInference mSolution;
	ArrayList<NodeMarginalSparse>[] nodeNbrMarginals;
	float[][] expectedFeatures;
	UDGM model;
	SGInference solution;
	int BFS_LEVEL = 2;	
	boolean isNodePotActive = false;
	
	public NCL(SocialGraph snGraph, Properties options, Random random, int debugLvl,
			boolean isNodePotActive, double[] initWeights, double C) throws Exception {
		super(snGraph, options, random, debugLvl, initWeights, C);
		this.isNodePotActive = isNodePotActive;
	}
	
	@Override
	public void learnWeights() throws Exception {
		System.out.println("Initial Weights: " + ArrayUtils.toStrArr(weights));
		model = new UDGM(snGraph.graph, snGraph.numLabels, new PottsPotential(0.5, 0), 2);
		miniGMFactory = new MiniGraphModelFactory();
		miniGMs  = new MiniGraphModelFactory.MiniGM[snGraph.numObservedNodes];
		for(int node=snGraph.observedNodes.nextSetBit(0), index=0; node>=0; node=snGraph.observedNodes.nextSetBit(node+1), ++index){
			MiniGraphModelFactory.MiniGM miniGM = miniGMFactory.createMiniGraph(snGraph, node, BFS_LEVEL);
			miniGMs[index] = miniGM;
		}
		
		expectedFeatures = new float[snGraph.numObservedNodes][weights.length];
		double oldNorm = 0, newNorm = 0;
		int iter = 0;
		MAX_ITERS = 100;
		do{
			SGModelUtils.setPotentials(snGraph, this.model, this.weights, nodeFeatureGenerator, edgeFeatureGenerator);
			solution = new SGInference(model, 1, Method.BP, false, snGraph.observedNodes);
			solution.setIsDegreeInvGamma(true);
			solution.getTopK(1);

			double[] marginal = new double[snGraph.numLabels];
			for(int user=snGraph.observedNodes.nextSetBit(0), index=0; user>=0; user=snGraph.observedNodes.nextSetBit(user+1), ++index){
				MiniGraphModelFactory.MiniGM miniGM = miniGMs[index];
				Arrays.fill(expectedFeatures[index], 0);
				TIntIntIterator itr = miniGM.nodeIndex.iterator();
				while(itr.hasNext()){
					itr.advance();
					int node = itr.key();
					int nodeIdx = itr.value();
					if(nodeFeatureGenerator != null && !snGraph.isNodeObserved(node)){
						solution.getNodeSumMarginal(node, marginal, null);
						nodeFeatureGenerator.getFeatureValues(node, nodeFeatureIds, nodeFeatureLabels);
						for (int k = 0; k < nodeFeatureIds.length && nodeFeatureIds[k] != -1; k++) {
							int fNum = nodeFeatureIds[k];
							int lbl = nodeFeatureLabels[k];
							expectedFeatures[index][fNum] += nodeFeatureGenerator.nFVal * marginal[lbl];
						}
					}
					for(int j=miniGM.sGraph.getNumNeighbours(nodeIdx)-1; j>=0; --j){
						int nbrIdx = miniGM.sGraph.getNeighbour(nodeIdx, j);
						if(nodeIdx < nbrIdx){	
							int nbr = miniGM.reverseNodeIndex[nbrIdx];
							double mu = solution.getEdgeSameLabelProbability(node, nbr);
							for(int fNum = 0; fNum < numEdgeFeatures; ++fNum)
								expectedFeatures[index][numNodeFeatures + fNum] += mu * edgeFeatureGenerator.getFeatureValue(node, nbr, fNum);
						}
					}		
				}
			}
			
			oldNorm = ArrayUtils.norm2(weights);
			maximizationStep();
			newNorm = ArrayUtils.norm2(weights);
			
			if(debugLvl > 0) System.out.println("E-Step "+iter+" Weights: "+ArrayUtils.toStrArr(weights)+" || "+newNorm);
			++iter;
		}while(iter < MAX_ITERS && Math.abs(newNorm - oldNorm) > GRAD_DELTA);
	}
	
	LambdaInference[] lmbdInference;
	MiniGraphModelFactory.MiniGM[] miniGMs;
	double[] miniMarginal;
	
	public void maximizationStep() throws Exception{
		lmbdInference = new LambdaInference[snGraph.numObservedNodes];
		miniMarginal = new double[snGraph.numLabels];
		double oldNorm = 0, newNorm = 0;
		int iter = 0;
		do{	
			for(int user=snGraph.observedNodes.nextSetBit(0), index = 0; user >= 0; user = snGraph.observedNodes.nextSetBit(user+1), ++index){
				MiniGraphModelFactory.MiniGM miniGM = miniGMs[index];
				MiniGraphModelFactory.setMiniGraphPotentials(user, snGraph, nodeFeatureGenerator, edgeFeatureGenerator, miniGM, weights, null);
				mSolution = new SGInference(miniGM.gModel, 0, Method.BP, false, miniGM.observedNodes);
				mSolution.setIsDegreeInvGamma(true);
				mSolution.setMaxIters(1);
				mSolution.getTopK(1);
				lmbdInference[index] = new LambdaInference(miniGM.sGraph, snGraph.numLabels, nodeFeatureGenerator, edgeFeatureGenerator);
				lmbdInference[index].storeLambdaValues(mSolution);
			}
			
			oldNorm = ArrayUtils.norm2(weights);
			trainer.optimize(weights, this);
			newNorm = ArrayUtils.norm2(weights);
			if(debugLvl>0) System.out.println("M-Step "+iter+" Weights: "+ArrayUtils.toStrArr(weights)+" || "+newNorm);
			++iter;		
		}while(iter < 100 && Math.abs(newNorm - oldNorm) > GRAD_DELTA);
	}
	
	@Override
	public double computeFunctionGradient(double[] lambda, double[] grad) throws Exception {
		double dataLL = 0;
		Arrays.fill(grad, 0);
		for(int user=snGraph.observedNodes.nextSetBit(0), index = 0; user >= 0; user = snGraph.observedNodes.nextSetBit(user+1), ++index){
			MiniGraphModelFactory.MiniGM miniGM = miniGMs[index];
			MiniGraphModelFactory.setMiniGraphPotentials(user, snGraph, nodeFeatureGenerator, edgeFeatureGenerator, miniGM, lambda, null);
			
			TIntIntIterator itr = miniGM.nodeIndex.iterator();
			while(itr.hasNext()){
				itr.advance();
				int node = itr.key();
				int nodeIdx = itr.value();
				if(nodeFeatureGenerator != null && !miniGM.observedNodes.get(nodeIdx)){
					lmbdInference[index].getLambdaNodeMarginal(miniGM.gModel, nodeIdx, miniMarginal, null);
					nodeFeatureGenerator.getFeatureValues(node, nodeFeatureIds, nodeFeatureLabels);
					for (int k = 0; k < nodeFeatureIds.length && nodeFeatureIds[k] >= 0; k++) {
						int fNum = nodeFeatureIds[k];
						int lbl = nodeFeatureLabels[k];
						grad[fNum] += nodeFeatureGenerator.nFVal * miniMarginal[lbl];
					}
				}
				for(int j=miniGM.sGraph.getNumNeighbours(nodeIdx)-1; j>=0; --j){
					int nbrIdx = miniGM.sGraph.getNeighbour(nodeIdx, j);
					if(nodeIdx < nbrIdx){
						double muMini = lmbdInference[index].getLambdaEdgeSameLabelProbability(miniGM.gModel, nodeIdx, nbrIdx);
						for(int fNum = 0; fNum < numEdgeFeatures; ++fNum)
							grad[numNodeFeatures + fNum] += muMini * edgeFeatureGenerator.getFeatureValue(node, miniGM.reverseNodeIndex[nbrIdx], fNum);
					}
				}
			}
			double logZ = lmbdInference[index].getConvexifiedLogZ(miniGM.gModel);
			dataLL += ArrayUtils.dotProduct(lambda, expectedFeatures[index]) - logZ ;
			for(int fNum = 0; fNum < grad.length; ++fNum)
				grad[fNum] += -expectedFeatures[index][fNum];
			
		}
		// Regularization stuff
		dataLL += -C_VALUE * ArrayUtils.norm2Square(lambda);
		for(int k=0; k<grad.length; ++k){ 
			grad[k] += 2 * C_VALUE * lambda[k];
		}
		return -dataLL;
	}	
}
