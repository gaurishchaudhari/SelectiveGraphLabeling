package iitb.sgl.learning;

import java.util.Arrays;
import java.util.BitSet;
import java.util.HashMap;

import iitb.sgl.data.EdgeFeatureGenerator;
import iitb.sgl.data.NodeFeatureGenerator;
import iitb.shared.RobustMath;
import iitb.shared.gm.NodeMarginalSparse;
import iitb.shared.gm.NodePotentialSparse;
import iitb.shared.gm.UDGM;
import iitb.shared.gm.inference.LabelScoreArraySparse;
import iitb.shared.graphs.UDGraph;
import iitb.sgl.inference.SGInference;

/**
 * Tamir Hazan NIPS-09: Approximate logZ computation.
 * Utilities for LAMBDA messages.
 */

public class LambdaInference {
	UDGraph graph;
	int numLabels;
	NodeFeatureGenerator nodeFeatureGenerator;
	EdgeFeatureGenerator edgeFeatureGenerator;
	HashMap<Integer, NodeMarginalSparse>[] lambda;
	
	public LambdaInference(UDGraph graph, int numLabels, NodeFeatureGenerator nodeFeatureGenerator, EdgeFeatureGenerator edgeFeatureGenerator){
		this.graph = graph;
		this.numLabels = numLabels;
		this.nodeFeatureGenerator = nodeFeatureGenerator;
		this.edgeFeatureGenerator = edgeFeatureGenerator;
	}

	@SuppressWarnings("unchecked")
	public void storeLambdaValues(SGInference solution){
		LabelScoreArraySparse tempBuf = new LabelScoreArraySparse(numLabels);
		lambda = new HashMap[graph.getNumNodes()];
		for(int node=graph.getNumNodes()-1; node>=0; --node){
			lambda[node] = new HashMap<Integer, NodeMarginalSparse>(graph.getNumNeighbours(node));
			for(int j=graph.getNumNeighbours(node)-1; j>=0; --j){
				int nbr = graph.getNeighbour(node, j);
				float dV = (float)solution.getEdgeLambda(tempBuf, node, nbr);
				NodeMarginalSparse lmbd = new NodeMarginalSparse(tempBuf.getNumActive());
				lmbd.defVal = dV;
				for(int i=tempBuf.getNumActive()-1; i>=0; --i){
					int lbl = tempBuf.getActiveLabel(i);
					lmbd.activeLabels.set(lbl);
					lmbd.setScore(lbl, tempBuf.get(lbl));
				}
				lambda[node].put(nbr, lmbd);
			}
		}
	}
	
	public NodeMarginalSparse getEdgeLambda(int src, int dest){
		return lambda[src].get(dest);
	}
	
	NodeMarginalSparse lmbdSrc;
	NodeMarginalSparse lmbdDst;
	
	/**
	 * Convex Approximation of LogZ for general Graphs.
	 */
	public double getConvexifiedLogZ(UDGM model) throws Exception{
		int  arity = model.getMaxArity();
		double nodeTerm = 0;
		double edgeTerm = 0;
		double[] aggrScores = new double[arity];
		for(int node = model.getGraph().getNumNodes() - 1; node >= 0; --node){
			if(nodeFeatureGenerator != null){
				NodePotentialSparse nPotS = (NodePotentialSparse)model.getNodePotentialTable(node);
				Arrays.fill(aggrScores, 0);
				for(int j = model.getGraph().getNumNeighbours(node) - 1; j >= 0; --j){
					int nbr = model.getGraph().getNeighbour(node, j);
					lmbdSrc = getEdgeLambda(node, nbr);
					for(int k=0; k<arity; ++k)
						aggrScores[k] += lmbdSrc.getScore(k);
				}
				double term = RobustMath.LOG0;
				for(int k = arity - 1; k >= 0; --k)
					term = RobustMath.logSumExp(term, nPotS.getPotentialValue(k) - aggrScores[k]);
				nodeTerm += term;
			}
			for(int j = model.getGraph().getNumNeighbours(node) - 1; j >= 0; --j){
				int nbr = model.getGraph().getNeighbour(node, j);
				if(node > nbr) continue;
				lmbdSrc = getEdgeLambda(node, nbr); 
				lmbdDst = getEdgeLambda(nbr, node);
				
				double term = RobustMath.LOG0;
				int numCommonActiveLbl = 0;
				for(int x = lmbdSrc.activeLabels.nextSetBit(0); x >= 0; x = lmbdSrc.activeLabels.nextSetBit(x+1)){
					for(int y = lmbdDst.activeLabels.nextSetBit(0); y >= 0; y = lmbdDst.activeLabels.nextSetBit(y+1)){
						term = RobustMath.logSumExp(term, model.getEdgePotential(node, x, nbr, y) + lmbdSrc.getScore(x) + lmbdDst.getScore(y));
						if(x==y) ++numCommonActiveLbl;
					}
				}
				for(int x = lmbdSrc.activeLabels.nextSetBit(0); x >= 0; x = lmbdSrc.activeLabels.nextSetBit(x+1)){
					if(lmbdDst.activeLabels.get(x))
						term = RobustMath.logSumExp(term, Math.log(arity - lmbdDst.activeLabels.cardinality()) + lmbdSrc.getScore(x) + lmbdDst.defVal);
					else{
						term = RobustMath.logSumExp(term, Math.log(arity - lmbdDst.activeLabels.cardinality() - 1)  + lmbdSrc.getScore(x) + lmbdDst.defVal);
						term = RobustMath.logSumExp(term, model.getEdgePotential(node, 0, nbr, 0) + lmbdSrc.getScore(x) + lmbdDst.defVal);
					}
				}
				
				for(int y = lmbdDst.activeLabels.nextSetBit(0); y >= 0; y = lmbdDst.activeLabels.nextSetBit(y+1)){
					if(lmbdSrc.activeLabels.get(y))
						term = RobustMath.logSumExp(term, Math.log(arity - lmbdSrc.activeLabels.cardinality()) + lmbdSrc.defVal + lmbdDst.getScore(y));
					else{
						term = RobustMath.logSumExp(term, Math.log(arity - lmbdSrc.activeLabels.cardinality() - 1)  + lmbdSrc.defVal + lmbdDst.getScore(y));
						term = RobustMath.logSumExp(term, model.getEdgePotential(node, 0, nbr, 0) + lmbdSrc.defVal + lmbdDst.getScore(y));
					}
				}
				int remArity = arity + numCommonActiveLbl - lmbdSrc.activeLabels.cardinality() - lmbdDst.activeLabels.cardinality();
				term = RobustMath.logSumExp(term, Math.log(remArity) + (model.getEdgePotential(node, 0, nbr, 0) + lmbdSrc.defVal + lmbdDst.defVal) );
				term = RobustMath.logSumExp(term, Math.log((arity-lmbdSrc.activeLabels.cardinality())*(arity-lmbdDst.activeLabels.cardinality())-remArity) + lmbdSrc.defVal + lmbdDst.defVal );
				
				edgeTerm += term;
			}
		}
		if(nodeFeatureGenerator!= null)
			return (nodeTerm + edgeTerm);
		else
			return edgeTerm;
	}
	
	public void getLambdaNodeMarginal(UDGM model, int nodeId, double[] marginal, BitSet activeNodes){
		int arity = model.getMaxArity();
		NodePotentialSparse nPotS = (NodePotentialSparse) model.getNodePotentialTable(nodeId);
		Arrays.fill(marginal, nPotS.getDefaultPotentialValue());
		for(int k = nPotS.getNodeArity() - 1; k >= 0; --k){
			int lbl = nPotS.getLabel(k);
			marginal[lbl] = nPotS.getPotentialValue(lbl);
		}
		for(int j = model.getGraph().getNumNeighbours(nodeId) - 1; j >= 0; --j){
			int nbr = model.getGraph().getNeighbour(nodeId, j);
			lmbdSrc = getEdgeLambda(nodeId, nbr);
			for(int k = 0; k < arity; ++k){
				marginal[k] += -lmbdSrc.getScore(k);
			}
		}
		double normalizer = RobustMath.LOG0;
		for(int i = 0; i < marginal.length; ++i)
			normalizer = RobustMath.logSumExp(normalizer, marginal[i]);
		for(int i = 0; i < marginal.length; ++i)
			marginal[i] = Math.exp(marginal[i] - normalizer);
	}
	
	public double getLambdaEdgeSameLabelProbability(UDGM model, int src, int dest) throws Exception{
		int arity = model.getNodeArity(src);
		assert(arity==model.getNodeArity(dest));
		// Compute Incoming Messages to SRC and DEST
		lmbdSrc = getEdgeLambda(src, dest);
		lmbdDst = getEdgeLambda(dest, src); 

		// COMPUTATIONS
		double potts = model.getEdgePotential(src, 0, dest, 0);
		double logNr = RobustMath.LOG0;
		double logTerm1 = RobustMath.LOG0;
		double logTerm2 = RobustMath.LOG0;
		double logRepeated =  RobustMath.LOG0;
		
		BitSet activeLabels = new BitSet();
		activeLabels.or(lmbdSrc.activeLabels);
		activeLabels.or(lmbdDst.activeLabels);
		
		for(int x=activeLabels.nextSetBit(0); x>=0; x=activeLabels.nextSetBit(x+1)){
			logNr = RobustMath.logSumExp(logNr, potts + lmbdSrc.getScore(x) + lmbdDst.getScore(x));
			logRepeated = RobustMath.logSumExp(logRepeated, lmbdSrc.getScore(x) + lmbdDst.getScore(x));
		}
		logNr = RobustMath.logSumExp(logNr, Math.log(arity-activeLabels.cardinality()) + potts + lmbdSrc.defVal + lmbdDst.defVal);
		logRepeated = RobustMath.logSumExp(logRepeated, Math.log(arity-activeLabels.cardinality()) + lmbdSrc.defVal + lmbdDst.defVal);

		for(int x=lmbdSrc.activeLabels.nextSetBit(0); x>=0; x=lmbdSrc.activeLabels.nextSetBit(x+1))
			logTerm1 = RobustMath.logSumExp(logTerm1, lmbdSrc.getScore(x));
		logTerm1 = RobustMath.logSumExp(logTerm1, Math.log(arity-lmbdSrc.activeLabels.cardinality()) + lmbdSrc.defVal);
		for(int y=lmbdDst.activeLabels.nextSetBit(0); y>=0; y=lmbdDst.activeLabels.nextSetBit(y+1))
			logTerm2 = RobustMath.logSumExp(logTerm2, lmbdDst.getScore(y));
		logTerm2 = RobustMath.logSumExp(logTerm2, Math.log(arity-lmbdDst.activeLabels.cardinality()) + lmbdDst.defVal);
		
		double logTerm = logTerm1 + logTerm2;
		logTerm = RobustMath.logMinusExp(logTerm, logRepeated);// + Math.log(arity*(arity-1));
		double logDn = RobustMath.logSumExp(logNr, logTerm);
		return Math.exp(logNr-logDn);
	}
}