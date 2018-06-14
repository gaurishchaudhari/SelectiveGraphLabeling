package iitb.sgl.inference;

import iitb.sgl.data.EdgeFeatureGenerator;
import iitb.sgl.data.NodeFeatureGenerator;
import iitb.sgl.data.SocialGraph;
import iitb.shared.gm.NodePotentialSparse;
import iitb.shared.gm.PottsPotential;
import iitb.shared.gm.UDGM;

/**
 * Builds a Weighted / Uniform Potts-Potential Graphical Model. 
 * Sets node and edge potentials.
 */

public class SGModelUtils {
	
	public static void setPotentials(SocialGraph snGraph, UDGM model, double[] weights, NodeFeatureGenerator nodeFeatureGenerator, EdgeFeatureGenerator edgeFeatureGenerator) throws Exception{
		int numNodeFeatures = 0;
		int[] nodeFeatureIds = null;
		int[] nodeFeatureLabels = null;
		if(nodeFeatureGenerator != null){
			numNodeFeatures = nodeFeatureGenerator.numFeatures();
			nodeFeatureIds = new int[numNodeFeatures];
			nodeFeatureLabels = new int[nodeFeatureIds.length];
		}
		int numEdgeFeatures = edgeFeatureGenerator.numFeatures();
		double avgPotts = 0;
		int nearHalfCount = 0;
		int count = 0;
		NodePotentialSparse nonActiveNodePot = new NodePotentialSparse(0);
		for(int node = snGraph.graph.getNumNodes()-1; node >= 0; --node){
			if(snGraph.isNodeObserved(node)){
				NodePotentialSparse activeNodePot = new NodePotentialSparse(0);
				activeNodePot.setPotential(snGraph.nodeLabels.get(node), 30);
				model.setNodePotentialTable(node, activeNodePot);
			}
			else{
				if(nodeFeatureGenerator == null){
					model.setNodePotentialTable(node, nonActiveNodePot);
				} else{
					//TODO: Set nodePotentials using nodeFeatures
					NodePotentialSparse nPotS = new NodePotentialSparse(0);
					nodeFeatureGenerator.getFeatureValues(node, nodeFeatureIds, nodeFeatureLabels);
					for (int k = 0; k < nodeFeatureIds.length && nodeFeatureIds[k] >= 0; ++k) {
						int fNum = nodeFeatureIds[k];
						int lbl = nodeFeatureLabels[k];
						nPotS.addPotential(lbl, weights[fNum] * nodeFeatureGenerator.nFVal);
					}
					//System.out.println(nPotS.toString());
					model.setNodePotentialTable(node, nPotS);
				}
			}
			for(int j = snGraph.graph.getNumNeighbours(node)-1; j >= 0; --j){
				int nbr = snGraph.graph.getNeighbour(node, j);
				if(node < nbr){
					double potts = 0;
					for(int fNum = 0; fNum < numEdgeFeatures; ++fNum){
						assert(!Double.isInfinite(edgeFeatureGenerator.getFeatureValue(node, nbr, fNum))):""+node+" "+nbr+" "+fNum;
						potts += weights[numNodeFeatures + fNum] * edgeFeatureGenerator.getFeatureValue(node, nbr, fNum);
					}
					//if(potts < 0) potts = 10e-5; 
					avgPotts += potts;
					nearHalfCount += (Math.abs(potts - 0.5) < 0.1) ? 1 : 0;
					++count;
					model.setEdgePotentialTable(node, nbr, new PottsPotential(potts, 0));
				}
			}
		} 
		System.out.println("Average Potts = "+avgPotts+"/"+count+" = "+String.format("%.4f",(avgPotts*1.0/count)) + " Fraction of edges with ~0.5 potts = "+(nearHalfCount*1.0/count));
	}
	
	public static void setNodePotentials(SocialGraph snGraph, UDGM model, double[] weights, NodeFeatureGenerator nodeFeatureGenerator, EdgeFeatureGenerator edgeFeatureGenerator) throws Exception{
		int numNodeFeatures = 0;
		int[] nodeFeatureIds = null;
		int[] nodeFeatureLabels = null;
		if(nodeFeatureGenerator != null){
			numNodeFeatures = nodeFeatureGenerator.numFeatures();
			nodeFeatureIds = new int[numNodeFeatures];
			nodeFeatureLabels = new int[nodeFeatureIds.length];
		}
		NodePotentialSparse nonActiveNodePot = new NodePotentialSparse(0);
		for(int node = snGraph.graph.getNumNodes()-1; node >= 0; --node){
			if(snGraph.isNodeObserved(node)){
				NodePotentialSparse activeNodePot = new NodePotentialSparse(0);
				activeNodePot.setPotential(snGraph.nodeLabels.get(node), 30);
				model.setNodePotentialTable(node, activeNodePot);
			}
			else{
				if(nodeFeatureGenerator == null || weights == null){
					model.setNodePotentialTable(node, nonActiveNodePot);
				} else{
					//TODO: Set nodePotentials using nodeFeatures
					NodePotentialSparse nPotS = new NodePotentialSparse(0);
					nodeFeatureGenerator.getFeatureValues(node, nodeFeatureIds, nodeFeatureLabels);
					for (int k = 0; k < nodeFeatureIds.length && nodeFeatureIds[k] >= 0; ++k) {
						int fNum = nodeFeatureIds[k];
						int lbl = nodeFeatureLabels[k];	
						nPotS.addPotential(lbl, weights[fNum] * nodeFeatureGenerator.nFVal);
					}
					model.setNodePotentialTable(node, nPotS);
				}
			}
		} 
	}
}

