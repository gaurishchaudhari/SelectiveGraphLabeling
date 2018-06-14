package iitb.sgl.learning;

import java.util.BitSet;

import gnu.trove.iterator.TIntIntIterator;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.map.hash.TIntIntHashMap;
import iitb.sgl.data.EdgeFeatureGenerator;
import iitb.sgl.data.NodeFeatureGenerator;
import iitb.sgl.data.SocialGraph;
import iitb.shared.gm.NodePotentialSparse;
import iitb.shared.gm.PottsPotential;
import iitb.shared.gm.UDGM;
import iitb.shared.gm.inference.LabelScoreArraySparse;
import iitb.shared.graphs.UDGraph;
import iitb.shared.graphs.UDGraphAdjList;
import iitb.sgl.inference.SGInference;

/**
 * Mini-graph Utilities
*/

public class MiniGraphModelFactory {
	
	public class MiniGM{
		public TIntIntHashMap nodeIndex;
		public int[] reverseNodeIndex;
		public UDGraph sGraph;
		public UDGM gModel;
		public BitSet observedNodes;
	}
	
	public MiniGraphModelFactory(){
	}
	
	public static int getNumActiveNodes(MiniGM miniGM){
		return miniGM.observedNodes.cardinality();
	}
	
	public MiniGM createMiniGraph(SocialGraph snGraph, int startNode, int level){
		MiniGM miniGM = new MiniGM();
		int idx=0;
		TIntIntHashMap nodeIndex = new TIntIntHashMap();
		nodeIndex.put(startNode, idx++);
		TIntArrayList current = new TIntArrayList();
		TIntArrayList neighbours = new TIntArrayList();
		current.add(startNode);
		for(int lvl = 0; lvl < level;){
			if(!current.isEmpty()) {
				int node = current.removeAt(0);
				// Do not expand activeNodes except for startNode
				if(node != startNode && snGraph.isNodeObserved(node)) 
					continue;
				for(int j=snGraph.graph.getNumNeighbours(node)-1; j>=0; --j){
					int nbr = snGraph.graph.getNeighbour(node, j);
					if(!nodeIndex.containsKey(nbr)){	// Ensured nbr!= startNode 
						//Ignore leaf nodes which are non-active
						if(snGraph.graph.getNumNeighbours(nbr)>1 || snGraph.isNodeObserved(nbr) ){
							neighbours.add(nbr);
							nodeIndex.put(nbr, idx++);
						}
					}
				}
			}
			else{
				// TODO: Verify it !!
				if(neighbours.size() == 0)
					break;
				current = neighbours;
				neighbours = new TIntArrayList();
				++lvl;
			}
			if(idx >= snGraph.graph.getNumNodes())
				break;
		}
		miniGM.reverseNodeIndex = new int[snGraph.graph.getNumNodes()];
		
		/* Number the nodes for tree schedule of messages */
		TIntIntHashMap nodeIndex2 = new TIntIntHashMap(nodeIndex.size());
		TIntIntIterator iterTemp = nodeIndex.iterator();
		while(iterTemp.hasNext()){
			iterTemp.advance();
			int user = iterTemp.key();
			int val = idx - 1 - iterTemp.value();
			nodeIndex2.put(user, val);
			miniGM.reverseNodeIndex[val] = user;
		}
		nodeIndex = nodeIndex2;
		
		miniGM.nodeIndex = nodeIndex;
		miniGM.sGraph = new UDGraphAdjList(nodeIndex.size());
		miniGM.observedNodes = new BitSet(nodeIndex.size());
		TIntIntIterator iter = nodeIndex.iterator();
		while(iter.hasNext()){
			iter.advance();
			int user = iter.key();
			int userIdx = iter.value();
			if(snGraph.isNodeObserved(user)){
				miniGM.observedNodes.set(userIdx);
			}
			for(int j=snGraph.graph.getNumNeighbours(user)-1; j>=0; --j){
				int nbr = snGraph.graph.getNeighbour(user, j);
				if(nodeIndex.containsKey(nbr)) 
					miniGM.sGraph.addEdge(userIdx, nodeIndex.get(nbr));
			}
		}
		miniGM.observedNodes.clear(nodeIndex.get(startNode));
		miniGM.gModel = new UDGM(miniGM.sGraph, snGraph.numLabels, new PottsPotential(0.5, 0), 2);
		
		return miniGM;
	}
	
	public static void setMiniGraphPotentials(int startNode, SocialGraph snGraph, NodeFeatureGenerator nodeFeatureGenerator, EdgeFeatureGenerator edgeFeatureGenerator, MiniGM miniGM, double[] weights, SGInference solution){
		setMiniGraphEdgePotentials(snGraph, nodeFeatureGenerator, edgeFeatureGenerator, miniGM, weights);
		setMiniGraphNodePotentials(snGraph, nodeFeatureGenerator, miniGM, startNode, weights, solution);
	}
	
	public static void setMiniGraphEdgePotentials(SocialGraph snGraph, NodeFeatureGenerator nodeFeatureGenerator, EdgeFeatureGenerator edgeFeatureGenerator, MiniGM miniGM, double[] weights){
		int numNodeFeatures = 0;
		if(nodeFeatureGenerator != null) numNodeFeatures = nodeFeatureGenerator.numFeatures();
		int numEdgeFeatures = edgeFeatureGenerator.numFeatures();
		
		TIntIntIterator iter = miniGM.nodeIndex.iterator();
		while(iter.hasNext()){
			iter.advance();
			int user = iter.key();
			int userIdx = iter.value();
			for(int j=miniGM.sGraph.getNumNeighbours(userIdx)-1; j>=0; --j){
				int nbrIdx = miniGM.sGraph.getNeighbour(userIdx, j);
				if(userIdx < nbrIdx){
					double potts = 0;
					for(int fNum = 0; fNum < numEdgeFeatures; ++fNum)
						potts += weights[numNodeFeatures + fNum] * edgeFeatureGenerator.getFeatureValue(user, miniGM.reverseNodeIndex[nbrIdx], fNum);
					miniGM.gModel.setEdgePotentialTable(userIdx, nbrIdx, new PottsPotential(potts, 0));
				}
			}
		} 
	}
	
	public static void setMiniGraphNodePotentials(SocialGraph snGraph, NodeFeatureGenerator nodeFeatureGenerator, MiniGM miniGM, int startNode, double[] weights, SGInference solution){
		int numNodeFeatures = 0;
		int[] nodeFeatureIds = null;
		int[] nodeFeatureLabels = null;
		if(nodeFeatureGenerator != null){
			numNodeFeatures = nodeFeatureGenerator.numFeatures();
			nodeFeatureIds = new int[numNodeFeatures];
			nodeFeatureLabels = new int[nodeFeatureIds.length];
		}
		
		TIntIntIterator iter = miniGM.nodeIndex.iterator();
		NodePotentialSparse nonActiveNodePotential = new NodePotentialSparse(0);
		while(iter.hasNext()){
			iter.advance();
			int user = iter.key();
			int userIdx = iter.value();
			
			if(miniGM.observedNodes.get(userIdx) && user != startNode){	// user!=startNode is redundant (?)
				NodePotentialSparse activeNodePot = new NodePotentialSparse(0);
				activeNodePot.setPotential(snGraph.nodeLabels.get(user), 30);
				miniGM.gModel.setNodePotentialTable(userIdx, activeNodePot);
			}
			else{
				if(solution != null){
					System.out.println("Getting beliefs from outer graph");
					LabelScoreArraySparse partialBelief = new LabelScoreArraySparse(snGraph.numLabels);
					TIntArrayList minusNbrs = new TIntArrayList();
					for(int j=snGraph.graph.getNumNeighbours(user)-1; j>=0; --j){
						int nbr = snGraph.graph.getNeighbour(user, j);
						if(miniGM.nodeIndex.containsKey(nbr))	// startNode is also in minusNbrs
							minusNbrs.add(nbr);
					}
					double defVal = solution.getPartialBelief(user, minusNbrs, partialBelief);
					NodePotentialSparse partialActiveNodePot = new NodePotentialSparse(defVal - defVal);
					for(int k=partialBelief.getNumActive()-1; k>=0; --k){
						int lbl = partialBelief.getActiveLabel(k);
						partialActiveNodePot.setPotential(lbl, partialBelief.get(lbl) - defVal);
					}
					if(nodeFeatureGenerator != null){
						nodeFeatureGenerator.getFeatureValues(user, nodeFeatureIds, nodeFeatureLabels);
						for (int k = 0; k < nodeFeatureIds.length && nodeFeatureIds[k] != -1; k++) {
							int fNum = nodeFeatureIds[k];
							int lbl = nodeFeatureLabels[k];
							partialActiveNodePot.addPotential(lbl, nodeFeatureGenerator.nFVal * weights[fNum]);
						}
					}
					miniGM.gModel.setNodePotentialTable(userIdx, partialActiveNodePot);
				} else {
					if(nodeFeatureGenerator != null && weights != null){
						NodePotentialSparse nPotS = new NodePotentialSparse(0);
						nodeFeatureGenerator.getFeatureValues(user, nodeFeatureIds, nodeFeatureLabels);
						for (int k = 0; k < nodeFeatureIds.length && nodeFeatureIds[k] != -1; k++) {
							int fNum = nodeFeatureIds[k];
							int lbl = nodeFeatureLabels[k];
							nPotS.addPotential(lbl, weights[fNum] * nodeFeatureGenerator.nFVal);
						}
						miniGM.gModel.setNodePotentialTable(userIdx, nPotS);
					} else
						miniGM.gModel.setNodePotentialTable(userIdx, nonActiveNodePotential);
				}
			}
		}
	}
	
}
