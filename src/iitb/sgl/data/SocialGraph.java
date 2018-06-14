package iitb.sgl.data;

import java.util.ArrayList;
import java.util.BitSet;

import gnu.trove.list.array.TIntArrayList;
import iitb.shared.graphs.UDGraph;

/**
 * A Base class representing directed Social graph.
 */

public class SocialGraph {
	
	public UDGraph graph;
	public TIntArrayList nodeLabels;
	public int numLabels;
	public BitSet observedNodes;
	public int numObservedNodes;
	public ArrayList<TIntArrayList> inEdgesMap;		// followers
	public ArrayList<TIntArrayList> outEdgesMap;	// friends
	public ArrayList<TIntArrayList> inOutEdgesMap;	// both followers and friends
	
	public boolean isNodeObserved(int node){
		return observedNodes.get(node);
	}
	
	public int totalDegree(int node){
		return graph.getNumNeighbours(node);
	}
	
	public int inDegree(int node) {
		return inEdgesMap.get(node).size() + inOutEdgesMap.get(node).size();
	}
	
	public int outDegree(int node) {
		return outEdgesMap.get(node).size() + inOutEdgesMap.get(node).size();
	}
	
	public int observedDegree(int node){
		int obsDegree = 0;
		for(int j = graph.getNumNeighbours(node)-1; j >= 0; --j){
			int nbr = graph.getNeighbour(node, j);
			if(observedNodes.get(nbr)) 
				++obsDegree;
		}
		return obsDegree;
	}
}

