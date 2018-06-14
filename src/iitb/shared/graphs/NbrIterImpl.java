package iitb.shared.graphs;

import java.io.Serializable;

import iitb.shared.graphs.UDGraph.NbrIterator;

@SuppressWarnings("serial")
public class NbrIterImpl implements NbrIterator, Serializable {
	
	int nbrNum=-1;
	int node;
	UDGraph graph;
	
	public NbrIterImpl(int node2, UDGraph graph) {
		init(node2,graph);
	}
	
	public NbrIterImpl init(int node2, UDGraph graph) {
		node = node2;
		this.graph = graph;
		nbrNum=-1;
		return this;
	}
	
	@Override
	public int nextNbr() {
		nbrNum++;
		if (nbrNum < graph.getNumNeighbours(node))
			return graph.getNeighbour(node, nbrNum);
		return -1;
	}
}
