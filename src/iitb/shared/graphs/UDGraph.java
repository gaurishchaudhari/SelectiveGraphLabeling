package iitb.shared.graphs;

import java.io.Serializable;

import gnu.trove.list.array.TIntArrayList;

public interface UDGraph extends Serializable {
	
	public  interface NbrIterator {
		int nextNbr();
	}
	
	public int getNumNodes();
	
	public int getNumEdges();			
	
	public boolean isAdj(int node1, int node2);
	
	public void addEdge(int node1, int node2);		// add only if doesn't exist
	
	public TIntArrayList getNeighbours(int node);
	
	public int getNumNeighbours(int node);
    
	public int getNeighbour(int node, int nbrNum);
    
    public NbrIterator getNeighbourIterator(int node, NbrIterator iter);
	
    public String toString();
    
}
