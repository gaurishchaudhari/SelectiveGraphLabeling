package iitb.shared.graphs;

import java.util.ArrayList;
import gnu.trove.list.array.TIntArrayList;
import iitb.shared.graphs.NbrIterImpl;

@SuppressWarnings("serial")
public class UDGraphAdjList implements UDGraph{

    private int numNodes;
    private int numEdges;
    ArrayList<TIntArrayList> adjList;

    public UDGraphAdjList(int numNodes) {
        this.numNodes = 0;
        this.numEdges = 0;
        adjList = new ArrayList<TIntArrayList>();
        for (int i=0; i < numNodes; i++) {
            addNode();
        }
    }

    public void addNode() {
        numNodes++;
        adjList.add(new TIntArrayList());
    }

    public UDGraphAdjList(UDGraph graph){
        this(graph.getNumNodes());
        for(int i=0; i<numNodes; i++){
            for(int j=0; j<numNodes; j++){
                if(graph.isAdj(i, j))
                    this.addEdge(i, j);
            }
        }
    }

    public int getNumNodes() {
        return numNodes;
    }

    public int getNumEdges() {
        return numEdges;
    }

    public boolean isAdj(int u, int v) {
        return adjList.get(u).contains(v);
    }

    public void addEdge(int u, int v) {
    	if(adjList.get(u).contains(v))
    		return;
        adjList.get(u).add(v);
        adjList.get(v).add(u);
        numEdges++;
    }

    public TIntArrayList getNeighbours(int u) {
        return adjList.get(u);
    }
    
    public int getNeighbour(int node, int nbrNum) {
        return adjList.get(node).get(nbrNum);
    }

    public int getNumNeighbours(int node) {
        return adjList.get(node).size(); 
    }

    public String toString(){
        StringBuffer sb = new StringBuffer(); 
        sb.append("(UnDirected) nV = " + numNodes + ", nE = " + numEdges + "\n");
        for (int i = 0; i < numNodes; i++) {
            sb.append(i + " :");
            TIntArrayList ni = getNeighbours(i);
            for (int j=0; j<ni.size(); j++)
                sb.append(" " + ni.get(j));
            sb.append("\n");
        }
        return sb.toString();
    }

	@Override
	public NbrIterator getNeighbourIterator(int node, NbrIterator iter) {
		if (iter != null)
			return ((NbrIterImpl)iter).init(node,this);
		return new NbrIterImpl(node,this);
	}
}
