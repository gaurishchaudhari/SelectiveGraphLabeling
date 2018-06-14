package iitb.shared.gm;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Arrays;

import iitb.shared.ArrayUtils;
import iitb.shared.gm.PotentialInterface.PotType;
import iitb.shared.graphs.UDGraph;
import sun.reflect.generics.reflectiveObjects.NotImplementedException;

public class UDGM implements GraphicalModel{
	
	public String name;

	protected UDGraph graph;

	protected NodePotential nodePotentials[];

	protected EdgePotentialStore edgePotentials;

	protected int[] nodeArities;

	protected PotentialInterface defaultCliquePotential;

	protected int maxNodeArity;

	protected int maxCliqueSize;

	protected int numOfNodes;
	
	private double constantPot;
	
	PotType potType = PotType.DENSE;
	
	public UDGM(UDGraph graph, int uniformNodeArity, PotentialInterface defCliquePotential, int cliqueSize){
		init(graph,defCliquePotential,graph.getNumNodes(),uniformNodeArity,cliqueSize);
	}
	
	void init(UDGraph graph, PotentialInterface defCliquePotential, int numOfNodes, int uniformNodeArity, int maxCliqueSize){
		nodeArities = new int[numOfNodes];
		Arrays.fill(nodeArities, uniformNodeArity);
		assert(graph == null || graph.getNumNodes()==numOfNodes);
		init(graph,nodeArities,defCliquePotential,maxCliqueSize);
	}
	
	public void init(UDGraph graph, int[] nodeArities, PotentialInterface defaultPotential, int maxCliqueSize) {
		this.graph = graph;
		this.numOfNodes = nodeArities.length;
		this.maxCliqueSize = maxCliqueSize;
		defaultCliquePotential = defaultPotential;
		this.nodeArities = nodeArities;
		this.maxNodeArity = ArrayUtils.max(nodeArities);
		
		nodePotentials = new NodePotential[numOfNodes];
		edgePotentials = new EdgePotentialStore(); 
	}

	@Override
	public UDGraph getGraph() {
		return graph;
	}

	@Override
	public int getMaxArity() {
		 return maxNodeArity;
	 }

	@Override
	public int getNodeArity(int num) {
		return nodeArities[num];
	}

	@Override
	public int getMaxCliqueSize() {
		return maxCliqueSize;
	}
	
	@Override
	public double getConstantPotential() {
		return constantPot;
	}

	@Override
	public double getNodePotential(int node, int label){
		 if(nodePotentials[node] != null)
			 return nodePotentials[node].getPotentialValue(label);
		 else
			 return 0;
	 }

	@Override
	public double getEdgePotential(int node1, int label1, int node2, int label2){
		if(!graph.isAdj(node1, node2))
			return 0;
		
		int[] nodes = new int[2];
		int[] labels = new int[2];
		nodes[0] = (node1 < node2) ? node1 : node2;
		nodes[1] = (node1 < node2) ? node2 : node1;
		labels[0] = (node1 < node2) ? label1 : label2;
		labels[1] = (node1 < node2) ? label2 : label1;
		
		PotentialInterface pot =  edgePotentials.getCliquePotentialTable(nodes[0], nodes[1]);
		return (pot != null) ? pot.getPotential(labels) : 0;
	}
	
	public void setNodePotentialTable(int node, PotentialInterface potential) {
		nodePotentials[node] = (NodePotential) potential;
	}
	
	public void setEdgePotentialTable(int j, int i, PairwisePotentialInterface edgePot) {
		if (j > i) {
			System.err.println("Potential tables need to be in increasing order of nodeIds");
			throw new NotImplementedException();
		}
		edgePotentials.setCliquePotentialTable(j, i, edgePot);
	}
	
	@Override
	public NodePotential getNodePotentialTable(int node) {
		return nodePotentials[node];
	}

	@Override
	public PairwisePotentialInterface getEdgePotentialTable(int node1, int node2){
		if(!graph.isAdj(node1, node2))
			return null;
		
		int[] nodes = new int[2];
		nodes[0] = (node1 < node2) ? node1 : node2;
		nodes[1] = (node1 < node2) ? node2 : node1;
		
		return (PairwisePotentialInterface)edgePotentials.getCliquePotentialTable(nodes[0], nodes[1]);
	}
	
	@Override
	public double getScore(int[] labeling) {
		double score = getConstantPotential();
		for(int i=0;i < numOfNodes;i++) {           
			if (labeling[i] >= 0) score += getNodePotential(i, labeling[i]);
		}
		//System.out.println("Sum of node potentials "+score);
		for(int i=0;i < numOfNodes;i++) {
			int sz = graph.getNumNeighbours(i);
			for(int j=0; j<sz; j++){
				int nbr = graph.getNeighbour(i,j);
				if(nbr > i) {
					score += getEdgePotential(i, labeling[i], nbr, labeling[nbr]);
				}
			}
		}
		
		return score;
	}
	
	@SuppressWarnings("serial")
	class EdgePotentialStore implements Serializable {
		ArrayList<PairwisePotentialInterface>[] edgePotentials;
		
		public PairwisePotentialInterface getCliquePotentialTable(int n0, int n1) {
			if (!graph.isAdj(n0, n1)) return null;
			if (edgePotentials==null) allocate(); 
			int j = graph.getNeighbours(n0).indexOf(n1);
			return edgePotentials[n0].get(j);
		}
		
		public void setCliquePotentialTable(int n0, int n1, PotentialInterface potential) {
			assert (graph.isAdj(n0, n1));
			if (edgePotentials==null) allocate();
			int j = graph.getNeighbours(n0).indexOf(n1);
			edgePotentials[n0].set(j,(PairwisePotentialInterface) potential);
		}
		
		@SuppressWarnings("unchecked")
		private void allocate() {
			edgePotentials = (ArrayList<PairwisePotentialInterface>[]) new ArrayList[graph.getNumNodes()];
			for (int i = 0; i < nodeArities.length; ++i) {
				edgePotentials[i] = new ArrayList<PairwisePotentialInterface>(graph.getNumNeighbours(i));
				for(int j = edgePotentials[i].size() - 1; j >= 0; --j) {
					edgePotentials[i].set(j, null);
				}
			}
		}
	}
	
}
