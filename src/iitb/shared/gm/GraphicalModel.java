package iitb.shared.gm;

import iitb.shared.graphs.UDGraph;


/**
 * This abstract class represents template for any graphical model representation.
 * All classes defining graphical model representation must inherit this class. 
 */
public interface GraphicalModel {
	
    @SuppressWarnings("serial")
	public class ParameterException extends Exception {
        public ParameterException(String string) {
            super(string);
        }
    }

	public int getMaxArity();
	
    public int getNodeArity(int num);

    abstract public int getMaxCliqueSize();
    
	public double getNodePotential(int node, int label);
    
    public double getEdgePotential(int node1, int label1, int node2, int label2);

    public NodePotential getNodePotentialTable(int i);
    
    public PairwisePotentialInterface getEdgePotentialTable(int n1, int n2);
    
    public UDGraph getGraph();
 	
    public double getScore(int labeling[]);

    public double getConstantPotential();

    abstract public String toString();
}
