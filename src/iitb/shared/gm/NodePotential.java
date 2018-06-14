package iitb.shared.gm;

import java.util.Arrays;

import iitb.shared.RobustMath;

/**
 * This class represents node potentials stored in an array.
 */
@SuppressWarnings("serial")
public class NodePotential implements PotentialInterface {
	
	protected double[] potentialValues;
	
	public NodePotential(){
	}
	
	public NodePotential(int nodeArity){
        potentialValues = new double[nodeArity];
	}
	
	public NodePotential(int nodeArity, double initValue){
        potentialValues = new double[nodeArity];
        Arrays.fill(potentialValues,initValue);
    }
	
	public NodePotential(NodePotential np){
		double[] from = np.getPotentialTable();
		potentialValues = new double[from.length];
		System.arraycopy(from, 0, potentialValues, 0, from.length);
	}
	
	private double[] getPotentialTable() {
		return potentialValues;
	}
    
	/** Calls generic method to get clique potential value.*/
	public double getPotentialValue(int label){
        return potentialValues[label];
	}
	
    public int getCliqueSize() {
        return 1;
    }
    public int getMaxNodeArity() {
        return potentialValues.length;
    }
    
    public double getPotential(int[] labels) {
        return potentialValues[labels[0]];
    }
    
    public void setPotential(int[] labels, double potentialValue) {
        potentialValues[labels[0]] = potentialValue;
    }
    
    public void setPotential(int label, double potentialValue) {
        potentialValues[label] = potentialValue;
    }
    
    public void addPotential(int label, double potVal) {
        potentialValues[label] += potVal;
    }
    
    public void zeroAllPotentials() {
        Arrays.fill(potentialValues, 0);
    }
    
    public String toString(){
       String str="";
       for (int i = 0; i < potentialValues.length; i++) {
           str += potentialValues[i] + " ";
       }
       return str;
    }
    
    public double[] getPotentialsArray() {
        return potentialValues;
    }
    
    public NodePotential clone() {
        return new NodePotential(this);
    }
    
	public void normalize() {
		double norm = RobustMath.LOG0;
		for (int i = 0; i < potentialValues.length; i++) {
			norm = RobustMath.logSumExp(norm, potentialValues[i]);
		}
		for (int i = 0; i < potentialValues.length; i++) {
			potentialValues[i] = potentialValues[i]-norm;
		}
	}
	
	public double getPotential(int[] vars, int[] fullLabels) {
		return getPotentialValue(fullLabels[vars[0]]);
	}
}
