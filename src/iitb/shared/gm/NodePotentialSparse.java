package iitb.shared.gm;

import java.util.Arrays;

import sun.reflect.generics.reflectiveObjects.NotImplementedException;

/**
 * This class represents sparse node potentials. 
 * Only active labels and their potentials are to be stored.
 * For rest of the labels, defPot is assumed.
 */
@SuppressWarnings("serial")
public class NodePotentialSparse extends NodePotential {
	
	private short labels[];	// Gaurish: Changed data type from byte(1-byte) to short(2 bytes)
	private double defPot;
	
	public NodePotentialSparse(){
		super(0);
	}
	public NodePotentialSparse(double defValue){
		this();
		defPot = defValue;
	}
	public NodePotentialSparse(NodePotentialSparse nps){
		super(nps);
		this.defPot = nps.defPot;
		short[] lbl = nps.getLabelsArray();
		labels = null;
		if(lbl != null){
			labels = new short[lbl.length];
			System.arraycopy(lbl, 0, this.labels, 0, lbl.length);
		}
	}
	/** Calls generic method to get clique potential value.*/
	public double getPotentialValue(int label){
		if (labels == null) return defPot;
		for (int i = 0; i < labels.length; i++) {
			if (labels[i]==label) return potentialValues[i];
		}
		return defPot;
	}
	public int getCliqueSize() {
		return 1;
	}
	public double getPotential(int[] labels) {
		return getPotentialValue(labels[0]);
	}
	public void setPotential(int[] labels, double potentialValue) {
		setPotential(labels[0], potentialValue);
	}
	public void setPotential(int label, double potentialValue) {
		if (labels != null) {
			for (int i = 0; i < labels.length; i++) {
				if (labels[i]==label) {
					potentialValues[i] = potentialValue;
					return;
				}
			}
		}
		short oldLabels[] = labels;
		double oldPots[] = potentialValues;
		labels = new short[labels==null?1:labels.length+1];
		potentialValues = new double[labels.length];
		if (oldLabels != null) {
			for (int i = 0; i < oldPots.length; i++) {
				labels[i] = oldLabels[i];
				potentialValues[i] = oldPots[i];
			}
		}
		labels[labels.length-1] = (short) label;
		potentialValues[potentialValues.length-1] = potentialValue;
		oldLabels = null;
		oldPots = null;
	}
	
	/**
	 * Similar to setPotential except that we have to add potVal 
	 * instead of replacing it for lbl.
	 */
	public void addPotential(int label, double potVal) {
		if (labels != null) {
			for (int i = 0; i < labels.length; i++) {
				if (labels[i]==label) {
					potentialValues[i] += potVal;
					return;
				}
			}
		}
		short oldLabels[] = labels;
		double oldPots[] = potentialValues;
		labels = new short[labels==null?1:labels.length+1];
		potentialValues = new double[labels.length];
		if (oldLabels != null) {
			for (int i = 0; i < oldPots.length; i++) {
				labels[i] = oldLabels[i];
				potentialValues[i] = oldPots[i];
			}
		}
		labels[labels.length-1] = (short) label;
		potentialValues[potentialValues.length-1] = potVal;
		oldLabels = null;
		oldPots = null;
	}
	
	public void zeroAllPotentials() {
		Arrays.fill(potentialValues, 0);
		defPot=0;
	}
	public String toString(){
		String str="";
		if (labels != null) {
			for (int i = 0; i < potentialValues.length; i++) {
				str += (labels[i] + ":" + potentialValues[i]) + " ";
			}
		}
		return str;
	}
	
	public double[] getPotentialsArray() {
		return potentialValues;
	}
	
	public int getNodeArity() {
		if(labels==null)
			return 0;
		return labels.length;
	}
	
	public short[] getLabelsArray() {
		return labels;
	}
	
	public short getLabel(int i) {
		return labels[i];
	}
	
	public double getDefaultPotentialValue(){
		return defPot;
	}
	
	public NodePotentialSparse clone() {
		return new NodePotentialSparse(this);
	}
	
	public void normalize() {
		throw new NotImplementedException();
	}
	
	public double getPotential(int[] vars, int[] fullLabels) {
		return getPotentialValue(fullLabels[vars[0]]);
	}
}