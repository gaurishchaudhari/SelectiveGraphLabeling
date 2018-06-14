package iitb.shared.gm;

import sun.reflect.generics.reflectiveObjects.NotImplementedException;

@SuppressWarnings("serial")
public class PottsPotential implements PairwisePotentialInterface {
    double vSame;
    double vDiff;
    /*
     * @same is the potentials when the two labels are equal and 
     * @diff is the potential for when the labels are not equal. 
     */
    public PottsPotential(double same, double diff) {
        vSame = same;
        vDiff = diff;
    }
    public double getPotentialValue(int label1, int label2) {
        return (label1==label2)?vSame:vDiff;
    }
    public void setSamePotential(double same) {
    	vSame = same;
    }
    public void setDiffPotential(double val) {
    	vDiff = val;
    }
    public void setPotentialValue(int label1, int label2, double val) {
    	throw new NotImplementedException();
    }

    public int getCliqueSize() {
        return 2;
    }

    public int getMaxNodeArity() {
        throw new NotImplementedException();
    }

    public double getPotential(int[] labels) {
        return getPotentialValue(labels[0], labels[1]);
    }

    public void setPotential(int[] labels, double potentialValue) {
        setPotentialValue(labels[0],labels[1],potentialValue);
    }

    public void zeroAllPotentials() {
        vSame = vDiff = 0;
    }
    public PottsPotential clone() {
        return new PottsPotential(vSame,vDiff);
    }
    public PotentialInterface transpose() {
        return this;
    }
    
    public double getPotential(int[] vars, int[] fullLabels) {
		return getPotentialValue(fullLabels[vars[0]],fullLabels[vars[1]]);
	}

	public String toString() {return vSame + ":"+vDiff;}
}
