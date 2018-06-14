package iitb.shared.gm;

import java.io.Serializable;

public interface PotentialInterface extends Serializable {
    
    public static enum PotType {DENSE, SPARSE, POTTS};
    
    public double getPotential(int[] labels);
    
    public double getPotential(int vars[], int[] fullLabels);
  
    public void setPotential(int[] labels, double potentialValue);

    /** Returns clique size for which this potential is defined.*/
    public int getCliqueSize();
    
    /** Returns maximum node arity used in defining clique potential.*/
    public int getMaxNodeArity();
    
    public PotentialInterface clone();
    
    
}
