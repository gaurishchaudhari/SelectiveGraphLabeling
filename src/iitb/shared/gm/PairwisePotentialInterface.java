package iitb.shared.gm;

public interface PairwisePotentialInterface extends PotentialInterface {
	
    public double getPotentialValue(int label1, int label2);
    
    public void setPotentialValue(int label1, int label2, double val);
    
}
