package iitb.sgl.learning;

public interface Learner {	
	public void learnWeights() throws Exception;
	public double[] getWeights();
}
