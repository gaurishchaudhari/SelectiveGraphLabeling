package iitb.shared.gm.inference;

public class SolutionWithBounds extends Solution {
    public double bound;
    public SolutionWithBounds() {
        super();
    }

    public SolutionWithBounds(int n) {
        super(n);
    }

    public SolutionWithBounds(Solution s, double bound) {
        super(s);
        this.bound = bound;
        next = s.next;
    }
    
    public String toString(){
    	return "Bound: "+bound + " " + super.toString();
    }
    public SolutionWithBounds(int[] labelings, double score) {
        this(new Solution(labelings, score),Double.MAX_VALUE);
    }

}
