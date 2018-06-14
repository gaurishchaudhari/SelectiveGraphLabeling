package iitb.shared.gm.inference;

import java.util.Arrays;

public class Labeling {
    public int[] labeling;
    public double score;
    public Labeling(){}
    public Labeling(int n) {
        labeling = new int[n];
    }
    public Labeling(int[] labelings, double score) {
        this.labeling = labelings;
        this.score = score;
    }
    public String signature() {
        return Arrays.toString(labeling);
    }
    public String toString() {
        return "Score: " + score + " " + signature();
    }
}
