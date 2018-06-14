package iitb.shared.gm.inference;

import java.util.Arrays;

import iitb.shared.ArrayUtils;

public class LabelScoreArray {
	double activeScore[];
	public LabelScoreArray(int sz) {
		activeScore = new double[sz];
	}
	public void clear() {
		Arrays.fill(activeScore, 0);
	}
	public void add(int label, double score) {
		activeScore[label] += score;
	}
	public double get(int label) {
		return activeScore[label];
	}
	public void add(LabelScoreArray src) {
		ArrayUtils.arrayCopy(src.activeScore, activeScore);
	}
	public double[] toNativeArray() {
		return activeScore;
	}
}
