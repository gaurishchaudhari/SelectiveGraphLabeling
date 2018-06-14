package iitb.shared.gm.inference;

import java.util.Arrays;

public class LabelScoreArraySparse extends LabelScoreArray {
	int labelPtrs[]; 
	int activeLabels[];
	int numActive = 0;
	public LabelScoreArraySparse(int sz) {
		super(sz);
		labelPtrs = new int[sz];
		Arrays.fill(labelPtrs, -1);
		activeLabels = new int[sz];
		// Added on 14-Nov Gaurish
		Arrays.fill(activeLabels, -1);
	}
	public void clear() {
		numActive = 0;
		for (int i = 0; i < activeLabels.length && activeLabels[i] >= 0; i++) {
			labelPtrs[activeLabels[i]] = -1;
			activeLabels[i] = -1;
			activeScore[i] = 0;
		}
	}
	public int getPutPtr(int label) {
		if (labelPtrs[label] < 0) { 
			numActive++;
			labelPtrs[label] = numActive-1;
			activeLabels[numActive-1] = label;
		}
		return labelPtrs[label];
	}
	public void add(int label, double score) {
		int ptr = getPutPtr(label);
		activeScore[ptr] += score;
	}
	
	public double get(int label) {
		if (labelPtrs[label] >= 0) {
			return activeScore[labelPtrs[label]];
		}
		return 0;
	}
	public void add(LabelScoreArraySparse src) {
		numActive = src.numActive;
		for (int i = 0; i < src.activeLabels.length && src.activeLabels[i] >= 0; i++) {
			labelPtrs[src.activeLabels[i]]=src.labelPtrs[src.activeLabels[i]];
			activeLabels[i] = src.activeLabels[i];
			activeScore[i] = src.activeScore[i];
		}
	}
	public double[] toNativeArray() {
		double ar[] = new double[labelPtrs.length];
		for (int i = 0; i < ar.length; i++) {
			if (labelPtrs[i] > 0) ar[i] = activeScore[labelPtrs[i]];
		}
		return ar;
	}
	
	public boolean isLabelActive(int lbl){
		return (labelPtrs[lbl] != -1);
	}
	
	public void sortTopKByScore(int K){
		boolean swapped;
		K = Math.min(numActive, K);
		for (int i = 0; i < K; ++i) {
			swapped = false;
			for(int j=numActive-2-i; j >=0 ; --j){
				if(activeScore[j] < activeScore[j+1]){
					//swap
					double tS = activeScore[j];
					activeScore[j] = activeScore[j+1];
					activeScore[j+1] = tS;
					
					int tL = activeLabels[j];
					activeLabels[j] = activeLabels[j+1];
					activeLabels[j+1] = tL;
					
					int tP = labelPtrs[activeLabels[j]];
					labelPtrs[activeLabels[j]] = labelPtrs[activeLabels[j+1]];
					labelPtrs[activeLabels[j+1]] = tP;
					
					swapped = true;
				}
			}
			if(!swapped) break;
		}
	}
	
	public int getNumActive(){
		return numActive;
	}
	
	public int getActiveLabel(int i) {
		return activeLabels[i];
	}
	
	public int[] getActiveLabelsArray() {
		return activeLabels;
	}
	
	@Override
	public String toString(){
		StringBuilder sb = new StringBuilder();
		for(int i=0; i<activeLabels.length; ++i)
			sb.append(activeLabels[i]+":"+activeScore[i]+", ");
		return sb.toString();
	}
}