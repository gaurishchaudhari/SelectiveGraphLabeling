package iitb.shared.gm;

import java.util.BitSet;

/**
 * Store Node Marginal Table of a node in sparse format. Space and Time Efficient.
 */

public class NodeMarginalSparse {
	
	private short[] labels;
	private float[] scores;
	private short idx;
	public BitSet activeLabels;
	public float defVal;
	public NodeMarginalSparse(int sz){
		scores = new float[sz];
		labels = new short[sz];
		idx = 0;
		activeLabels = new BitSet(sz);
		defVal = 0;
	}
	public boolean isLabelActive(int lbl){
		return activeLabels.get(lbl);
	}
	public int numActive(){
		return activeLabels.cardinality();
	}
	public float getScore(int lbl){
		if(activeLabels.get(lbl)){
			for(int i = 0; i < idx; ++i){
				if(labels[i] == lbl) return scores[i];
			}
		}
		return defVal;
	}
	public void setScore(int lbl, double score){
		activeLabels.set(lbl);
		labels[idx] = (short)lbl;
		scores[idx] = (float)score;
		++idx;
	}
}
