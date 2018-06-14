package iitb.sgl.data;

import java.util.Arrays;

public class Example {
	public float[] features;
	public int lbl;

	public Example(int numFeatures, int lbl) {
		this.features = new float[numFeatures];
		Arrays.fill(this.features, 0);
		this.lbl = lbl;
	}
}
