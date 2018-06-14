package iitb.shared;

import java.util.*;

public class ArrayUtils{
	
	public static int maxIndex(double [] scores) {
		int maxI=0;
		for (int i = 0; i < scores.length; i++) {
			if (scores[i] > scores[maxI]) {
				maxI=i;
			}
		}
		return maxI;
	}
	
	public static int minIndex(double [] scores, BitSet bitset) {
		int maxI=-1;
		for (int i = 0; i < scores.length; i++) {
			if (bitset.get(i)) continue;
			if (maxI == -1 || scores[i] < scores[maxI]) {
				maxI=i;
			}
		}
		return maxI;
	}
	
	public static int maxIndex(float [] scores) {
		int maxI=0;
		for (int i = 0; i < scores.length; i++) {
			if (scores[i] > scores[maxI]) {
				maxI=i;
			}
		}
		return maxI;
	}

	public static float sum(float[] fs) {
		float sumI=0;
		for (int i = 0; i < fs.length; i++) {
			sumI += fs[i];
		}
		return sumI;
	}
	
	public static double norm2Square(double ar[]) {
		double v = 0;
		for (int f = 0; f < ar.length; f++)
			v += ar[f]*ar[f];
		return v;
	}
	
	public static double norm2(double ar[]) {
		return Math.sqrt(norm2Square(ar));
	}
	
	public static double norm2(float[] ar) {
		return Math.sqrt(norm2Square(ar));
	}
	
	public static double norm2Square(float[] ar) {
		double v = 0;
		for (int f = 0; f < ar.length; f++)
			v += ar[f]*ar[f];
		return v;
	}
	
	public static float max(float[] fs) {
		float max=fs[0];
		for (int i = 1; i < fs.length; i++) {
			max = Math.max(max,fs[i]);
		}
		return max;
	}
	
	public static int min(int[] fs) {
		if (fs.length == 0) return 0;
		int minV=fs[0];
		for (int i = 1; i < fs.length; i++) {
			minV = Math.min(minV,fs[i]);
		}
		return minV;
	}
	
	public static float min(float[] fs) {
		float minV=fs[0];
		for (int i = 1; i < fs.length; i++) {
			minV = Math.min(minV,fs[i]);
		}
		return minV;
	}
	public static double[] parseDoubleArray(String str) {
		String[] a = (str.contains(",")?str.split(","):str.split(" "));
		double c[] = new double[a.length];
		for(int i=0;i < a.length;i++) {
			c[i] = Double.parseDouble(a[i]);
		}
		return c;
	}
	
	public static int[] parseIntArray(String str) {
		str = str.trim();
		if (str.length()==0) return new int[0];
		String[] a = str.split(" ");
		int c[] = new int[a.length];
		for(int i=0;i < a.length;i++) {
			c[i] = Integer.parseInt(a[i]);
		}
		return c;
	}
	
	public static float[] parseFloatArray(String str) {
		String[] a = str.split(" ");
		float c[] = new float[a.length];
		for(int i=0;i < a.length;i++) {
			c[i] = Float.parseFloat(a[i]);
		}
		return c;
	}
	
	public static BitSet smallestK(double vals[], int k) {
		BitSet bitset = new BitSet(vals.length);
		for (int i = 0; i < k; i++) {
			int pos = minIndex(vals,bitset);
			bitset.set(pos);
		}
		return bitset;
	}
	
	public static int find(int[] array, int i) {
		for (int j = 0; j < array.length; j++) {
			if (array[j]==i) return j;
		}
		return -1;
	}
	
	public static int max(int[] array) {
		if (array==null || array.length==0) return 0;
		int max = array[0];
		for (int i = 0; i < array.length; i++) {
			max = Math.max(max, array[i]);
		}
		return max;
	}
	
	public static boolean isSorted(int[] is) {
		for (int i = 1; i < is.length; i++) {
			if (is[i-1] > is[i])
				return false;
		}
		return true;
	}
	
	public static int sum(int[] labeling) {
		int s = 0;
		for (int i = 0; i < labeling.length; i++) {
			s += labeling[i];
		}
		return s;
	}
	
	public static double dotProduct(double[] vec1, double[] vec2) {
		double dotp = 0;
		for (int i = 0; i < vec2.length; i++) {
			dotp += vec1[i]*vec2[i];
		}
		return dotp;
	}
	public static double dotProduct(double[] vec1, float[] vec2) {
		double dotp = 0;
		for (int i = 0; i < vec2.length; i++) {
			dotp += vec1[i]*vec2[i];
		}
		return dotp;
	}
	
	public static double normInfty(double[] array) {
		if (array==null || array.length==0) return 0;
		double max = array[0];
		for (int i = 0; i < array.length; i++) {
			max = Math.max(max, Math.abs(array[i]));
		}
		return max;
	}
	
	public static double norm2Square(double[] ar, int from, int to) {
		double v = 0;
		to = Math.min(to, ar.length);
		for (int f = from; f < to; f++)
			v += ar[f]*ar[f];
		return v;
	}
	
	public static String concat(String[] labels, String del) {
		String str = "";
		for (int i = 0; i < labels.length; i++) {
			if (i > 0) str += del;
			str += labels[i];
		}
		return str;
	}
	
	public static int[] subsetEntries(int[] clique, int[] labels) {
		int a[] = new int[clique.length];
		for (int i = 0; i < a.length; i++) {
			a[i] = labels[clique[i]];
		}
		return a;
	}
	
	public static double sum(double[] labeling) {
		double s = 0;
		for (int i = 0; i < labeling.length; i++) {
			s += labeling[i];
		}
		return s;
	}
	
	public static void arrayCopy(double[] src, double[] dest) {
		System.arraycopy(src, 0, dest, 0, src.length);
	}
	
	public static void arrayCopy(double[] src, double[] dest, int firstN) {
		System.arraycopy(src, 0, dest, 0, firstN);
	}
	
	public static double KLDistance(double dist1[], double dist2[]) {
		double dist = 0;
		for (int i = 0; i < dist1.length; i++) {
			dist += dist1[i]*Math.log(dist1[i]/dist2[i]);
		}
		return dist;
	}
	
	public static double L1Distance(double[] dist1, double[] dist2) {
		double dist = 0;
		for (int i = 0; i < dist1.length; i++) {
			dist += Math.abs(dist1[i]-dist2[i]);
		}
		return dist;
	}
	
	public static void arrayCopy(float[] src, float[] dest) {
		System.arraycopy(src, 0, dest, 0, src.length);
	}
	
	public static void arrayCopy(int[] src, int[] dest) {
		System.arraycopy(src, 0, dest, 0, src.length);
	}
	
	public static String toStrArr(double []arr){
		StringBuilder sb = new StringBuilder("[");
		int i;
		for(i=0; i<arr.length-1; ++i)
			sb.append(String.format("%.4f", arr[i])+", ");
		sb.append(String.format("%.4f", arr[i])+"]");
		return sb.toString();
	}
}
