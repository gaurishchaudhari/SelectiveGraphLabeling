package iitb.sgl.data;

import java.util.Arrays;

/**
 * Generates edge features given a social graph.
 */

public class EdgeFeatureGenerator {
	
	public static float expectedUniDir = 1.50f;
	public static float expectedBiDir = 2.14f;
	
	/* Turn ON and OFF log edge features */
	public static boolean isLogFeature = true;
	
	static int lambda = 1;
	
	SocialGraph snGraph;
	
	public static class FeatureType {
		SocialGraph snGraph;
		float scale = 1.0f;
		float getFeature(int i, int j) {
			return 1.0f / scale; 
		}
	}
	
	public static class InFeatureType extends FeatureType {
		float getFeature(int i, int j) {
			if(!isLogFeature){
				int inDeg = snGraph.inDegree(i);
				return (inDeg > 0) ? 1.0f / inDeg : 0f;
			}
			else
				return (float)(1.0 / Math.log(snGraph.inDegree(i) + lambda)); 
		}
	}
	
	public static class OutFeatureType extends FeatureType {
		float getFeature(int i, int j) {
			if(!isLogFeature){
				int outDeg = snGraph.outDegree(j);
				return (outDeg > 0) ? 1.0f / outDeg : 0f;
			}
			else
				return (float)(1.0 / Math.log(snGraph.outDegree(j) + lambda));
		}
	}
	
	public static class TotalTotalFeatureType extends FeatureType {
		float getFeature(int i, int j) {
			if(!isLogFeature){
				int degI = snGraph.totalDegree(i);
				int degJ = snGraph.totalDegree(j);
				return (degI > 0 ? 1.0f / degI : 0) + (degJ > 0 ? 1.0f / degJ : 0); 
			}
			else
				return (float)(1.0 / Math.log(snGraph.totalDegree(i) + lambda) + 1.0 / Math.log(snGraph.totalDegree(j) + lambda));
		}
	}
	
	public static class InInFeatureType extends FeatureType {
		float getFeature(int i, int j) {
			if(!isLogFeature){
				int inDegI = snGraph.inDegree(i);
				int inDegJ = snGraph.inDegree(j);
				return (inDegI > 0 ? 1.0f / inDegI : 0) + (inDegJ > 0 ? 1.0f / inDegJ : 0);
			}
			else
				return (float)(1.0 / Math.log(snGraph.inDegree(i) + lambda) + 1.0 / Math.log(snGraph.inDegree(j) + lambda));
		}
	}
	
	public static class OutOutFeatureType extends FeatureType {
		float getFeature(int i, int j) {
			if(!isLogFeature){
				int outDegI = snGraph.outDegree(i);
				int outDegJ = snGraph.outDegree(j);
				return (outDegI > 0 ? 1.0f / outDegI : 0) + (outDegJ > 0 ? 1.0f / outDegJ : 0);
			}
			else
				return (float)(1.0 / Math.log(snGraph.outDegree(i) + lambda) + 1.0 /Math.log(snGraph.outDegree(j) + lambda));
		}
	}
	
	public static class IsObservedFeatureType extends FeatureType {
		float getFeature(int i, int j){
			return (snGraph.isNodeObserved(i) ? 1 : 0) + (snGraph.isNodeObserved(j) ? 1 : 0);	
			/*if((snGraph.observedDegree(i) == 1 && snGraph.observedNodes.get(j)) || (snGraph.observedDegree(j) == 1 && snGraph.observedNodes.get(i)))
				return 1;
			return 0;*/
		}
	}
	
	public static class PurityFeatureType extends FeatureType {
		short[] labelCount;
		float getFeature(int i, int j){
			Arrays.fill(labelCount, (short)0);
			int total = 0;
			for(int nbrNum = snGraph.graph.getNumNeighbours(i)-1; nbrNum>=0; --nbrNum){
				int nbr = snGraph.graph.getNeighbour(i, nbrNum);
				if(snGraph.observedNodes.get(nbr)) {
					++labelCount[snGraph.nodeLabels.get(nbr)];
					++total;
				}
			}
			float entropyI = 0;
			for(int k = 0; k < labelCount.length; ++k){
				if(labelCount[k] > 0)
					entropyI -= (labelCount[k]*1.0f/total)*Math.log(labelCount[k]*1.0f/total);
			}
			
			Arrays.fill(labelCount, (short)0);
			total = 0;
			for(int nbrNum = snGraph.graph.getNumNeighbours(j)-1; nbrNum>=0; --nbrNum){
				int nbr = snGraph.graph.getNeighbour(j, nbrNum);
				if(snGraph.observedNodes.get(nbr)) {
					++labelCount[snGraph.nodeLabels.get(nbr)];
					++total;
				}
			}
			float entropyJ = 0;
			for(int k = 0; k < labelCount.length; ++k){
				if(labelCount[k] > 0) 
					entropyJ -= (labelCount[k]*1.0f/total)*Math.log(labelCount[k]*1.0f/total);
			}
			return 0.5f * (2 - ((entropyI + entropyJ)/(float)Math.log(snGraph.numLabels)));
		}
	}
	
	private int useNumFeatures = 7;
	
	FeatureType featureTypes[] = null;
	public EdgeFeatureGenerator(SocialGraph snGraph) {
		this.snGraph = snGraph;
		if(snGraph instanceof TwitterGraph){
			useNumFeatures = 7;
			featureTypes = new FeatureType[useNumFeatures];
			featureTypes[0] = new FeatureType();
			featureTypes[1] = new FeatureType();
			featureTypes[2] = new InFeatureType();
			featureTypes[3] = new OutFeatureType();
			featureTypes[4] = new TotalTotalFeatureType();
			featureTypes[5] = new InInFeatureType();
			featureTypes[6] = new OutOutFeatureType();
		}
		else{
			useNumFeatures = 9;
			featureTypes = new FeatureType[useNumFeatures];
			featureTypes[0] = new FeatureType();
			featureTypes[1] = new FeatureType();
			featureTypes[2] = new InFeatureType();
			featureTypes[3] = new OutFeatureType();
			featureTypes[4] = new TotalTotalFeatureType();
			featureTypes[5] = new InInFeatureType();
			featureTypes[6] = new OutOutFeatureType();
			featureTypes[7] = new PurityFeatureType();
			featureTypes[8] = new IsObservedFeatureType();
			((PurityFeatureType)featureTypes[7]).labelCount = new short[snGraph.numLabels];
		}
		
		featureTypes[0].scale = expectedUniDir;
		featureTypes[1].scale = expectedBiDir;
		for(int i = 0; i < featureTypes.length; ++i)
			featureTypes[i].snGraph = snGraph;
	}
	
	public double getFeatureValue(int i, int j, int featureNumber) {
		/* featureNumber 1 is redundant. Instead put featureNumber 1 for R-Friends */	
		if(featureNumber == 0 || featureNumber == 1 || featureNumber==4 || featureNumber == 7 || featureNumber == 8){
			return featureTypes[featureNumber].getFeature(i, j);
		}
		if(snGraph.inEdgesMap.get(i).contains(j)){
			if(featureNumber==2 || featureNumber==3)
				return featureTypes[featureNumber].getFeature(i, j);
		}
		else if(snGraph.outEdgesMap.get(j).contains(i)){
			if(featureNumber==2 || featureNumber==3)
				return featureTypes[featureNumber].getFeature(j, i);
		}
		else if(snGraph.inOutEdgesMap.get(i).contains(j) || snGraph.inOutEdgesMap.get(j).contains(i)){
			if(featureNumber==5 || featureNumber==6)
				return featureTypes[featureNumber].getFeature(i, j);
		}
		return 0;
	}
	
	public int numFeatures() {
		return featureTypes.length;
	}
}
