package iitb.sgl.data;

public class NodeFeatureGenerator {
	
	SocialGraph snGraph;
	int numNodeAttributes;
	int[] tmpNodeIdAttributes;
	NodeFeatured nfSocialGraph;
	public float nFVal;
	
	public NodeFeatureGenerator(SocialGraph snGraph) {
		this.snGraph = snGraph;
		assert(snGraph instanceof PokecGraph);
		nfSocialGraph = (NodeFeatured) snGraph;
		numNodeAttributes = nfSocialGraph.numNodeAttributes();
		tmpNodeIdAttributes = new int[numNodeAttributes];
		// Set Default Node Feature Value
		nFVal = 1.0f / snGraph.numLabels; 
	}
	
	public void getFeatureValues(int i, int featureIds[], int featureLbls[]) {
		// find the valid binary attributes for user i.		
		nfSocialGraph.getNodeAttributes(i,tmpNodeIdAttributes);
		int fIter = 0;
		for (int aIter = 0; aIter < tmpNodeIdAttributes.length && tmpNodeIdAttributes[aIter] >= 0; aIter++) {
			int aId = tmpNodeIdAttributes[aIter];
			for (int lbl = 0; lbl < snGraph.numLabels; lbl++) {
				featureIds[fIter] = aId * snGraph.numLabels + lbl;
				featureLbls[fIter] = lbl;
				++fIter;
			}
		}
		if (fIter < featureIds.length) {
			featureIds[fIter] = -1;
			featureLbls[fIter] = -1;
		}
	}
	
	public int numFeatures() {
		return numNodeAttributes * snGraph.numLabels;
	}
}
