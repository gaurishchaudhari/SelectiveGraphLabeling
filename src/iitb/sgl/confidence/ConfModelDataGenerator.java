package iitb.sgl.confidence;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.Arrays;
import java.util.BitSet;
import java.util.Random;
import java.util.StringTokenizer;

import iitb.sgl.data.EdgeFeatureGenerator;
import iitb.sgl.data.Example;
import iitb.sgl.data.NodeFeatureGenerator;
import iitb.sgl.data.NodeFeatured;
import iitb.sgl.data.SocialGraph;
import iitb.sgl.inference.SGInference;
import iitb.sgl.inference.SGModelUtils;
import iitb.sgl.learning.MiniGraphModelFactory;
import iitb.sgl.learning.MiniGraphModelFactory.MiniGM;
import iitb.sgl.main.RunGraphLabeler;
import iitb.shared.ArrayUtils;
import iitb.shared.gm.NodeMarginalSparse;
import iitb.shared.gm.PottsPotential;
import iitb.shared.gm.UDGM;
import iitb.shared.gm.inference.ScalableMessagePassing.Method;
import weka.core.Attribute;
import weka.core.FastVector;
import weka.core.Instance;
import weka.core.Instances;

/**
 * Our Approach: Training and Test Set Generator.
 * Performs inference on full graph and then mini-graph inference at each observed node
 * to generate positive and negative examples in training and test data.
 * ALso provides some meta-data in non-ARFF files produced as part of output.
 * @author gaurish
 */

public class ConfModelDataGenerator {
	final static String bag = "/mnt/bag/wwt/twitter/";
	static int numClassifierFeatures = 22;

	String trainingFile;
	String testFile;
	
	SocialGraph snGraph;
	NodeFeatureGenerator nodeFeatureGenerator;
	EdgeFeatureGenerator edgeFeatureGenerator;
	int numNodeFeatures;
	int numEdgeFeatures;
	int numFeatures;
	int[] nodeFeatureIds;
	int[] nodeFeatureLabels;
	int[] labeling;
	Example[] data;
	int[] miniGraphLabels;
	float[] miniGraphConfidences;
	double[] modelWeights;
	SGInference solution;
	NodeMarginalSparse[] nodeMarginals;
	MiniGM miniGM;
	SGInference mSolution;
	double[] marginal;
	BitSet candidateLabels;
	
	Instances trainingSet;
	Instances testSet;
	FastVector fvWekaAttributes;
	
	public ConfModelDataGenerator(SocialGraph snGraph, String modelWeightsFile, String trainingFile, String testFile) throws Exception{
		this.snGraph = snGraph;
		this.nodeFeatureGenerator = null;
		this.numNodeFeatures = 0;
		if(snGraph instanceof NodeFeatured){
			this.nodeFeatureGenerator = new NodeFeatureGenerator(snGraph);
			this.numNodeFeatures = nodeFeatureGenerator.numFeatures();
			this.nodeFeatureIds = new int[numNodeFeatures];
			this.nodeFeatureLabels = new int[nodeFeatureIds.length];
		}
		this.edgeFeatureGenerator = new EdgeFeatureGenerator(snGraph);
		this.numEdgeFeatures = edgeFeatureGenerator.numFeatures();
		this.numFeatures = numNodeFeatures + numEdgeFeatures;
		this.labeling = new int[snGraph.graph.getNumNodes()];
		Arrays.fill(labeling, -1);
		for(int i=snGraph.graph.getNumNodes()-1; i>=0; --i)
			labeling[i] = snGraph.isNodeObserved(i) ? snGraph.nodeLabels.get(i) : -1;
		this.data = new Example[snGraph.graph.getNumNodes()];
		this.marginal = new double[snGraph.numLabels];
		this.candidateLabels = new BitSet();
		
		if(modelWeightsFile != null){
			this.modelWeights = RunGraphLabeler.loadWeightsFromFile(modelWeightsFile);	// Read Weights Learnt for Graphical Model
			System.out.println("Optimal Graphical Model Feature Weights = "+ArrayUtils.toStrArr(modelWeights));
		}
		this.trainingFile = trainingFile;
		this.testFile = testFile;
	}
	
	public void evaluate() throws Exception{
		int total = 0;
		int correct = 0;
		Random r = new Random(5);
		for(int user=snGraph.graph.getNumNodes()-1; user>=0; --user){
			if(snGraph.isNodeObserved(user)) continue;			
			int trueIndex = snGraph.nodeLabels.get(user);
			int predictedIndex = solution.getSolution(user);
			if(predictedIndex == -1)
				predictedIndex = r.nextInt(snGraph.numLabels);		
			labeling[user] = predictedIndex;
			++total;
			if (trueIndex == predictedIndex) ++correct;
		}
		System.out.println("Full-Graph Inference Accuracy : " + correct + " / " + total + " = "+ (correct * 100.0) / total);
	}
	
	public void createTrainSet() throws Exception{
		int index = 0;
		int numPosInst = 0;
		int numNegInst = 0;
		miniGraphLabels = new int[data.length];
		miniGraphConfidences = new float[data.length];
		Arrays.fill(miniGraphLabels, -99);
		Arrays.fill(miniGraphConfidences, -9999);
		for(int user=snGraph.observedNodes.nextSetBit(0); user>=0; user=snGraph.observedNodes.nextSetBit(user+1)){
			addExample(user);
			if(data[user] == null)
				continue;
			if(data[user].lbl == 1) ++numPosInst;
			else ++numNegInst;
			if(index % 1000 == 0)
				System.out.println(index+" Done!");
			++index;
		}
		System.gc();
		System.out.println("# Positive Examples "+numPosInst +" ; # Negative Examples "+numNegInst);		
	}
	
	public void addExample(int user) throws Exception{
		MiniGraphModelFactory miniGMFactory = new MiniGraphModelFactory();
		miniGM = miniGMFactory.createMiniGraph(snGraph, user, 2);
		if(modelWeights != null)
			MiniGraphModelFactory.setMiniGraphPotentials(user, snGraph, nodeFeatureGenerator, edgeFeatureGenerator, miniGM, modelWeights, null);
		else
			MiniGraphModelFactory.setMiniGraphNodePotentials(snGraph, nodeFeatureGenerator, miniGM, user, null , null);
		
		mSolution = new SGInference(miniGM.gModel, 0, Method.BP, false, miniGM.observedNodes);
		mSolution.getTopK(1);
		int predictedLabel = mSolution.getSolution(miniGM.nodeIndex.get(user));
		int trueLabel = snGraph.nodeLabels.get(user);
		int lbl = (predictedLabel == trueLabel) ? 1 : 0;
		 
		if(predictedLabel == -1) return;
		
		miniGraphLabels[user] = predictedLabel;
		miniGraphConfidences[user] = (float)mSolution.getConfidence(miniGM.nodeIndex.get(user));
		data[user] = new Example(numClassifierFeatures, lbl);
		
		getFeatures(user, data[user]);
	}
	
	public void getFeatures(int userId, Example dt){
		BitSet lvl1Nbrs = new BitSet(snGraph.graph.getNumNodes());
		BitSet lvl2Nbrs = new BitSet(snGraph.graph.getNumNodes());
		for(int j=snGraph.graph.getNumNeighbours(userId)-1; j>=0; --j){
			int nbr = snGraph.graph.getNeighbour(userId, j);
			lvl1Nbrs.set(nbr);
		}
		for(int nbr2 = lvl1Nbrs.nextSetBit(0); nbr2 >= 0; nbr2 = lvl1Nbrs.nextSetBit(nbr2+1)){
			for(int j=snGraph.graph.getNumNeighbours(nbr2)-1; j>=0; --j){
				int nbr = snGraph.graph.getNeighbour(nbr2, j);
				if(nbr != userId && !lvl1Nbrs.get(nbr))
					 lvl2Nbrs.set(nbr);
			}
		}
		int predictedLbl = snGraph.isNodeObserved(userId)
				? mSolution.getSolution(miniGM.nodeIndex.get(userId))
				: labeling[userId];
				
		double predConfidence = snGraph.isNodeObserved(userId)
				? mSolution.getConfidence(miniGM.nodeIndex.get(userId))
				: solution.getConfidence(userId);
		int numObsNodesLvl1 = 0;
		int numObsNodesLblLvl1 = 0;
		int numObsNodesLblLvl2 = 0;
		int numUnObsNodesLvl1 = 0;
		int numUnObsNodesLblLvl1 = 0;
		
		short[] obsLblCount = new short[snGraph.numLabels];
		short[] allLblCount = new short[snGraph.numLabels];
		float sumMarginalPredLbl = 0;
		float sumMarginalLbl = 0;
		float sumMarginalWinningLbl = 0;
		
		for(int nbr1=lvl1Nbrs.nextSetBit(0); nbr1>=0; nbr1=lvl1Nbrs.nextSetBit(nbr1+1)){
			if(snGraph.isNodeObserved(nbr1)){
				++numObsNodesLvl1;
				++obsLblCount[labeling[nbr1]];
				if(labeling[nbr1] == predictedLbl)
					++numObsNodesLblLvl1;
			}
			else{
				++numUnObsNodesLvl1;
				if(labeling[nbr1] == predictedLbl)
					++numUnObsNodesLblLvl1;
			}
			++allLblCount[labeling[nbr1]];
			if(labeling[nbr1] == predictedLbl)
				sumMarginalPredLbl += nodeMarginals[nbr1].getScore(predictedLbl);
			sumMarginalLbl += nodeMarginals[nbr1].getScore(predictedLbl);
			sumMarginalWinningLbl += solution.getConfidence(nbr1);
		}
		for(int nbr2=lvl2Nbrs.nextSetBit(0); nbr2>=0; nbr2=lvl2Nbrs.nextSetBit(nbr2+1)){
			if(snGraph.isNodeObserved(nbr2)){
				//++numObsNodesLvl2;
				if(labeling[nbr2] == predictedLbl)
					++numObsNodesLblLvl2;
			}
		}
		
		int numOnlyFlwLbl = 0;
		int numOnlyFlw = 0;
		int numOnlyFrndLbl = 0;
		int numOnlyFrnd = 0;
		int numRFrndLbl = 0;
		int numRFrnd = 0;
		for(int nbr1=lvl1Nbrs.nextSetBit(0); nbr1>=0; nbr1=lvl1Nbrs.nextSetBit(nbr1+1)){
			if(snGraph.inEdgesMap.get(userId).contains(nbr1)){
				++numOnlyFlw;
				if(labeling[nbr1] == predictedLbl) ++numOnlyFlwLbl;
			}
			else if(snGraph.inEdgesMap.get(nbr1).contains(userId)){
				++numOnlyFrnd;
				if(labeling[nbr1] == predictedLbl) ++numOnlyFrndLbl;
			}
			else if(snGraph.inOutEdgesMap.get(userId).contains(nbr1) || snGraph.inOutEdgesMap.get(nbr1).contains(userId)){
				++numRFrnd;
				if(labeling[nbr1] == predictedLbl) ++numRFrndLbl;
			}
		}
		
		dt.features[0] = numObsNodesLblLvl1 / 1.5f;
		dt.features[1] = numObsNodesLblLvl2 / 1.5f;
		dt.features[2] = (float)predConfidence;
		dt.features[3] = sumMarginalPredLbl / (snGraph.graph.getNumNeighbours(userId) + 0.5f);
		dt.features[4] = sumMarginalPredLbl / (sumMarginalLbl + 0.2f);
		dt.features[5] = getPurity(obsLblCount);
		dt.features[6] = getPurity(allLblCount);
		dt.features[7] = sumMarginalPredLbl / (sumMarginalWinningLbl + 0.2f);
		dt.features[8] = getTop2Difference(userId);
		dt.features[9] = (snGraph.graph.getNumNeighbours(userId) == 1 && snGraph.isNodeObserved(snGraph.graph.getNeighbour(userId, 0))) ? 1 : 0;
		dt.features[10] = numObsNodesLvl1 / (snGraph.totalDegree(userId) + 0.5f);
		
		dt.features[11] = (float)(numOnlyFlwLbl / 1.5f);
		dt.features[12] = (float)(numOnlyFlwLbl / (numOnlyFlw +  0.5f));
		dt.features[13] = (float)(numOnlyFrndLbl / 1.5f);
		dt.features[14] = (float)(numOnlyFrndLbl / (numOnlyFrnd +  0.5f));
		dt.features[15] = (float)(numRFrndLbl / 2.14f);
		dt.features[16] = (float)(numRFrndLbl / (numRFrnd +  0.5f));
		
		dt.features[17] = (numObsNodesLblLvl1 + numUnObsNodesLblLvl1) / 2.0f;
		dt.features[18] = (numObsNodesLblLvl1 + numUnObsNodesLblLvl1) / (numObsNodesLvl1 + numUnObsNodesLvl1 + 0.5f);
		dt.features[19] = (numObsNodesLblLvl1 + 0.5f * numUnObsNodesLblLvl1) / (numObsNodesLvl1 + 0.5f * numUnObsNodesLvl1 + 0.5f);
		dt.features[20] = numObsNodesLblLvl1 / (numObsNodesLvl1 + 0.5f);
		dt.features[21] = numObsNodesLblLvl1 > 0 ? 1 : 0;
	}
	
	public float getPurity(short[] lblCount){
		short total = 0;
		for (int i = 0; i < lblCount.length; i++) 
			total += lblCount[i];
		float entropy = 0;
		for(int k = 0; k < lblCount.length; ++k){
			if(lblCount[k] > 0) 
				entropy += -(lblCount[k]*1.0f/total)*Math.log(lblCount[k]*1.0f/total);
		}
		return 1f - (entropy / (float)Math.log(lblCount.length));
	}
	
	public float getTop2Difference(int userId){
		float diff = 0;
		if(snGraph.isNodeObserved(userId)){
			candidateLabels.clear();
			float pr1 = (float)mSolution.getNodeSumMarginal(miniGM.nodeIndex.get(userId), marginal, candidateLabels);
			float pr2 = pr1; 
			for(int lbl = candidateLabels.nextSetBit(0); lbl >= 0; lbl = candidateLabels.nextSetBit(lbl+1)) {
				float val = (float)marginal[lbl];
				if(val > pr1){
					pr2 = pr1;
					pr1 = val;
				}
				else if(val > pr2)
					pr2 = val;
			}
			diff = (float)(pr1 - pr2);
		}
		else{
			candidateLabels = nodeMarginals[userId].activeLabels;
			float pr1 = nodeMarginals[userId].defVal;
			float pr2 = pr1; 
			for(int lbl = candidateLabels.nextSetBit(0); lbl >= 0; lbl = candidateLabels.nextSetBit(lbl+1)) {
				float val = nodeMarginals[userId].getScore(lbl);
				if(val > pr1){
					pr2 = pr1;
					pr1 = val;
				}
				else if(val > pr2)
					pr2 = val;
			}
			diff = (float)(pr1 - pr2);
		}
		return diff;
	}
	
	public void createTestSet(){
		for(int user=0; user < snGraph.graph.getNumNodes(); ++user){
			if(snGraph.isNodeObserved(user))
				continue;
			int predictedLbl = snGraph.nodeLabels.get(user) == labeling[user] ? 1 : 0;
			data[user] = new Example(numClassifierFeatures, predictedLbl);
			getFeatures(user, data[user]);
		}
	}
	
	public void generateFeatures() throws Exception{		
		UDGM model = new UDGM(snGraph.graph, snGraph.numLabels, new PottsPotential(0.5, 0), 2);
		SGModelUtils.setPotentials(snGraph, model, modelWeights, nodeFeatureGenerator, edgeFeatureGenerator);
		solution = new SGInference(model, 1, Method.BP, false, snGraph.observedNodes);
		solution.getTopK(1);
			
		nodeMarginals = new NodeMarginalSparse[snGraph.graph.getNumNodes()];
		for(int node=snGraph.graph.getNumNodes()-1; node>=0; --node){
			candidateLabels.clear();		
			double dV = solution.getNodeSumMarginal(node, marginal, candidateLabels);
			nodeMarginals[node] = new NodeMarginalSparse(candidateLabels.cardinality());
			nodeMarginals[node].defVal = (float)dV;
			for(int lbl=candidateLabels.nextSetBit(0); lbl>=0; lbl=candidateLabels.nextSetBit(lbl+1)){
				nodeMarginals[node].activeLabels.set(lbl);
				nodeMarginals[node].setScore(lbl, marginal[lbl]);
			}
		}
		evaluate();
		createTrainSet();
		createTestSet();
		
		// Write train and test files
		writeData(trainingFile, testFile, data);
	}

	public void writeData(String trainFile, String testFile, Example[] data) throws Exception{
		FileWriter fwTrn = new FileWriter(trainFile);
		FileWriter fwTst = new FileWriter(testFile);
		for(int ix = 0; ix < snGraph.graph.getNumNodes(); ++ix){
			if(data[ix] != null){
				if(snGraph.isNodeObserved(ix)){
					fwTrn.write(ix + " " + snGraph.nodeLabels.get(ix) + " " + miniGraphLabels[ix] + " " + 
							miniGraphConfidences[ix] + " "+ snGraph.isNodeObserved(ix) + " "+data[ix].lbl+" ");
					for(int f = 0; f < data[ix].features.length; ++f)
						fwTrn.write(data[ix].features[f]+" ");
					fwTrn.write("\n");
				}
				else{
					fwTst.write(ix + " " + snGraph.nodeLabels.get(ix) + " " + solution.getSolution(ix) + " " + 
							solution.getConfidence(ix) + " "+ snGraph.isNodeObserved(ix) + " "+data[ix].lbl+" ");
					for(int f = 0; f < data[ix].features.length; ++f)
						fwTst.write(data[ix].features[f]+" ");
					fwTst.write("\n");
				}
			}
		}
		fwTrn.close();
		fwTst.close();
	}
	
	public void generateWekaDataset() throws Exception {
		FastVector classLbl = new FastVector(2);
		classLbl.addElement("positive");
		classLbl.addElement("negative");
		Attribute classAttribute = new Attribute("classLbl", classLbl);

		fvWekaAttributes = new FastVector(
				ConfModelDataGenerator.numClassifierFeatures + 1);
		for (int fNum = 0; fNum < ConfModelDataGenerator.numClassifierFeatures; ++fNum)
			fvWekaAttributes.addElement(new Attribute("f" + fNum));
		fvWekaAttributes.addElement(classAttribute);

		trainingSet = new Instances("twiterTrain", fvWekaAttributes, 10);
		trainingSet.setClassIndex(fvWekaAttributes.size() - 1); // Set class index

		String line;
		BufferedReader br = new BufferedReader(new FileReader(trainingFile));
		while ((line = br.readLine()) != null) {
			StringTokenizer st = new StringTokenizer(line);
			Instance iExample = new Instance(fvWekaAttributes.size());
			st.nextToken();	// USER-ID		
			st.nextToken();	// TRUE-LABEL
			st.nextToken();	// PREDICTED-LABEL
			st.nextToken();	// CONFIDENCE
			st.nextToken();	// Is-OBSERVED
			int lbl = Integer.parseInt(st.nextToken());
			int fNum = 0;
			while (st.hasMoreTokens()) {
				iExample.setValue(fNum, Float.parseFloat(st.nextToken()));
				++fNum;
			}
			iExample.setValue(classAttribute, lbl == 1 ? "positive" : "negative");
			trainingSet.add(iExample);
		}
		br.close();

		System.out.println("Number of Attributes = " + trainingSet.numAttributes());
		System.out.println("Training Set Size = " + trainingSet.numInstances());

		testSet = new Instances("twiterTest", fvWekaAttributes, 10);
		testSet.setClassIndex(fvWekaAttributes.size() - 1);
		
		br = new BufferedReader(new FileReader(testFile));
		while ((line = br.readLine()) != null) {
			StringTokenizer st = new StringTokenizer(line);
			Instance iExample = new Instance(fvWekaAttributes.size());
			st.nextToken();	// USER-ID		
			st.nextToken();	// TRUE-LABEL
			st.nextToken();	// PREDICTED-LABEL
			st.nextToken();	// CONFIDENCE
			st.nextToken();	// Is-OBSERVED
			int lbl = Integer.parseInt(st.nextToken());
			int fNum = 0;
			while (st.hasMoreTokens()) {
				iExample.setValue(fNum, Float.parseFloat(st.nextToken()));
				++fNum;
			}
			iExample.setValue(classAttribute, lbl == 1 ? "positive" : "negative");
			testSet.add(iExample);
		}
		br.close();

		System.out.println("Test Set Size = " + testSet.numInstances());

		// Write ARFF
		FileWriter fw = new FileWriter(trainingFile+".arff");
		fw.write(trainingSet.toString());
		fw.close();
		fw = new FileWriter(testFile+".arff");
		fw.write(testSet.toString());
		fw.close();
		
		System.out.println(trainingFile+" AND "+testFile+" ARFF files are Written !");
	}
}

