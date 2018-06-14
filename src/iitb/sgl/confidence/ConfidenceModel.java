package iitb.sgl.confidence;

import java.io.FileReader;
import java.io.FileWriter;

import weka.classifiers.Classifier;
import weka.classifiers.Evaluation;
import weka.classifiers.functions.Logistic;
import weka.classifiers.functions.SMO;
import weka.core.Instance;
import weka.core.Instances;

/**
 * Run Second-Stage Classifier on the training data to estimate confidence values of predictions.
 * @author gaurish
 */

public class ConfidenceModel {

	public Classifier getSMO(String kernel) throws Exception{
		Classifier cModel = new SMO();
		double C = 1000.0;
		String[] options = null;
		if(kernel.equals("RBF")){
			double gamma = 0.1;
			options = weka.core.Utils.splitOptions("-M -C "+ C +" -L 0.001 -P 1.0E-12 -N 0 -V -1 -W 1 -K \"weka.classifiers.functions.supportVector.RBFKernel -C 250007 -G "+gamma+"\"");
		}
		else if(kernel.equals("Poly")){
			double exponent = 1.0;
			options = weka.core.Utils.splitOptions("-M -C "+ C +" -L 0.001 -P 1.0E-12 -N 0 -V -1 -W 1 -K \"weka.classifiers.functions.supportVector.PolyKernel -C 250007 -E "+exponent+"\"");			
		}
		cModel.setOptions(options);
		return cModel;
	}
	public Classifier getLogistic() throws Exception {
		Classifier cModel = new Logistic();
		return cModel;
	}

	
	public ConfidenceModel(String trainARFF, String testARFF, String outputFile) throws Exception {
		Instances trainingSet = new Instances(new FileReader(trainARFF));
		trainingSet.setClassIndex(ConfModelDataGenerator.numClassifierFeatures);
		
		for(int i=0; i<trainingSet.numInstances(); ++i){
			
			// Pokec
			for(int k: new int[]{7,8,3,4,2,19,5,22,9,20,6,1,10,13}){
				trainingSet.instance(i).setValue(k, 0);
			}
			
			/* //Twitter
			for(int k: new int[]{2,5,7,8}){
				trainingSet.instance(i).setValue(k, 0);
			}*/
			
			/*Change Instance Weights
			if(trainingSet.instance(i).value(2) >  0.3){	// fNum_2: marginal prob of prediction
				trainingSet.instance(i).setWeight(100);
			}*/	
		}
		
		//Classifier cModel = getSMO("Poly");
		Classifier cModel = getLogistic();
		cModel.buildClassifier(trainingSet);
		System.out.println("Using Model " + cModel.toString());

		Evaluation eTest = new Evaluation(trainingSet);
		eTest.evaluateModel(cModel, trainingSet);
		String strSummary = eTest.toSummaryString();
		System.out.println("Evaluation Summary");
		System.out.println(strSummary);

		// Get the confusion matrix
		double[][] cmMatrix = eTest.confusionMatrix();
		System.out.println("Confusion Matrix");
		for (int ix = 0; ix < cmMatrix.length; ++ix) {
			for (int iy = 0; iy < cmMatrix[ix].length; ++iy)
				System.out.print(cmMatrix[ix][iy] + "\t");
			System.out.println();
		}
		
		// Apply classifier on test set
		Instances testSet = new Instances(new FileReader(testARFF));
		testSet.setClassIndex(ConfModelDataGenerator.numClassifierFeatures);		
		int correct = 0;
		int total = 0;
		FileWriter fw = new FileWriter(outputFile);
		for(int ix = 0; ix < testSet.numInstances(); ++ix){
			Instance iUse = testSet.instance(ix);
			double[] distr = cModel.distributionForInstance(iUse);
			int lbl = (int)iUse.classValue();
			++total;
			lbl = (lbl == 1) ? 0 : 1;
			if(lbl == 1) ++correct;
			if(ix % 100000==0) System.out.println(distr[0]+" "+lbl);
			fw.write(distr[0]+" "+lbl+"\n");
		}
		fw.close();
		System.out.println("Binary Clasifier Test Accuracy = "+correct+" / "+total+" = "+(correct*100.0/total));
	}
	
	public static void main(String[] args) throws Exception {
		String path = "OurApproach/";
		String trainARFF = ConfModelDataGenerator.bag + path + "NCL_10_2_trainingSet.arff";
		String testARFF = ConfModelDataGenerator.bag + path + "NCL_10_2_testSet.arff";
		String outputFile = ConfModelDataGenerator.bag + "NCL_BinLR_10_r2";
		
//		String trainARFF = args[0];
//		String testARFF = args[1];
//		String outputFile = args[2];
		
		new ConfidenceModel(trainARFF, testARFF, outputFile);
		
		//PR_PlotGenerator.generatePRData(outputFile, null, true, null, null, null, 10);
		System.out.println("Done !");
	}
}
