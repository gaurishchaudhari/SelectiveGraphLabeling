package iitb.shared.gm.inference;

import java.io.Serializable;
import java.util.Arrays;
import java.util.BitSet;

import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.list.array.TIntArrayList;
import iitb.shared.RobustMath;
import iitb.shared.gm.GraphicalModel;
import iitb.shared.gm.inference.ScalableMessagePassing.FwdBwdGraph.EdgeIter;
import iitb.shared.gm.inference.SparseScalableMP.SparseMessageCache.MessageIterator;
import iitb.shared.gm.NodePotentialSparse;
import sun.reflect.generics.reflectiveObjects.NotImplementedException;

/**
 * Sparse Version of ScalableMessagePassing
 * @author gaurish
 */

public class SparseScalableMP extends ScalableMessagePassing {
	private static final long serialVersionUID = 1L;

	public SparseScalableMP() {
		super();
	}	
	public SparseScalableMP(boolean trws, boolean maxProd){		
		super(trws,maxProd);
	}	
	protected SparseMessageCache messageCacheSparse;	// Created a new reference for superclass's messageCache reference which can hold its subclass
	protected BitSet labelSet;
	protected MessageIterator msgIter;
	protected LabelScoreArraySparse tempBufSrcS;
	protected LabelScoreArraySparse tempBufDstS;
	protected NodePotentialSparse nPotS;
	protected double[] score = new double[1];
	
	static final boolean toBePruned = true;
	static final double THRESHOLD = 10e-5;
	
	public SparseScalableMP(GraphicalModel model) {
		this();
		if (model != null) {
			initialize(model);
			messageCacheSparse = (SparseMessageCache)messageCache;
			labelSet = new BitSet(model.getMaxArity());
		}
	}
	public SparseScalableMP(GraphicalModel model, int debug) {
		this(model);
		this.debugLvl = debug;
	}
	public SparseScalableMP(GraphicalModel model, int debug, Method method, boolean maxProd, BitSet activeNodes) {
		this(method, maxProd);
		if (model != null) {
			initialize(model);
			messageCacheSparse = (SparseMessageCache)messageCache;
			labelSet = new BitSet(model.getMaxArity());
		}
		this.debugLvl = debug;
		this.observedNodes = activeNodes;
	}
	public SparseScalableMP(Method method,  boolean maxProd) {
		super(method, maxProd);
	}
	public static class SparseMessageCache extends ScalableMessagePassing.MessageCache {
		private static final long serialVersionUID = 1L;
		/*int offsets[][];
		int currentEdge = 0; int currentDir = 0, currentOffset = 0;*/
		TIntArrayList labels[];
		TDoubleArrayList messages[];
		int prevOffset=0;
		double defMessages[][];
		int initialCap;
		
		public SparseMessageCache(){}
		
		public SparseMessageCache(GraphicalModel graph, FwdBwdGraph orderedGraph) {
			initialCap = orderedGraph.numEdges() * 10;
			labels = new TIntArrayList[2];
			labels[1] = new TIntArrayList(initialCap);
			labels[0] = new TIntArrayList(initialCap);
			messages = new TDoubleArrayList[2];
			messages[1] = new TDoubleArrayList(initialCap);
			messages[0] = new TDoubleArrayList(initialCap);
			offsets = new int[2][graph.getGraph().getNumEdges()];
			defMessages = new double[2][graph.getGraph().getNumEdges()];		
			currentOffset=-1;
			prevOffset=-1;
		}
		@Override
		public void clear() {	// It is not needed in this algorithm
			throw new NotImplementedException();
		}
		/* Message Iterator */
		public class MessageIterator implements Serializable {
			private static final long serialVersionUID = 1L;
			
			int startIndex;
			int endIndex;
			int cdir;
			public MessageIterator() {}
			public MessageIterator(int src, int dest, int edgeNum) {
				init(src, dest, edgeNum);
			}
			public MessageIterator init(int src, int dest, int edgeNum) {
				this.cdir = src < dest ? 1 : 0;
				this.startIndex = getOffset(edgeNum, cdir);
				this.endIndex = getOffset(edgeNum+1, cdir);
				return this;
			}
			public int nextLabelValue(double[] val) {
				val[0] = 0;
				if(startIndex < 0) return -1;
				if(startIndex < endIndex){
					val[0] = messages[cdir].get(startIndex);
					int lbl = labels[cdir].get(startIndex);
					++startIndex;
					return lbl;
				}
				else
					return -1;
			}
		}
		public MessageIterator getMessageIterator(int src, int dest, int edgeNum, MessageIterator msgIter) {
			if (msgIter == null) {
				msgIter = new MessageIterator(src, dest, edgeNum);
			} else {
				msgIter.init(src, dest, edgeNum);
			}
			return msgIter;
		}		
		@Override
		public double get(int sender, int recvr, int edgeNum, int recvLabel) {
			//throw new NotImplementedException();
			MessageIterator msgIter = null;
			double[] score = new double[1];
			int lbl = -1;
			for(msgIter = getMessageIterator(sender, recvr, edgeNum, msgIter),lbl=msgIter.nextLabelValue(score);lbl>=0;lbl = msgIter.nextLabelValue(score)){
				if(lbl == recvLabel)
					return score[0];
			}
			return getDefMessage(sender, recvr, edgeNum);
		}
		private int getOffset(int edgeNum, int dir) {
			// When we want messageSize of last edgeNum, we return end index upto which labels array is filled
			return (edgeNum < offsets[dir].length ? (offsets[dir][edgeNum]==-1?currentOffset+1:offsets[dir][edgeNum]) : (dir==currentDir?currentOffset+1:prevOffset));	
		}
		@Override
		public void startPass(int dir) {
			prevOffset=currentOffset+1;
			currentEdge = -1;
			currentDir = dir;
			currentOffset=-1;
			messages[dir].reset();
			labels[dir].reset();
			Arrays.fill(offsets[dir], -1);
			Arrays.fill(defMessages[dir], 0.0);
		}
		@Override
		public void nextEdge() {
			currentEdge++;
			offsets[currentDir][currentEdge] = currentOffset+1;
		}
		@Override
		public void addMessage(int recvLabel, double value) {
			currentOffset += 1;
			labels[currentDir].add(recvLabel);
			messages[currentDir].add(value);
		}
		public void addDefMessage(int sender, int recvr, int edgeNum, double defVal){
			defMessages[sender<recvr?1:0][edgeNum] = defVal;
		}	
		public double getDefMessage(int sender, int recvr, int edgeNum){
			return defMessages[sender<recvr?1:0][edgeNum];
		}
	}	
	@Override
	protected double computeAndSubstractMin(LabelScoreArray di, int i) {
		LabelScoreArraySparse diS=(LabelScoreArraySparse)di;
		int activeArity = diS.getNumActive();
		if(activeArity==0) {
			double copyVal = defaultValue;
			defaultValue = 0;
			return copyVal;
		}
		double vMin = Double.POSITIVE_INFINITY;
		if(activeArity < model.getMaxArity()) 
			vMin = defaultValue;
		for (int k=0; k<activeArity; k++){
			vMin = Math.min(diS.get(diS.getActiveLabel(k)), vMin); 
		}
		for (int k=0; k<activeArity; k++) {
			diS.add(diS.getActiveLabel(k), -vMin);
		}
		defaultValue -= vMin;
		return vMin;
	}
	@Override
	protected double computeMin(LabelScoreArray di, int nodeId) {
		LabelScoreArraySparse diS=(LabelScoreArraySparse)di;
		int activeArity = diS.getNumActive();
		int kMin = 0;
		double vMin = Integer.MAX_VALUE;		
		for (int k=0; k<activeArity; ++k){
			int l = diS.getActiveLabel(k);
			if (vMin > diS.get(l)){
				vMin = diS.get(l);
				kMin = l;
			} 
			else if (vMin == diS.get(l)) {
				if(kMin > -1) kMin=Math.min(kMin, l);
				//System.out.println("Tie at label "+k+ " and "+kMin);
			}
		}
		solutionLabel[nodeId] = (activeArity==0) ? 0 : kMin;
		return (activeArity==0) ? 0 : vMin;
	}	

	@Override
	protected void computeMarginalSolution(){
		LabelScoreArraySparse diS = (LabelScoreArraySparse)tempBufSrc;
		int numNodes = model.getGraph().getNumNodes();
		for (int nodeId=0; nodeId < numNodes; nodeId++) {
			diS.clear();
			activateNodeLabels(diS, nodeId);
			copyNodePotential(diS, nodeId);
			addMessagesToNode(nodeId, diS);
			double normalizeConst = RobustMath.LOG0;
			normalizeConst = RobustMath.logSumExp(normalizeConst, defaultValue + Math.log(model.getNodeArity(nodeId)-diS.getNumActive()));
			double logConf = defaultValue;
			int predLabel = -1;
			if(diS.getNumActive() > 0){
				for(int k=diS.getNumActive() - 1; k>=0; --k){
					int lbl = diS.getActiveLabel(k);
					double logPhi = diS.get(lbl);
					normalizeConst = RobustMath.logSumExp(normalizeConst, logPhi);
					if(logConf < logPhi){
						logConf = logPhi;
						predLabel = lbl;
					}
				}
				assert(logConf >= defaultValue);
				if(logConf < defaultValue){
					if(debugLvl>0) System.out.println("Assigned label coming out to be Inactive for "+nodeId+" with "+diS.getNumActive()+" active labels and confidence is "+logConf+"/"+normalizeConst);
					labelSet.clear();
					for(int k=diS.getNumActive() - 1; k>=0; --k){
						labelSet.set(diS.getActiveLabel(k));
					}
					predLabel = labelSet.nextClearBit(0);
					logConf = defaultValue;
				}
			}
			confidence[nodeId] = Math.exp(logConf - normalizeConst);
			solutionLabel[nodeId] = predLabel;
			//System.out.println(nodeId + " "+ predLabel+" "+confidence[nodeId]);
		}
	}
	@Override
	protected void copyList(LabelScoreArray dst, LabelScoreArray src){
		LabelScoreArraySparse dstS = (LabelScoreArraySparse)dst;
		LabelScoreArraySparse srcS = (LabelScoreArraySparse)src;
		dstS.clear();
		dstS.add(srcS);
	}	
	@Override
	protected double updateMessages(LabelScoreArray di, int edgeNum, int srcNode, int dstNode) {
		LabelScoreArraySparse diS = ( LabelScoreArraySparse)di;
		tempBufDstS = (LabelScoreArraySparse) tempBufDst;
		tempBufSrcS = (LabelScoreArraySparse) tempBufSrc;
		tempBufDstS.clear(); 
		tempBufSrcS.clear();
		
		double gamma = getGamma(edgeNum, srcNode, dstNode);
		int reverseEdgeNum = orderedGraph.getEdgeNumber(dstNode, srcNode);

		for(int j=diS.getNumActive()-1; j>=0; --j){
			tempBufSrcS.add(diS.getActiveLabel(j), 0);	// Active Labels for theta^hat
		}
		/*if (method != Method.MeanField) {
			int lbl=-1;
			for(msgIter = messageCacheSparse.getMessageIterator(dstNode, srcNode, reverseEdgeNum, msgIter), lbl = msgIter.nextLabelValue(score); lbl >=0; lbl = msgIter.nextLabelValue(score)){
				tempBufSrcS.add(lbl, 0);
			}
		}*/
		double defMsg = messageCacheSparse.getDefMessage(dstNode, srcNode, reverseEdgeNum);
		int arity = tempBufSrcS.getNumActive();
		for(int i=0; i<arity; i++){
			int label = tempBufSrcS.getActiveLabel(i);
			if(diS.isLabelActive(label))
				tempBufSrcS.add(label, gamma*diS.get(label));
			else
				tempBufSrcS.add(label, gamma*defaultValue);
			if (method != Method.MeanField)
				tempBufSrcS.add(label, -defMsg);
		}
		double inactiveLabelValue = gamma*defaultValue;
		if (method != Method.MeanField){
			int lbl = -1;
			for(msgIter = messageCacheSparse.getMessageIterator(dstNode, srcNode, reverseEdgeNum, msgIter), lbl = msgIter.nextLabelValue(score); lbl >=0; lbl = msgIter.nextLabelValue(score)){
				tempBufSrcS.add(lbl, -score[0] + defMsg);
			}
			inactiveLabelValue -= defMsg;
		}
		
		//gsc: Assume Edge potentials are Potts potential
		double vMin = Double.POSITIVE_INFINITY;
		defMsg = 0;
		if(maxProd){
			defMsg = marginalComputer.initValue();
			for (int ksource=0; ksource<arity; ++ksource) {
				double v =  tempBufSrcS.get(tempBufSrcS.getActiveLabel(ksource));
				// v =v + \theta{ij}(x,y); Assumed \theta{ij}(x,y) = 0; default message is for inactive labels, hence edge potential is zero for active labels.
				defMsg = marginalComputer.aggregate(defMsg, v);
			}
			double preComputedVal = defMsg;	
			defMsg = marginalComputer.aggregate(defMsg, inactiveLabelValue + (method==Method.TRWS ? -1 : 1) * model.getEdgePotential(srcNode, 0, dstNode, 0));
			vMin = Math.min(vMin, defMsg);
			
			preComputedVal = marginalComputer.aggregate(preComputedVal, inactiveLabelValue);	
			for (int kdest=0; kdest<arity; ++kdest) {
				int kdestLbl = tempBufSrcS.getActiveLabel(kdest);
				double pot = marginalComputer.aggregate(preComputedVal, tempBufSrcS.get(kdestLbl) + (method==Method.TRWS ? -1 : 1) * model.getEdgePotential(srcNode, kdestLbl, dstNode, kdestLbl));
				tempBufDstS.add(kdestLbl, pot);
				vMin = Math.min(vMin, pot);
			}
		}
		else{
			double preComputedVal = marginalComputer.initValue();
			for (int ksource=0; ksource<arity; ksource++) {
				double v =  tempBufSrcS.get(tempBufSrcS.getActiveLabel(ksource));
				preComputedVal = marginalComputer.aggregate(preComputedVal, v);
			}
			//double activePreComputedVal = preComputedVal;
			int numInActiveLabels = model.getMaxArity() - arity;
			defMsg = 0;
			if(numInActiveLabels > 0){
				preComputedVal = marginalComputer.aggregate(preComputedVal, Math.log(numInActiveLabels) + inactiveLabelValue);
				try {
					//defMsg = marginalComputer.aggregate(activePreComputedVal, Math.log(numInActiveLabels-1) + inactiveLabelValue);
					defMsg = RobustMath.logMinusExp(preComputedVal, inactiveLabelValue);
					defMsg = marginalComputer.aggregate(defMsg, inactiveLabelValue +(method==Method.TRWS ? -1 : 1) * model.getEdgePotential(srcNode, 0, dstNode, 0));
				} catch (Exception e) { e.printStackTrace(); System.exit(1);}
				vMin = Math.min(vMin, defMsg);
			}
			for (int kdest=0; kdest<arity; ++kdest) {
				int kdestLbl = tempBufSrcS.getActiveLabel(kdest);
				double temp = tempBufSrcS.get(kdestLbl);
				double v = 0.0;
				try {
					v = RobustMath.logMinusExp(preComputedVal, temp);
				} catch (Exception e) { e.printStackTrace(); }
				double pot = marginalComputer.aggregate(v, temp + (method==Method.TRWS ? -1 : 1) * model.getEdgePotential(srcNode, kdestLbl, dstNode, kdestLbl));
				tempBufDstS.add(kdestLbl, pot);
				vMin = Math.min(vMin, pot);
			}
			//vMin = Math.min(vMin, defMsg);
		}	
		//tempBufDstS.sortTopKByScore(K);
		if(defMsg - vMin < 0) defMsg = vMin;
		messageCacheSparse.addDefMessage(srcNode, dstNode, edgeNum, defMsg - vMin);
		for (int kdest=0; kdest<arity; ++kdest){
			int kdestLbl = tempBufDstS.getActiveLabel(kdest);
			double pot = tempBufDstS.get(kdestLbl);
			if(!toBePruned)
				messageCacheSparse.addMessage(kdestLbl, pot-vMin);
			else if(Math.abs(pot-defMsg) > THRESHOLD)
				 messageCacheSparse.addMessage(kdestLbl, pot-vMin);
		}
		
		return (arity>0) ? vMin : defMsg;
	}	
	
	/*protected double logExpSubtract(double v1, double v2){
		assert(v1 >= v2):v1+" "+v2;	// Code expects v1 to be greater than or equal to v2.
		if (Math.abs(v1-v2) < Double.MIN_VALUE)
	        return Double.NEGATIVE_INFINITY;
	    double vmin = Math.min(v1,v2);
	    double vmax = Math.max(v1,v2);
	    if ( vmax > vmin + 30 ) {
	        return vmax;
	    }
	    else{	
	    	double retval = vmax + Math.log(1.0 - Math.exp(vmin-vmax));
	    	return retval;
	    }
	}*/
	
	@Override
	protected void activateNodeLabels(LabelScoreArray di, int nodeId) {
		LabelScoreArraySparse diS = (LabelScoreArraySparse)di;
		nPotS = (NodePotentialSparse)model.getNodePotentialTable(nodeId);
		for(int j=nPotS.getNodeArity()-1; j>=0; --j){
			diS.add(nPotS.getLabel(j), 0);
		}
		int nbr = -1;
		for (nbrIter = model.getGraph().getNeighbourIterator(nodeId, nbrIter), nbr = nbrIter.nextNbr(); nbr != -1; nbr = nbrIter.nextNbr()) {
			int k=-1;
			for(msgIter = messageCacheSparse.getMessageIterator(nbr, nodeId, orderedGraph.getEdgeNumber(nbr, nodeId), msgIter), k=msgIter.nextLabelValue(score); k>=0; k=msgIter.nextLabelValue(score)){ 
				diS.add(k, 0);
			}
		}
	}	
	@Override
	protected void addEdgeMessages(LabelScoreArray di, int edgeNum, int dst, int src) {
		LabelScoreArraySparse diS=(LabelScoreArraySparse)di;
		double defMsg = messageCacheSparse.getDefMessage(src, dst, edgeNum);
		for(int j=diS.getNumActive()-1; j>=0; --j){
			diS.add(diS.getActiveLabel(j), defMsg);
		}
		int lbl = -1;
		for(msgIter = messageCacheSparse.getMessageIterator(src, dst, edgeNum, msgIter), lbl=msgIter.nextLabelValue(score); lbl>=0; lbl=msgIter.nextLabelValue(score)){
			diS.add(lbl, score[0] - defMsg);
		}
		defaultValue += (method==Method.TRWS ? -1 : 1) * defMsg;
	}
	@Override
	protected void findActiveLabels(LabelScoreArray di, int nodeId) {
		LabelScoreArraySparse diS=(LabelScoreArraySparse)di;
		nPotS = (NodePotentialSparse)model.getNodePotentialTable(nodeId);
		for(int j=nPotS.getNodeArity()-1; j>=0; --j){
			diS.add(nPotS.getLabel(j), 0);
		}
		EdgeIter edgeIter = null;
		int nbr = -1;
		for (edgeIter=orderedGraph.getEdgeIter(nodeId, 0, edgeIter), nbr = edgeIter.nextNbr(); nbr >= 0; nbr=edgeIter.nextNbr()) {
			diS.add(solutionLabel[nbr], 0);
		}
		edgeIter = null;
		nbr = -1;
		for (edgeIter=orderedGraph.getEdgeIter(nodeId, 1, edgeIter), nbr = edgeIter.nextNbr(); nbr >= 0; nbr=edgeIter.nextNbr()){
			int lbl=-1;
			for(msgIter = messageCacheSparse.getMessageIterator(nbr, nodeId, orderedGraph.getEdgeNumber(nbr, nodeId), msgIter), lbl = msgIter.nextLabelValue(score); lbl >=0; lbl = msgIter.nextLabelValue(score)){
				diS.add(lbl, 0);
			}
		}
	}

	protected double defaultValue = 0;
	@Override
	protected void copyNodePotential(LabelScoreArray di, int nodeId) {
		LabelScoreArraySparse diS=(LabelScoreArraySparse)di;
		nPotS = (NodePotentialSparse)model.getNodePotentialTable(nodeId);
		for(int j=diS.getNumActive()-1; j>=0; j--){
			int lbl = diS.getActiveLabel(j);
			diS.add(lbl, (method==Method.TRWS ? -1 : 1) * nPotS.getPotentialValue(lbl));	// getPotentialValue will return defPot if lbl is not active
		}
		defaultValue = (method==Method.TRWS ? -1 : 1) * nPotS.getDefaultPotentialValue();
	}

	@Override
	protected void addColumnToMessage(LabelScoreArray di, int source, int dest) {
		LabelScoreArraySparse diS=(LabelScoreArraySparse)di;
		int sourcelabel = solutionLabel[source];
		for(int k = diS.getNumActive()-1; k>=0; k--){
			int lbl = diS.getActiveLabel(k);
			diS.add(lbl, -model.getEdgePotential(source, sourcelabel, dest, lbl));
		}
	}
	
	@Override
	protected void addMessage(LabelScoreArray di, int dstNode, int srcNode, int edgeNum) {
		LabelScoreArraySparse diS = (LabelScoreArraySparse)di;
		double defMsg = messageCacheSparse.getDefMessage(srcNode, dstNode, edgeNum);
		for (int k=diS.getNumActive()-1; k>=0; k--) {
			diS.add(diS.getActiveLabel(k), defMsg);
		}
		int lbl = -1;
		for(msgIter = messageCacheSparse.getMessageIterator(srcNode, dstNode, edgeNum, msgIter), lbl = msgIter.nextLabelValue(score); lbl>=0; lbl=msgIter.nextLabelValue(score)){
			diS.add(lbl, score[0] - defMsg);
		}
	}
	
	@Override
	protected LabelScoreArray newLabelScoreArray(int sz) {
		return new LabelScoreArraySparse(sz);
	}
	
	@Override
	protected MessageCache newMessageCache() {
		return new SparseMessageCache(model, orderedGraph);
	}
	
	/**
	 * Get sum-marginal probability distribution table for nodeId. Also gives activeLabels 
	 * in the marginal[].
	 */
	public double getNodeSumMarginal(int nodeId, double[] marginal, BitSet activeLabels){
		LabelScoreArraySparse diS = (LabelScoreArraySparse)tempBufSrc;
		diS.clear();
		activateNodeLabels(diS, nodeId);
		copyNodePotential(diS, nodeId);
		addMessagesToNode(nodeId, diS);
		Arrays.fill(marginal, defaultValue);
		double normalizeConst = RobustMath.LOG0;
		normalizeConst = RobustMath.logSumExp(normalizeConst, defaultValue + Math.log(model.getNodeArity(nodeId)-diS.getNumActive()));
		if(diS.getNumActive() > 0){
			for(int k=diS.getNumActive() - 1; k>=0; --k){
				int lbl = diS.getActiveLabel(k);
				marginal[lbl] = diS.get(lbl);
				if(activeLabels!=null) activeLabels.set(lbl);
				normalizeConst = RobustMath.logSumExp(normalizeConst, marginal[lbl]);
			}
		}
		double defVal = 0;
		for(int i=0; i<marginal.length; ++i){
			marginal[i] = Math.exp(marginal[i] - normalizeConst);
			if(!diS.isLabelActive(i)){defVal = marginal[i];}
		}
		return defVal;
	}
}
