package iitb.sgl.inference;

import java.io.Serializable;
import java.util.BitSet;
import java.util.Properties;

import gnu.trove.list.array.TIntArrayList;
import iitb.shared.RobustMath;
import iitb.shared.gm.NodePotentialSparse;
import iitb.shared.gm.UDGM;
import iitb.shared.gm.inference.LabelScoreArraySparse;
import iitb.shared.gm.inference.SolutionWithBounds;
import iitb.shared.gm.inference.SparseScalableMP;

/**
 * Utilities supporting SparseScalableMP such as Entropy based logZ computation etc.
 * During learning, we use SocInference (extending SparseScalableMP) as a subroutine.
 * @author gaurish
 */

public class SGInference extends SparseScalableMP implements Serializable {
	private static final long serialVersionUID = 1L;
	
	LabelScoreArraySparse lambdaDi;
	public SGInference(UDGM model, int debugLvl, Method method, boolean maxProd,BitSet observedNodes) {
		super(model, debugLvl, method, maxProd, observedNodes);
		lambdaDi = (LabelScoreArraySparse)newLabelScoreArray(model.getMaxArity());
	}
	
	boolean isDegreeInvGamma = false;
	int maxIters = 100;
	
	public void setIsDegreeInvGamma(boolean val){
		this.isDegreeInvGamma = val;
	}
	public void setMaxIters(int iter){
		this.maxIters = iter;
	}
	
	@Override
	public SolutionWithBounds getTopK(int k) throws UnsupportedOperationException {
		Properties opt = new Properties();
		opt.setProperty("maxIters", ""+maxIters);
		run(opt, false);
		return getSolution();
	}
	
	@Override
	protected double getGamma(int edgeNum, int srcNode, int dstNode) {
		return (method==Method.TRWS ? 1.0/Math.max(orderedGraph.getFwdDegree(srcNode),orderedGraph.getBwdDegree(srcNode)) : (isDegreeInvGamma ? (1.0/(1 + model.getGraph().getNumNeighbours(srcNode))) : 1) );
	}
	
	public double getEdgeLambda(LabelScoreArraySparse diS, int srcNode, int dstNode){
		diS.clear(); 
		lambdaDi.clear();
		activateNodeLabels(lambdaDi, srcNode);
		copyNodePotential(lambdaDi, srcNode);
		addMessagesToNode(srcNode, lambdaDi);
		double gamma = 1.0 / (1 + model.getGraph().getNumNeighbours(srcNode)) ;
		int reverseEdgeNum = orderedGraph.getEdgeNumber(dstNode, srcNode);
		double defMsg = messageCacheSparse.getDefMessage(dstNode, srcNode, reverseEdgeNum);
		for(int j = lambdaDi.getNumActive()-1; j >= 0; --j){
			int lbl = lambdaDi.getActiveLabel(j);
			diS.add(lbl, (gamma * lambdaDi.get(lbl)) - defMsg);
		}
		int lbl = -1;
		for(msgIter = messageCacheSparse.getMessageIterator(dstNode, srcNode, reverseEdgeNum, msgIter), lbl = msgIter.nextLabelValue(score); lbl >=0; lbl = msgIter.nextLabelValue(score)){
			diS.add(lbl, -score[0] + defMsg);
		}	
		return (gamma * defaultValue - defMsg);
	}
	
	/**
	 * Get logZ for Trees with node potentials
	 */
	public double getTreeLogZ() throws Exception{
		int  arity = model.getMaxArity();
		double[] marginal = new double[arity];
		double logZ = 0;
		for(int node = model.getGraph().getNumNodes()-1; node>=0; --node){
			getNodeSumMarginal(node, marginal, null);
			logZ += (1 - model.getGraph().getNumNeighbours(node)) * entropy(marginal);
			NodePotentialSparse nPotS = (NodePotentialSparse)model.getNodePotentialTable(node);
			for(int m = nPotS.getNodeArity()-1 ; m>=0; --m)
				logZ += model.getNodePotential(node, nPotS.getLabel(m)) * marginal[nPotS.getLabel(m)];
			for(int j=model.getGraph().getNumNeighbours(node)-1; j>=0; --j){
				int nbr = model.getGraph().getNeighbour(node, j);
				if(node > nbr) continue;
				EdgeMarginal em = new EdgeMarginal(node, nbr);
				/*for(int x = model.getNodeArity(node)-1; x>=0; --x){
					for(int y = model.getNodeArity(nbr)-1; y>=0; --y){
						double p = em.getEdgeMarginalScore(x, y);
						if(p>0) logZ += - p * Math.log(p);
						logZ += p * model.getEdgePotential(node, x, nbr, y);
					}
				}*/
				for(int x = em.activeLabels.nextSetBit(0); x>=0; x = em.activeLabels.nextSetBit(x+1)){
					for(int y = 0; y < arity; ++y){
						double p = em.getEdgeMarginalScore(x, y);
						if(p > 0) logZ += - p * Math.log(p);
						logZ += p * model.getEdgePotential(node, x, nbr, y);
					}
				}
				for(int x = 0; x < arity; ++x){
					if(em.activeLabels.get(x)) continue;
					for(int y = em.activeLabels.nextSetBit(0); y>=0; y = em.activeLabels.nextSetBit(y+1)){
						double p = em.getEdgeMarginalScore(x, y);
						if(p > 0) logZ += - p * Math.log(p);
						logZ += p * model.getEdgePotential(node, x, nbr, y);
					}
				}
				int remArity = arity - em.activeLabels.cardinality();
				double p = em.getSameEdgeMarginalScore();
				if(p > 0) logZ += -remArity * p * Math.log(p);
				logZ += remArity * p * model.getEdgePotential(node, 0, nbr, 0);
				p = em.getDiffEdgeMarginalScore();
				if(p > 0) logZ += -(remArity*remArity - remArity) * p * Math.log(p);
				logZ += (remArity*remArity-remArity) * p * model.getEdgePotential(node, 0, nbr, 1);
			}
			// node processed
		}
		return logZ;
	}
	
	/**
	 * Get logZ for Trees without node potentials
	 */
	public double getNoObsTreeLogZ() throws Exception{
		int  k = model.getMaxArity();
		double[] marginal = new double[k];
		double logZ = 0;
		for(int node = model.getGraph().getNumNodes()-1; node>=0; --node){
			getNodeSumMarginal(node, marginal, null);
			logZ += (1 - model.getGraph().getNumNeighbours(node)) * entropy(marginal);
			for(int i=0; i<marginal.length; ++i){
				logZ += model.getNodePotential(node, i) * marginal[i];
			}
			for(int j=model.getGraph().getNumNeighbours(node)-1; j>=0; --j){
				int nbr = model.getGraph().getNeighbour(node, j);
				if(node < nbr){
					double mu = getEdgeSameLabelProbability(node, nbr);
					if(mu!=0) logZ += - mu * (Math.log(mu) - Math.log(k));
					if(mu!=1) logZ += -(1-mu) * (Math.log(1-mu) - Math.log(k*k-k));
					logZ += mu * model.getEdgePotential(node, 0, nbr, 0);
				}
			}
		}
		return logZ;
	}
	
	/**
	 * Returns Pr(x_i = x_j) for edge i---j
	 */
	public double getEdgeSameLabelProbability(int src, int dest) throws Exception{
		int arity = model.getNodeArity(src);
		assert(arity==model.getNodeArity(dest));
		tempBufSrcS = (LabelScoreArraySparse)tempBufSrc;
		tempBufDstS = (LabelScoreArraySparse)tempBufDst;
		// Compute Incoming Messages to SRC
		tempBufSrcS.clear();
		nPotS = (NodePotentialSparse)model.getNodePotentialTable(src);
		for(int j=nPotS.getNodeArity()-1; j>=0; --j)
			tempBufSrcS.add(nPotS.getLabel(j), 0);
		int nbr = -1;
		for (nbrIter = model.getGraph().getNeighbourIterator(src, nbrIter), nbr = nbrIter.nextNbr(); nbr != -1; nbr = nbrIter.nextNbr()) {
			if(nbr==dest) continue;
			int k = -1;
			for(msgIter = messageCacheSparse.getMessageIterator(nbr, src, orderedGraph.getEdgeNumber(nbr, src), msgIter), k=msgIter.nextLabelValue(score); k>=0; k=msgIter.nextLabelValue(score)) 
				tempBufSrcS.add(k, 0);
		}
		copyNodePotential(tempBufSrcS, src);
		nbr = -1;
		for (nbrIter = model.getGraph().getNeighbourIterator(src, nbrIter), nbr = nbrIter.nextNbr(); nbr != -1; nbr = nbrIter.nextNbr()) {
			if(nbr == dest) continue;
			addEdgeMessages(tempBufSrcS, orderedGraph.getEdgeNumber(nbr, src), src, nbr);
		}
		double defSrc = defaultValue;
		
		// Compute Incoming Messages to DEST
		tempBufDstS.clear();
		nPotS = (NodePotentialSparse)model.getNodePotentialTable(dest);
		for(int j=nPotS.getNodeArity()-1; j>=0; --j)
			tempBufDstS.add(nPotS.getLabel(j), 0);
		nbr = -1;
		for (nbrIter = model.getGraph().getNeighbourIterator(dest, nbrIter), nbr = nbrIter.nextNbr(); nbr != -1; nbr = nbrIter.nextNbr()) {
			if(nbr == src) continue;
			int k = -1;
			for(msgIter = messageCacheSparse.getMessageIterator(nbr, dest, orderedGraph.getEdgeNumber(nbr, dest), msgIter), k=msgIter.nextLabelValue(score); k>=0; k=msgIter.nextLabelValue(score)) 
				tempBufDstS.add(k, 0);
		}
		copyNodePotential(tempBufDstS, dest);
		nbr = -1;
		for (nbrIter = model.getGraph().getNeighbourIterator(dest, nbrIter), nbr = nbrIter.nextNbr(); nbr != -1; nbr = nbrIter.nextNbr()) {
			if(nbr == src) continue;
			addEdgeMessages(tempBufDstS, orderedGraph.getEdgeNumber(nbr, dest), dest, nbr);
		}
		double defDst = defaultValue;
		
		// COMPUTATIONS
		double potts = model.getEdgePotential(src, 0, dest, 0);
		double logNr = RobustMath.LOG0;
		double logTerm1 = RobustMath.LOG0;
		double logTerm2 = RobustMath.LOG0;
		double logRepeated =  RobustMath.LOG0;
		
		BitSet activeLabels = new BitSet();
		for(int i=tempBufSrcS.getNumActive()-1; i>=0; --i)
			activeLabels.set(tempBufSrcS.getActiveLabel(i));
		for(int i=tempBufDstS.getNumActive()-1; i>=0; --i)
			activeLabels.set(tempBufDstS.getActiveLabel(i));
		
		for(int x=activeLabels.nextSetBit(0); x>=0; x=activeLabels.nextSetBit(x+1)){
			logNr = RobustMath.logSumExp(logNr, potts + tempBufSrcS.get(x) + tempBufDstS.get(x));
			logRepeated = RobustMath.logSumExp(logRepeated, tempBufSrcS.get(x) + tempBufDstS.get(x));
		}
		logNr = RobustMath.logSumExp(logNr, Math.log(arity-activeLabels.cardinality()) + potts + defSrc + defDst);
		logRepeated = RobustMath.logSumExp(logRepeated, Math.log(arity-activeLabels.cardinality()) + defSrc + defDst);

		for(int i=tempBufSrcS.getNumActive()-1; i>=0; --i)
			logTerm1 = RobustMath.logSumExp(logTerm1, tempBufSrcS.get(tempBufSrcS.getActiveLabel(i)));
		logTerm1 = RobustMath.logSumExp(logTerm1, Math.log(arity-tempBufSrcS.getNumActive()) + defSrc);
		for(int i=tempBufDstS.getNumActive()-1; i>=0; --i)
			logTerm2 = RobustMath.logSumExp(logTerm2, tempBufDstS.get(tempBufDstS.getActiveLabel(i)));
		logTerm2 = RobustMath.logSumExp(logTerm2, Math.log(arity-tempBufDstS.getNumActive()) + defDst);
		
		double logTerm = logTerm1 + logTerm2;
		logTerm = RobustMath.logMinusExp(logTerm, logRepeated);// + Math.log(arity*(arity-1));
		double logDn = RobustMath.logSumExp(logNr, logTerm);
		
		return Math.exp(logNr-logDn);
	}
	
	/**
	 * It returns un-normalized partial belief at nodeId (i.e. product of incoming messages coming 
	 * from all neighbors except those in minusNbrs list)
	 */
	public double getPartialBelief(int nodeId, TIntArrayList minusNbrs, LabelScoreArraySparse belief){
		belief.clear();
		// Get Active Labels for Belief
		nPotS = (NodePotentialSparse)model.getNodePotentialTable(nodeId);
		for(int j = nPotS.getNodeArity()-1; j>=0; --j)
			belief.add(nPotS.getLabel(j), 0);
		int nbr = -1;
		for (nbrIter = model.getGraph().getNeighbourIterator(nodeId, nbrIter), nbr = nbrIter.nextNbr(); nbr != -1; nbr = nbrIter.nextNbr()) {
			if(!minusNbrs.contains(nbr)) {
				int lbl = -1;
				for(msgIter = messageCacheSparse.getMessageIterator(nbr, nodeId, orderedGraph.getEdgeNumber(nbr, nodeId), msgIter), lbl=msgIter.nextLabelValue(score); lbl>=0; lbl=msgIter.nextLabelValue(score))
					belief.add(lbl, 0);
			}
		}
		// Compute Belief
		copyNodePotential(belief, nodeId);
		nbr = -1;
		for (nbrIter = model.getGraph().getNeighbourIterator(nodeId, nbrIter), nbr = nbrIter.nextNbr(); nbr != -1; nbr = nbrIter.nextNbr()) {
			if(!minusNbrs.contains(nbr))
				addEdgeMessages(belief, orderedGraph.getEdgeNumber(nbr, nodeId), nodeId, nbr);
		}
		return defaultValue;
	}
	
	class EdgeMarginal{
		int src;
		int dest;
		double normalizer;
		LabelScoreArraySparse beliefSrc;
		LabelScoreArraySparse beliefDst;
		BitSet activeLabels;
		double defSrc;
		double defDst;
		
		public EdgeMarginal(int src, int dest) throws Exception{
			this.src = src;
			this.dest = dest;
			this.activeLabels = new BitSet(model.getNodeArity(src));
			this.normalizer = getNormalizer(src, dest); 
		}
		
		public double getEdgeMarginalScore(int label1, int label2){
			double p = model.getEdgePotential(src, label1, dest, label2) + beliefSrc.get(label1) + beliefDst.get(label2);
			return Math.exp(p - normalizer);
		}
		public double getSameEdgeMarginalScore(){
			double p = model.getEdgePotential(src, 0, dest, 0) + defSrc + defDst;
			return Math.exp(p - normalizer);
		}
		public double getDiffEdgeMarginalScore(){
			double p = model.getEdgePotential(src, 0, dest, 1) + defSrc + defDst;
			return Math.exp(p - normalizer);
		}
		
		public double getNormalizer(int src, int dest) throws Exception{
			int arity = model.getNodeArity(src); assert(arity==model.getNodeArity(dest));
			tempBufSrcS = (LabelScoreArraySparse)tempBufSrc;
			tempBufDstS = (LabelScoreArraySparse)tempBufDst;
			// Compute Incoming Messages to SRC
			tempBufSrcS.clear();
			nPotS = (NodePotentialSparse)model.getNodePotentialTable(src);
			for(int j=nPotS.getNodeArity()-1; j>=0; --j)
				tempBufSrcS.add(nPotS.getLabel(j), 0);
			int nbr = -1;
			for (nbrIter = model.getGraph().getNeighbourIterator(src, nbrIter), nbr = nbrIter.nextNbr(); nbr != -1; nbr = nbrIter.nextNbr()) {
				if(nbr==dest) continue;
				int k = -1;
				for(msgIter = messageCacheSparse.getMessageIterator(nbr, src, orderedGraph.getEdgeNumber(nbr, src), msgIter), k=msgIter.nextLabelValue(score); k>=0; k=msgIter.nextLabelValue(score)) 
					tempBufSrcS.add(k, 0);
			}
			copyNodePotential(tempBufSrcS, src);
			nbr = -1;
			for (nbrIter = model.getGraph().getNeighbourIterator(src, nbrIter), nbr = nbrIter.nextNbr(); nbr != -1; nbr = nbrIter.nextNbr()) {
				if(nbr == dest) continue;
				addEdgeMessages(tempBufSrcS, orderedGraph.getEdgeNumber(nbr, src), src, nbr);
			}
			this.defSrc = defaultValue;
			
			// Compute Incoming Messages to DEST
			tempBufDstS.clear();
			nPotS = (NodePotentialSparse)model.getNodePotentialTable(dest);
			for(int j=nPotS.getNodeArity()-1; j>=0; --j)
				tempBufDstS.add(nPotS.getLabel(j), 0);
			nbr = -1;
			for (nbrIter = model.getGraph().getNeighbourIterator(dest, nbrIter), nbr = nbrIter.nextNbr(); nbr != -1; nbr = nbrIter.nextNbr()) {
				if(nbr == src) continue;
				int k = -1;
				for(msgIter = messageCacheSparse.getMessageIterator(nbr, dest, orderedGraph.getEdgeNumber(nbr, dest), msgIter), k=msgIter.nextLabelValue(score); k>=0; k=msgIter.nextLabelValue(score)) 
					tempBufDstS.add(k, 0);
			}
			copyNodePotential(tempBufDstS, dest);
			nbr = -1;
			for (nbrIter = model.getGraph().getNeighbourIterator(dest, nbrIter), nbr = nbrIter.nextNbr(); nbr != -1; nbr = nbrIter.nextNbr()) {
				if(nbr == src) continue;
				addEdgeMessages(tempBufDstS, orderedGraph.getEdgeNumber(nbr, dest), dest, nbr);
			}
			this.defDst = defaultValue;
			
			this.beliefSrc = tempBufSrcS;
			this.beliefDst = tempBufDstS;
			/*double logDn = RobustMath.LOG0;
			for(int x=0; x<arity; ++x){
				for(int y=0; y<arity; ++y)
					logDn = RobustMath.logSumExp(logDn, model.getEdgePotential(src, x, dest, y) + tempBufSrcS.get(x) + tempBufDstS.get(y));
			}*/
			double potts = model.getEdgePotential(src, 0, dest, 0);
			double logNr = RobustMath.LOG0;
			double logTerm1 = RobustMath.LOG0;
			double logTerm2 = RobustMath.LOG0;
			double logRepeated =  RobustMath.LOG0;
			
			this.activeLabels = new BitSet();
			for(int i=tempBufSrcS.getNumActive()-1; i>=0; --i)
				activeLabels.set(tempBufSrcS.getActiveLabel(i));
			for(int i=tempBufDstS.getNumActive()-1; i>=0; --i)
				activeLabels.set(tempBufDstS.getActiveLabel(i));
			
			for(int x=activeLabels.nextSetBit(0); x>=0; x=activeLabels.nextSetBit(x+1)){
				logNr = RobustMath.logSumExp(logNr, potts + tempBufSrcS.get(x) + tempBufDstS.get(x));
				logRepeated = RobustMath.logSumExp(logRepeated, tempBufSrcS.get(x) + tempBufDstS.get(x));
			}
			logNr = RobustMath.logSumExp(logNr, Math.log(arity-activeLabels.cardinality()) + potts + defSrc + defDst);
			logRepeated = RobustMath.logSumExp(logRepeated, Math.log(arity-activeLabels.cardinality()) + defSrc + defDst);

			for(int i=tempBufSrcS.getNumActive()-1; i>=0; --i)
				logTerm1 = RobustMath.logSumExp(logTerm1, tempBufSrcS.get(tempBufSrcS.getActiveLabel(i)));
			logTerm1 = RobustMath.logSumExp(logTerm1, Math.log(arity-tempBufSrcS.getNumActive()) + defSrc);
			for(int i=tempBufDstS.getNumActive()-1; i>=0; --i)
				logTerm2 = RobustMath.logSumExp(logTerm2, tempBufDstS.get(tempBufDstS.getActiveLabel(i)));
			logTerm2 = RobustMath.logSumExp(logTerm2, Math.log(arity-tempBufDstS.getNumActive()) + defDst);
			
			double logTerm = logTerm1 + logTerm2;
			logTerm = RobustMath.logMinusExp(logTerm, logRepeated);// + Math.log(arity*(arity-1));
			double logDn = RobustMath.logSumExp(logNr, logTerm);
			return logDn;
		}
	}
	
	private double entropy(double[] marginal){
		double h=0;
		for(int i=0; i<marginal.length; ++i){
			if(marginal[i]>0) h += - marginal[i] * Math.log(marginal[i]);
		}
		return h;
	}
}