package iitb.shared.gm.inference;

import iitb.shared.gm.inference.LabelScoreArray;
import iitb.shared.gm.inference.SolutionWithBounds;
import iitb.shared.gm.inference.ScalableMessagePassing.FwdBwdGraph.EdgeIter;
import iitb.shared.RobustMath;
import iitb.shared.graphs.UDGraph;
import iitb.shared.graphs.UDGraph.NbrIterator;
import iitb.shared.gm.GraphicalModel;
import iitb.shared.gm.PotentialInterface;

import java.io.Serializable;
import java.util.Arrays;
import java.util.BitSet;
import java.util.Properties;

import sun.reflect.generics.reflectiveObjects.NotImplementedException;

public class ScalableMessagePassing implements Serializable {
	private static final long serialVersionUID = 1L;

	public static class FwdBwdGraph implements Serializable {
		private static final long serialVersionUID = 1L;
		UDGraph graph;
		int fwdEdgeOffsets[];
		int backEdgeOffsets[];
		public FwdBwdGraph(UDGraph graph) {
			this.graph = graph;
			fwdEdgeOffsets = new int[graph.getNumNodes()+1];
			backEdgeOffsets = new int[graph.getNumNodes()];
			for (int i = 0; i < fwdEdgeOffsets.length-1; i++) {
				fwdEdgeOffsets[i+1] = fwdEdgeOffsets[i];
				int sz = graph.getNumNeighbours(i);
				for (int n = 0; n < sz; n++) {
					int nbr = graph.getNeighbour(i, n);
					if (nbr > i) fwdEdgeOffsets[i+1]++;
				}
			}
			for (int i = backEdgeOffsets.length-1; i > 0; i--) {
				backEdgeOffsets[i-1] = backEdgeOffsets[i];
				int sz = graph.getNumNeighbours(i);
				for (int n = 0; n < sz; n++) {
					int nbr = graph.getNeighbour(i, n);
					if (nbr < i) backEdgeOffsets[i-1]++;
				}
			}
		}
		public int numEdges() {return fwdEdgeOffsets[fwdEdgeOffsets.length-1];}
		public class EdgeIter implements NbrIterator {
			NbrIterator nbrIter;
			int edgeNum=-1;
			int srcNode;
			int forward;
			public EdgeIter(int node, int dir) {
				init(node,dir);
			}
			public EdgeIter init(int node,int dir) {
				this.forward=dir;
				this.srcNode = node;
				edgeNum = -1;
				nbrIter = graph.getNeighbourIterator(node, nbrIter);
				return this;
			}
			public int nextNbr() {
				do {
					int nbr = nbrIter.nextNbr();
					if (nbr < 0) return -1;
					if (valid(nbr)) {
						edgeNum++;
						return nbr;
					}
				} while (true);
			}
			protected boolean valid(int nbr) {
				return (forward==1?nbr > srcNode:nbr < srcNode);
			}
			public int edgeNumber() {
				return (forward==1?fwdEdgeOffsets[srcNode]:backEdgeOffsets[srcNode])+edgeNum;
			}
		}
		public EdgeIter getEdgeIter(int node, int dir, EdgeIter iter) {
			if (iter==null) {
				iter = new EdgeIter(node,dir);
			} else {
				iter.init(node,dir);
			}
			return iter;
		}
		public int getEdgeNumber(int srcNode, int destNode) {
			if (srcNode < destNode) {
				return fwdEdgeOffsets[srcNode] + getNbrNumber(srcNode,destNode,true);
			} else {
				return backEdgeOffsets[srcNode]  + getNbrNumber(srcNode,destNode,false);
			}
		}
		NbrIterator tmpIter;
		private int getNbrNumber(int srcNode, int destNode, boolean forward) {
			tmpIter = graph.getNeighbourIterator(srcNode, tmpIter);
			for (int nbr = tmpIter.nextNbr(), nbrNum = 0; nbr >= 0; nbr = tmpIter.nextNbr()) {
				if (nbr == destNode) 
					return nbrNum;
				if ((nbr < srcNode && !forward) ||  (nbr > srcNode && forward)) nbrNum++;
			}
			return -Integer.MAX_VALUE;
		}
		public int getFwdDegree(int srcNode) {
			return fwdEdgeOffsets[srcNode+1] - fwdEdgeOffsets[srcNode];
		}
		public int getBwdDegree(int srcNode) {
			return (srcNode > 0?backEdgeOffsets[srcNode-1]:graph.getNumEdges()) - backEdgeOffsets[srcNode]; 
		}
	}
	public static class MessageCache implements Serializable {
		private static final long serialVersionUID = 1L;
		double msgs[][]; 
		int offsets[][];

		int currentEdge = 0; int currentDir = 0, currentOffset = 0;
		
		// gsc: Added a default constructor to extend this class in future subclasses
		public MessageCache(){		
		}
		
		public MessageCache(GraphicalModel graph, FwdBwdGraph orderedGraph) {
			int numMsgsFwd = 0, numMsgsBwd=0;
			for (int i = graph.getGraph().getNumNodes()-1; i >= 0; i--) {
				numMsgsFwd += orderedGraph.getBwdDegree(i)*graph.getNodeArity(i);
				numMsgsBwd += orderedGraph.getFwdDegree(i)*graph.getNodeArity(i);
			}
			msgs = new double[2][];
			msgs[1] = new double[numMsgsFwd];
			msgs[0] = new double[numMsgsBwd];
			offsets = new int[2][graph.getGraph().getNumEdges()];
		}

		public void clear() {
			Arrays.fill(msgs[0], 0);
			Arrays.fill(msgs[1], 0);
		}

		public double get(int sender, int recvr, int edgeNum, int recvLabel) {
			int dir = sender<recvr?1:0;
			if (offsets[dir][edgeNum] < 0) return 0;
			int pos = getOffset(edgeNum, dir)+recvLabel;
			return msgs[dir][pos];
			//return (pos < getOffset(edgeNum+1, dir)?msgs[dir][pos]:0);
		}

		private int getOffset(int edgeNum, int dir) {
			return (edgeNum < offsets[dir].length?offsets[dir][edgeNum]:msgs[dir].length);
		}
		public void startPass(int dir) {
			currentEdge = -1;
			currentDir = dir;
			currentOffset=-1;
			Arrays.fill(msgs[dir],0);
			Arrays.fill(offsets[dir], -1);
		}
		public void nextEdge() {
			currentEdge++;
			offsets[currentDir][currentEdge] = currentOffset+1;
		}
		public void addMessage(int recvLabel, double value) {
			msgs[currentDir][getOffset(currentEdge,currentDir)+recvLabel] += value;
			currentOffset = Math.max(getOffset(currentEdge,currentDir)+recvLabel,currentOffset);
		}
	}
	MessageCache messageCache;
	protected FwdBwdGraph orderedGraph=null;
	protected int m_iterMax;
	protected double m_eps;
	int debugLvl=0;
	protected GraphicalModel model;
	int solutionLabel[];
	double energy;
	protected double lowerBoundPrev;
	protected double lowerBound;
	public static enum Method {BP,TRWS,MeanField};
	protected SolutionWithBounds solution;
	boolean maxProd=true;
	//gsc: array of winning marginal probabilities of all nodes.
	protected double confidence[];
	protected BitSet observedNodes;
	
	static class MarginalComputer implements Serializable{
		private static final long serialVersionUID = 1L;
		double initValue() {
			return Double.POSITIVE_INFINITY;
		}
		double aggregate(double aggValue, double elem) {
			return Math.min(aggValue, elem);
		}
	}
	static class SumMarginalComputer extends MarginalComputer {
		private static final long serialVersionUID = 1L;
		double initValue() {
			return RobustMath.LOG0;
		}
		double aggregate(double aggValue, double elem) {
			return RobustMath.logSumExp(aggValue, elem);
		}
	}
	MarginalComputer marginalComputer;
	protected Method method=Method.TRWS;
	/*
	 * To run BP call this with trws=false and maxProd=true or false depending on max product or sum product message passing.
	 */
	public ScalableMessagePassing(boolean trws, boolean maxProd){
		this(trws ? Method.TRWS : Method.BP, maxProd);
		
	}
	/*
	 * To run BP call this with trws=false and maxProd=true or false depending on max product or sum product message passing.
	 */
	public ScalableMessagePassing(Method method, boolean maxProd){
		this.method=method;
		this.maxProd=maxProd;
		if (maxProd) {
			marginalComputer=new MarginalComputer();
		} else {
			marginalComputer = new SumMarginalComputer();
		}
	}
	public ScalableMessagePassing(){
		this(true,true);
	}
	public ScalableMessagePassing(GraphicalModel model){
		this();
		if (model != null) initialize(model);
	}

	public ScalableMessagePassing(GraphicalModel model,int debug){
		this(model);
		this.debugLvl = debug;
	}
	
	// gsc: Added a single constructor to set trws and maxProd flags and call initialize(model)
	public ScalableMessagePassing(GraphicalModel model, int debug, Method method, boolean maxProd, BitSet activeNodes) {
		this(method, maxProd);
		if (model != null)
			initialize(model);
		this.debugLvl = debug;
		this.observedNodes = activeNodes;
	}
	
	//////////////////////////////////////////////////////////
	//                 Energy construction                  //
	//////////////////////////////////////////////////////////

	//////////////////////////////////////////////////////////
	//                Energy minimization starts            //
	//////////////////////////////////////////////////////////
	// Returns number of iterations. Sets lowerBound and energy.
	private int minimize_TRW_S(Properties options){
		lowerBoundPrev=lowerBound = Double.NEGATIVE_INFINITY;
		m_iterMax = Integer.parseInt(options.getProperty("maxIters","100"));
		m_eps=Double.parseDouble(options.getProperty("eps","1e-10"));
		return doMessagePassing(0);
	}
	protected int doMessagePassing(double initLB) {
		double vMin;
		int iter =0;
		LabelScoreArray Di = newLabelScoreArray(model.getMaxArity());
		// main loop
		
		EdgeIter edgeIter = null;
		int numNodes = model.getGraph().getNumNodes();
		for (iter=1; ; iter++) {
			lowerBoundPrev = lowerBound;
			lowerBound = initLB;
			for (int dir = 1; dir >= 0; dir--) {
				messageCache.startPass(dir);
				for (int node=0; node < numNodes; node++) {
					int nodeId = (dir==1)? node : numNodes-node-1;
					//gsc: Added following two lines (which avoids overriding doMessagePassing in SparseScalableMP)
					Di.clear(); 
					activateNodeLabels(Di, nodeId); // Find all active labels
					copyNodePotential(Di, nodeId);
					addMessagesToNode(nodeId, Di);
					
					/*LabelScoreArraySparse diS = (LabelScoreArraySparse)Di;
					for(int i=diS.getNumActive()-1; i>=0 ;--i)
						assert(!Double.isInfinite(diS.get(diS.getActiveLabel(i)))):""+nodeId;*/
					
					if (dir==0) {
						vMin = computeAndSubstractMin(Di, nodeId);
						lowerBound += vMin;
					}
					/*if(debugLvl > 1){
						System.out.print("1 Iter "+iter+" dir "+(dir==1?"FWD":"BKD")+" "+nodeId+" [");
						for(double a: Di.toNativeArray()){
							System.out.print(a+", ");
						}
						System.out.println("]");
					}*/		
					int nbr = -1;
					for (edgeIter=orderedGraph.getEdgeIter(nodeId, dir, edgeIter), nbr = edgeIter.nextNbr(); nbr >= 0; nbr=edgeIter.nextNbr()){
						int edgeNum = edgeIter.edgeNumber();
						messageCache.nextEdge();
						assert(messageCache.currentEdge==edgeNum);
						
						if(observedNodes != null && observedNodes.get(nbr))
							continue;
						
						vMin = updateMessages(Di, edgeNum, nodeId, nbr);//for forward pass dir is 0.		
						/*if(debugLvl > 0 && iter ==2){
							System.out.print("2 Iter "+iter+" dir "+(dir==1?"FWD":"BKD")+" "+nodeId+" --> "+nbr+" [");
							for(int i=0; i<model.getMaxArity();i++){
								System.out.print(messageCache.get(nodeId, nbr, edgeNum, i)+", ");
							}
							System.out.println("]");
						}*/
						if (dir==0) {lowerBound += vMin; } 
					}
				}
			}

			////////////////////////////////////////////////
			//          check stopping criterion          //
			////////////////////////////////////////////////
			boolean finishFlag = false;
			if (iter >= m_iterMax) 	{
				finishFlag = true;
			}
			if (m_eps >= 0)  {
//				System.out.println("Iter "+iter+" "+lowerBound+" "+lowerBoundPrev);
				if (iter > m_iterMax || Math.abs(lowerBound - lowerBoundPrev) <= m_eps) {
					finishFlag = true;
				}
				lowerBoundPrev = lowerBound;
			}

			if (debugLvl > 0) System.out.println("Iteration "+iter + " lower bound "+lowerBoundPrev);
			if ((finishFlag && method==Method.TRWS) || (finishFlag && method!=Method.TRWS)) {
				// gsc: Compute Marginal if BP is run.
				if (maxProd)
					energy = computeSolutionAndEnergy();
				else
					computeMarginalSolution();
				break;
			}
		}
		return iter;
	}
	
	protected void activateNodeLabels(LabelScoreArray di, int nodeId) {
		// Required only for Sparse Implementation
	}
	protected void findActiveLabels(LabelScoreArray di, int nodeId) {
		// Required only for Sparse Implementation
	}
	protected LabelScoreArray newLabelScoreArray(int sz) {
		return new LabelScoreArray(sz);
	}

	public double computeSolutionAndEnergy(){
		double E=0;
		LabelScoreArray DiBackward = newLabelScoreArray(model.getMaxArity());
		LabelScoreArray Di = newLabelScoreArray(model.getMaxArity());
		EdgeIter edgeIter = null;
		int numNodes = model.getGraph().getNumNodes();
		for (int nodeId=0; nodeId < numNodes; nodeId++) {
			//gsc: Added following two lines (which avoids overriding doMessagePassing in SparseScalableMP)
			DiBackward.clear();
			findActiveLabels(DiBackward, nodeId);
			copyNodePotential(DiBackward, nodeId);
			
			int nbr = -1;
			for (edgeIter=orderedGraph.getEdgeIter(nodeId, 0, edgeIter), nbr = edgeIter.nextNbr(); nbr >= 0; nbr=edgeIter.nextNbr()) {
				addColumnToMessage(DiBackward, nbr, nodeId);
			}
			// add forward edges
			copyList(Di, DiBackward);
			
			for (edgeIter=orderedGraph.getEdgeIter(nodeId, 1, edgeIter), nbr = edgeIter.nextNbr(); nbr >= 0; nbr=edgeIter.nextNbr()){
				addMessage(Di, nodeId, nbr, orderedGraph.getEdgeNumber(nbr, nodeId));
			}
			computeMin(Di,nodeId);
			if(solutionLabel[nodeId]>0)
				E += DiBackward.get(solutionLabel[nodeId]);
		}

		return E;
	}
	//gsc: Added following method for obtaining marginal solution
	protected void computeMarginalSolution() {
		LabelScoreArray Di = tempBufSrc;
		int numNodes = model.getGraph().getNumNodes();
		for (int nodeId = 0; nodeId < numNodes; nodeId++) {
			Di.clear();
			copyNodePotential(Di, nodeId);
			addMessagesToNode(nodeId, Di);
			double normalizeConst = RobustMath.LOG0;
			double logConf = Double.NEGATIVE_INFINITY;
			int predLabel = 0;
			double[] sumMargs = Di.toNativeArray();
			for (int j = 0; j < sumMargs.length; ++j) {
				double pot = sumMargs[j];
				normalizeConst = RobustMath.logSumExp(normalizeConst, pot);
				if (logConf < pot) {
					logConf = pot;
					predLabel = j;
				}
			}
			confidence[nodeId] = Math.exp(logConf - normalizeConst);
			solutionLabel[nodeId] = predLabel;
		}
	}

	//////////////////////////////////////////////////////////
	//                Energy minimization ends            //
	//////////////////////////////////////////////////////////

	protected double computeAndSubstractMin(LabelScoreArray di, int i) {
		double vMin = di.get(0);
		int nodeArity = model.getNodeArity(i);
		for (int k=1; k<nodeArity; k++)		{
			vMin = Math.min(di.get(k),vMin); 
		}
		for (int k=0; k<nodeArity; k++) {
			di.add(k, -vMin);
		}
		return vMin;
	}

	protected double computeMin(LabelScoreArray di, int nodeId) {
		double vMin = di.get(0);
		int nodeArity = model.getNodeArity(nodeId);
		int kMin = 0;
		for (int k=1; k<nodeArity; k++){
			if (vMin > di.get(k)){
				vMin = di.get(k);
				kMin = k;

			} else if (vMin==di.get(k)) {
				//System.out.println("Tie at label "+k+ " and "+kMin);
			}
		}
		solutionLabel[nodeId] = kMin;	
		return vMin;
	}

	protected LabelScoreArray tempBufSrc;
	protected LabelScoreArray tempBufDst;
	
	protected double updateMessages(LabelScoreArray di, int edgeNum, int srcNode, int dstNode) {
		int srcArity = model.getNodeArity(srcNode);
		int dstArity = model.getNodeArity(dstNode);
		tempBufDst.clear(); tempBufSrc.clear();
		double gamma = getGamma(edgeNum,srcNode,dstNode);
		
//		System.out.println(srcNode+":::"+dstNode+"=="+gamma);
		
		if (method != Method.MeanField) {
			int reverseEdgeNum = orderedGraph.getEdgeNumber(dstNode, srcNode);
			for(int i=0; i<srcArity; i++){
				tempBufSrc.add(i, gamma*di.get(i) - messageCache.get(dstNode, srcNode, reverseEdgeNum,  i) );
			}
		} else {
			for(int i=0; i<srcArity; i++){
				tempBufSrc.add(i, gamma*di.get(i));
			}
		}
		
		double vMin = Double.POSITIVE_INFINITY;
		for (int kdest=0; kdest<dstArity; kdest++) {
			double pot = marginalComputer.initValue();
			for (int ksource=0; ksource<srcArity; ksource++) {
				// gsc: TRWS / BP sign flag
				double v =  tempBufSrc.get(ksource) + (method==Method.TRWS?-1:1)*model.getEdgePotential(srcNode, ksource, dstNode, kdest);
				pot = marginalComputer.aggregate(pot, v);
			}
			tempBufDst.add(kdest,pot);
			vMin = Math.min(vMin, pot);
		}
		
		for (int kdest=0; kdest<dstArity; kdest++)
		{
			messageCache.addMessage(kdest, (tempBufDst.get(kdest)-vMin));
			//System.out.println(srcNode+"-->"+dstNode + ":"+kdest+" = "+(tempBufDst.get(kdest)-vMin)+ " src "+di.get(kdest));
		}
		
		return vMin;
	}

	// gsc: Changed the method mode from private to protected
	protected double getGamma(int edgeNum, int srcNode, int dstNode) {
		return (method==Method.TRWS?1.0/Math.max(orderedGraph.getFwdDegree(srcNode),orderedGraph.getBwdDegree(srcNode)):1);
	}

	protected void copyList(LabelScoreArray dst, LabelScoreArray src){
		dst.clear();
		dst.add(src);
	}

	protected void copyNodePotential(LabelScoreArray di, int nodeId) {
		// gsc: Commented following line because it to be always called either before 'activeNodeLabels' or 'findActiveLabels'. 
		/* di.clear(); */
		for(int j=0; j<model.getNodeArity(nodeId); j++){
			// gsc: TRWS / BP sign flag
			di.add(j, (method==Method.TRWS?-1:1)*model.getNodePotential(nodeId, j));
		}
	}

	// Returns an integer in [0,Ki). Can be called only after Minimize().
	public int getSolution(int i){
		int label = -1;
		label = solutionLabel[i];
		return label;
	}

	protected NbrIterator nbrIter=null;
	protected void addMessagesToNode(int i, LabelScoreArray Di) {
		int nbr = -1;
		for (nbrIter = model.getGraph().getNeighbourIterator(i, nbrIter), nbr = nbrIter.nextNbr(); nbr != -1; nbr = nbrIter.nextNbr()) {
			addEdgeMessages(Di, orderedGraph.getEdgeNumber(nbr, i), i, nbr);
		}
	}
	protected void addEdgeMessages(LabelScoreArray di, int edgeNum, int dst, int src) {
		for (int k=model.getNodeArity(dst)-1; k >= 0; k--)	{
			//System.out.println("adding " + di.get(k) + " to "+  messageCache.get(src,dst,edgeNum,  k));
			di.add(k, messageCache.get(src,dst,edgeNum,  k));
		}
	}

	protected void addMessage(LabelScoreArray v, int dstNode, int srcNode, int edgeNum) {
		int arity = model.getNodeArity(dstNode);
		for (int k=0; k<arity; k++)
		{
			v.add(k, messageCache.get(srcNode, dstNode, edgeNum, k));// to add message values.
		}
	}

	protected void addColumnToMessage(LabelScoreArray data, int source, int dest) {
		int k;
		int arity = model.getNodeArity(dest);
		int sourcelabel = solutionLabel[source];
		for(k = 0; k<arity; k++){
			data.add(k, (-1 *model.getEdgePotential(source, sourcelabel, dest,k)));
		}
	}
	
	public SolutionWithBounds getTopK(int k) throws UnsupportedOperationException {
		run(new Properties(), false);
		return getSolution();
	}

	public void initialize(GraphicalModel model) {
		this.model = model;
		assert(model.getMaxCliqueSize() <= 2);
		UDGraph graph = model.getGraph();
		orderedGraph = new FwdBwdGraph(graph);
		messageCache = newMessageCache();
		solutionLabel = new int[model.getGraph().getNumNodes()];
		// gsc: Initialize confidence array
		confidence = new double[solutionLabel.length];
		tempBufDst=newLabelScoreArray(model.getMaxArity());
		tempBufSrc=newLabelScoreArray(model.getMaxArity());
	}

	protected MessageCache newMessageCache() {
		return new MessageCache(model, orderedGraph);
	}

	public SolutionWithBounds getSolution() {
		return solution;
	}

	//gsc: Change visibility from private to protected. Calling getTopK() do not allow setting Options object and needs to be overriden
	protected int run(Properties options, boolean toPrint){
		minimize_TRW_S(options);
		int[] labeling = new int[model.getGraph().getNumNodes()];
		for(int i=0; i<labeling.length; i++){
			labeling[i] = getSolution(i);
		}
		solution = new SolutionWithBounds();
		solution.labeling = labeling;
		solution.score = model.getScore(labeling);
		solution.bound = -lowerBound+model.getConstantPotential();
		return 0;
	}

	public static double computeMAP(GraphicalModel model, int[] labeling, Properties options) {
		ScalableMessagePassing trws = new ScalableMessagePassing(model);
		if(options == null) {
			options = new Properties();
			options.put("maxIters","5");
		}
		trws.run(options, false);
		if(labeling != null) {
			System.arraycopy(trws.solution.labeling, 0, labeling, 0, trws.solution.labeling.length);
		}
		return -1.0*trws.lowerBound;
	}
	public double[] nodeLogMaxMarginal(int var) throws UnsupportedOperationException {
		LabelScoreArray Di = tempBufSrc;
		copyNodePotential(Di, var);
		addMessagesToNode(var, Di);
		double minMargs[] = Di.toNativeArray();
		for (int j = 0; j < minMargs.length; j++) {
			minMargs[j] *= -1;
		}
		return minMargs;
	}
	
	public double getLogPartitionValue() throws UnsupportedOperationException {
		throw new NotImplementedException();
	}

	public double[] nodeLogSumMarginal(int var) throws UnsupportedOperationException {
		LabelScoreArray Di = tempBufSrc;
		copyNodePotential(Di, var);
		addMessagesToNode(var, Di);
		double minMargs[] = Di.toNativeArray();
		for (int j = 0; j < minMargs.length; j++) {
			minMargs[j] *= (method == Method.TRWS) ? -1: 1;
		}
		return minMargs;
	}
	
	//gsc: Returns winning marginal proability of a node
	public double getConfidence(int nodeId){
		return confidence[nodeId];
	}
	
	public PotentialInterface logSumMarginal(int[] variables, PotentialInterface margClassToReuse)
	throws UnsupportedOperationException {
		throw new UnsupportedOperationException();
	}
}
