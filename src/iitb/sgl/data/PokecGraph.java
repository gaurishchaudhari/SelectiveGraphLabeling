package iitb.sgl.data;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.StringTokenizer;

import gnu.trove.list.array.TIntArrayList;
import iitb.shared.graphs.UDGraphAdjList;

public class PokecGraph extends SocialGraph implements NodeFeatured {
	
	final int NUM_ATTRIBUTES = 18;
	int[][] nodeFeatures;

	public PokecGraph(String path) throws Exception {
		String line;
		nodeLabels = new TIntArrayList() ;
		inEdgesMap = new ArrayList<TIntArrayList>();
		outEdgesMap = new ArrayList<TIntArrayList>();
		inOutEdgesMap = new ArrayList<TIntArrayList>();
		int maxAgeIndex = 0;
		
		BufferedReader br = new BufferedReader(new FileReader(path + "node_labels"));
		while((line=br.readLine())!=null){
			StringTokenizer st = new StringTokenizer(line, " ");
			st.nextToken();	// user
			int age = Integer.parseInt(st.nextToken());
			nodeLabels.add(age);
			maxAgeIndex = Math.max(maxAgeIndex, age);
			inEdgesMap.add(new TIntArrayList());
			outEdgesMap.add(new TIntArrayList());
			inOutEdgesMap.add(new TIntArrayList());
		}
		br.close();
		
		this.numLabels = maxAgeIndex + 1;
		nodeLabels.trimToSize();
		inEdgesMap.trimToSize();
		outEdgesMap.trimToSize();
		inOutEdgesMap.trimToSize();
		
		br = new BufferedReader(new FileReader(path + "edges_followersOnly"));
		while((line=br.readLine())!=null){
			StringTokenizer st = new StringTokenizer(line, " ");
			int u = Integer.parseInt(st.nextToken());
			int f = Integer.parseInt(st.nextToken());
			inEdgesMap.get(u).add(f);
		}
		br.close();
		
		br = new BufferedReader(new FileReader(path + "edges_friendsOnly"));
		while((line=br.readLine())!=null){
			StringTokenizer st = new StringTokenizer(line, " ");
			int u = Integer.parseInt(st.nextToken());
			int f = Integer.parseInt(st.nextToken());
			outEdgesMap.get(u).add(f);
		}
		br.close();
		
		br = new BufferedReader(new FileReader(path+"edges_rfriends"));
		while((line=br.readLine())!=null){
			StringTokenizer st = new StringTokenizer(line, " ");
			int u = Integer.parseInt(st.nextToken());
			int f = Integer.parseInt(st.nextToken());
			inOutEdgesMap.get(u).add(f);
		}
		br.close();
		
		this.graph = new UDGraphAdjList(nodeLabels.size());
		br = new BufferedReader(new FileReader(path + "edges"));
		while((line=br.readLine())!=null){
			StringTokenizer st = new StringTokenizer(line, " ");
			int u = Integer.parseInt(st.nextToken());
			int f = Integer.parseInt(st.nextToken());
			graph.addEdge(u, f);
		}
		br.close();
		
		System.out.println("Number of Nodes  = "+graph.getNumNodes());
		System.out.println("Number of Edges  = "+graph.getNumEdges());
		System.out.println("Number of Labels = "+numLabels);
		
		// Read Node Feature Vectors
		nodeFeatures = new int[graph.getNumNodes()][NUM_ATTRIBUTES];
		br = new BufferedReader(new FileReader(path + "node_features"));
		while((line=br.readLine())!=null){
			StringTokenizer st = new StringTokenizer(line, " ");
			int node = Integer.parseInt(st.nextToken());
			int aId = 0;
			while(st.hasMoreTokens())
				nodeFeatures[node][aId++] = Integer.parseInt(st.nextToken());
		}
		br.close();
	}
	
	public PokecGraph(String path, int percentObserved, int seed) throws Exception {
		this(path);
		observedNodes = new BitSet(graph.getNumNodes());
		observedNodes.clear();
		if(percentObserved == -1){	// set all nodes as observed
			for(int node = graph.getNumNodes()-1; node>=0; --node)
				observedNodes.set(node);	
		}
		else{
			BufferedReader br = new BufferedReader(new FileReader(path + "observedNodes_" + percentObserved + "_"+seed));
			String line;
			while((line=br.readLine())!=null){
				int an = Integer.parseInt(line);
				observedNodes.set(an);
			}
			br.close();
		}
		
		this.numObservedNodes = observedNodes.cardinality();
		System.out.println("Number of ObservedNodes = "+numObservedNodes);
	}
	
	@Override
	public void getNodeAttributes(int nodeId, int[] tmpNodeIdAttributes){
		int j = 0;
		for(int aId = 0; aId < tmpNodeIdAttributes.length; ++aId){
			if(nodeFeatures[nodeId][aId] == 1)
				tmpNodeIdAttributes[j++] = aId;
		}
		if(j < tmpNodeIdAttributes.length)
			tmpNodeIdAttributes[j] = -1;
	}
	
	@Override
	public int numNodeAttributes(){
		return nodeFeatures[0].length;
	}
}

