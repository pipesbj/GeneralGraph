
/*	Author: Benjamin Pipes
 * 	
 * 	This is a general graph structure that can provide
 * 	Traveling salesman, Transitive closure, and Bellman-Ford
 * 	data.
*/


import java.io.File;
import java.io.FileNotFoundException;
import java.util.Scanner;
import java.util.Stack;
import java.util.regex.Pattern;


public class GeneralGraph {
	private static Object[] nodes;	//node names
	private static int[][] weightMatrix; //stores weight values between two points
	private static int[][] weightMatrixun; //undirected, for TSP
	private static int inf = (int) Double.POSITIVE_INFINITY; //infinity value
	private static int size = 0; //node (vertex) count
	private static String name = ""; //local name of graph. for ease of understanding
	private static Object[] nodesTSP; //stores local node names for TSP (default: a, b, c, d...)
	private static int[][] adjMatrix; //stores adjacency values
	private static int[][] warshallsMatrix; //stored transitive closure values
	private static Object[] bfpath; //holds the desired cycle or path a user may be asking for
	
	
	/*	Constructor
	 * 	Input: 	number of nodes (size of graph)
	 * 			local name of graph. for display purposes
	 * 	Constructs a weight matrix, loads the
	 * 	node names into array, sets all values to infinity
	 * 
	 */
	public GeneralGraph(int n, String s){                                                                                                                                                                                                                                                                                                                     
		
		  size = n;
		  name = s;
		  nodes = new Object[n];
		  nodesTSP = new Object[n];
		  bfpath = new Object[n];
		  weightMatrix = new int[n][n];
		  weightMatrixun = new int[n][n];
		  adjMatrix = new int[size][size];
		  warshallsMatrix = new int[size][size];
		  
		  buildTSPnodeSet();
		  
		    for (int i = 0; i < n; i++){
		    	for (int j = 0; j < n; j++){
		    		weightMatrix[i][j] = inf;
		    		weightMatrixun[i][j] = inf;
		    	}
		    }
	}//constructor 
	
	/*	Constructor
	 * 	Input: 	number of nodes (size of graph)
	 * 	Constructs a weight matrix, loads the
	 * 	node names into array, sets all values to infinity
	 * 
	 */
	public GeneralGraph(int n){                                                                                                                                                                                                                                                                                                                     
		
		  size = n;
		  name = "NO_NAME";
		  nodes = new Object[n];
		  nodesTSP = new Object[n];
		  bfpath = new Object[n];
		  weightMatrix = new int[n][n];
		  weightMatrixun = new int[n][n];
		  adjMatrix = new int[size][size];
		  warshallsMatrix = new int[size][size];
		  
		  buildTSPnodeSet();
		  
		    for (int i = 0; i < n; i++){
		    	for (int j = 0; j < n; j++){
		    		weightMatrix[i][j] = inf;
		    		weightMatrixun[i][j] = inf;
		    	}
		    }
	}//constructor 
	
	/*	input: source, destination, value
	 * 	output: adds an edge to all structures
	 * 			that Bellman-Ford uses
	 *
	 * int input for node positions (0 = node1, 1 = node2 ...)
	 */
	public void bfa(int from, int to, int val){	
		weightMatrix[from][to] = val;
	}
	
	/*	input: source, destination, value
	 * 	output: adds an edge to all structures
	 * 			the transitive closure algorithms use
	 * int input for node positions (0 = node1, 1 = node2 ...)
	 */
	public void tca(int from, int to, int val){
		adjMatrix[from][to] = val;
		warshallsMatrix[from][to] = val;
	}
	
	/*	input: source, destination, value
	 * 	output: adds an edge to all structures
	 * 			the TSP algorithm uses
	 * int input for node positions (0 = node1, 1 = node2 ...)
	 */
	public void tspa(int from, int to, int val){
		weightMatrixun[from][to] = val;
		weightMatrixun[to][from] = val;
	}
	
	/*  input: source, destination, value, mode
	 * 			*mode indicates which set of secondary matrices
	 * 			*are to be handled
	 * 			*0 - tc, bf (primary: tsp)
	 * 			*1 - tc, tsp (primary: bf)
	 * 			*2 - tsp, bf (primary: tc)
	 * 				
	 * 	output: adds an edge to all structures
	 * 			the secondary operations use
	 */
	public void bldSec(int f, int t, int n, int mode){
		switch(mode){
		case 0 :
			tca(f, t, 1);
			bfa(f, t, n);
		case 1 :
			tca(f, t, 1);
			tspa(f, t, n);
		case 2 :
			tspa(f, t, n);
			bfa(f, t, n);
		}
		
	}
	
	/*	For the user to call for the number of nodes in the graph.
	 * 
	 */
	public int size(){
		return size;
	}
	
	/*	for the user to call for the weightmatrix
	 * 	0 - directed matrix 
	 * 	1 - unndirected matrix
	 */
	public int[][] wgtMtrx(int matrix){
		if (matrix == 0)
			return weightMatrix.clone();
		else
			return weightMatrixun.clone();
	}
	
	/*	For Bellman Ford implementations
	 * 	implemented to build from a file that is formatted as:
	 * 	<node count>
	 * 	<node1> <weight1.1> <weight1.n>
	 * 	<noden> <weightn.1> <weightn.n>
	 * 	
	 * builds a weight matrix
	 */
	public static GeneralGraph buildBf(String f, String nm) throws FileNotFoundException{
		
		int from;
		int to;
		int wgt = inf;
		int k = 0;
		String line;		
		Scanner in = new Scanner(new File(f));
		
		//read in first element, integer number of nodes
		GeneralGraph graph = new GeneralGraph(in.nextInt(), nm);
		in.nextLine();
		
		//loads the node names into an array
			while(in.hasNext()){
				line = in.nextLine();
				for(int i = 0; i < line.length(); i++){
					if(Pattern.matches("[0-9]|-", String.valueOf(line.charAt(i)))){
						nodes[k] = line.substring(0, i - 1);
						k++;
						break;
					}
				}
			}
			
		//reset file.
		//proceed through lines to get edge values
		in = new Scanner(new File(f));
		in.nextLine();
		
		int j = 0;
		String evStr;
		String edgevals[] = new String[size]; //array to hold edge values for each line.
		k = 0;
		
		//creates the edges of the graph
		while (in.hasNext()){
			line = in.nextLine();
			
			for(int i = 0; i < line.length(); i++){
				
				if(Pattern.matches("[0-9]|-", String.valueOf(line.charAt(i)))){
					//if end of node name is found, load rest of string as values for edges
					evStr = line.substring(i);
					edgevals = evStr.split(" ");
					break;
				}
			}
					from = j; 	//set from node to index corresponding to town. (0 = node1, 1 = node2, etc)
					to = k;		//possible destination from source node
					for(int g = 0; g < edgevals.length; g++){
						to = k; //change destination node
						if(edgevals[g].compareTo("-") != 0){
							//if weight is a valid weight
							wgt = graph.strToInt(edgevals[g]);
							
							//add edge to primary matrices
							graph.bfa(from, to, wgt);
							
							//add edge to secondary matrices
							graph.bldSec(from, to, wgt, 1);
							
						}
						else{
							wgt = inf;
							graph.bfa(from, to, wgt);
						}
						k++; //increment next possible destination counter.
					}			
			j++; //new source node
			k = 0; //reset destination node counter
			
		}//while not EOF
		
		in.close();
		return graph;
	}
	
	/*	For TSP implementations
	 * 	implemented to build from a file that is formatted as:
	 *	<node count>
	 * 	<node1> <weight1.1> <weight1.n>
	 * 	<noden> <weightn.1> <weightn.n>
	 */
	public static GeneralGraph buildTsp(String f, String nm) throws FileNotFoundException{
		
		int from;
		int to;
		int wgt = inf;
		int k = 0;
		String line;		
		Scanner in = new Scanner(new File(f));
	
		GeneralGraph graph = new GeneralGraph(in.nextInt(), nm);
		in.nextLine();

		//loads the node names into an array
		while(in.hasNext()){
			line = in.nextLine();
			for(int i = 0; i < line.length(); i++){
				if(Pattern.matches("[0-9]|-", String.valueOf(line.charAt(i)))){
					nodes[k] = line.substring(0, i - 1);	//load node name into respective array index
					k++;
					break;
				}
			}
		}
		
		//reset file.
		//proceed through lines to get edge values
		in = new Scanner(new File(f));
		in.nextLine();
		
		int j = 0;
		String evStr;
		String edgevals[] = new String[size]; //array to hold edge values for each line. String data type
		k = 0;
			
		//creates the edges of the graph
				while (in.hasNext()){
					line = in.nextLine();
					
					for(int i = 0; i < line.length(); i++){
						
						if(Pattern.matches("[0-9]|-", String.valueOf(line.charAt(i)))){
							//if end of node name is found, load rest of string as values for edges
							evStr = line.substring(i);
							edgevals = evStr.split(" "); //split rest of string into array
							break;
						}
					}
							from = j; 	//set from node to index corresponding to town. (0 = node1, 1 = node2, etc)
							to = k;		//possible destination from source node
							for(int g = 0; g < edgevals.length; g++){
								to = k; //change destination node
								if(edgevals[g].compareTo("-") != 0){
									//if weight is a valid weight
									wgt = graph.strToInt(edgevals[g]);
									
									//add edge to primary matrices
									graph.tspa(from, to, wgt);
									
									//add edge to secondary matrices
									graph.bldSec(from, to, wgt, 0);
									
								}
								else{
									wgt = inf;
									graph.tspa(from, to, wgt);
								}
								k++; //increment next possible destination counter.
							}			
					j++; //new source node
					k = 0; //reset destination node counter
					
				}//while not EOF
				in.close();
		return graph;

	}
	
	/* For Warshalls transitive closure implementations
	 * implemented to build an adjacency matrix
	 * from a file formatted as:
	 * 0010
	 * 1001
	 * 0010
	 * 1000
	 * 
	 */
	public static GeneralGraph buildTc(String f, String nm) throws FileNotFoundException{
		
		return buildAdjMat(f, nm);
	}//constructor 
	
	/*	Prints the node names to the console window
	 * 
	 */
	public void displayNodes(){
		System.out.println("\nNodes values for graph \"" + name + "\":");
		for(int i = 0; i < size; i++)
			System.out.print("\"" + nodes[i] + "\" ");
		System.out.println("\n");
	}
	
	/*	Returns a String array of node names
	 * 
	 */
	public Object[] getNodes(int mode){
		if (mode == 0)
			return nodes.clone();
		else
			return nodesTSP.clone();
		
	}
	
	/*	Returns initial node set for TSP operations
	 * 
	 */
 	public Object[] getNodeSet(){
		//returns a clone of the master nodes array
		//otherwise both are in
		//same address.
		return nodesTSP.clone();
	}
	
 	/*	Returns the value for a given edge in a given matrix to access
 	 * 
 	 */
	public int getEdgeWeight(int from, int to, int mode){
		//mode: 0 - TSP call, unweighted graph
		//		1 - the other calls, weighed graph
		switch(mode){
		case 0 :
			if (weightMatrixun[from][to] != inf){
				return weightMatrixun[from][to];
				}
			else
				return inf;
		case 1 :
			return weightMatrix[from][to];
		case 2 :
			return warshallsMatrix[from][to];
		case 4 :
			return adjMatrix[from][to];
		default:
			return inf;
		}
			
		
	}

	/*	Input: 	int, the source node to shortest paths from.
	 * 	output:	void, prints the paths from source to each destination, if possible.
	 * 
	 */
	public void performBellmanFord(int source){
		
		
		@SuppressWarnings("unchecked")
		Stack<String>[] pathSeq = new Stack[size];		//array of paths, from source node to each destination node
		Object[] temp = nodes.clone();	//holds the list of node names
		int currWgt = 0;							//the current weight being checked
		Object[] predNodes = nodes.clone(); //predecessor nodes
		int chainDist = 0; 				//distance between a sequence of nodes
		int[] distances = new int[size];	//distances from source node to each other node.
	
		for (int j = 0; j < distances.length; j++){
			//fill source distances with values
			if(weightMatrix[source][j] >= 0)
				distances[j] = weightMatrix[source][j];	//fill distances from source to dest based on source row values
			else
				distances[j] = inf;						//if no values present, set to infinity
	
			pathSeq[j] = new Stack<String>();	//create array of stacks
			pathSeq[j].add((String) temp[source]);		//add source node for each path sequence
			
		}
		
		for (int i = 0; i < weightMatrix.length; i++){
			for (int j = 0; j < weightMatrix[0].length; j++){
				//traverse through weight matrix
				currWgt = weightMatrix[i][j];
				
				if (currWgt != inf){
					//if edge has a valid weight
					if(i == source && currWgt >= 0){
						//if we are on the source node's row, add neighboring distances
						distances[j] = currWgt;
						pathSeq[j].push((String) temp[j]);	//add next node in path to respective path list
					}
					else{
						if(distances[i] != inf){
							chainDist = distances[i] + currWgt;
							
							if(chainDist < distances[j] && j != source && currWgt >= 0){
								//if new chain distance is better than existing distance
								//we aren't referring to the distance to the source node
								//the weight being assessed is not negative
								String pred = (String) temp[j];
								distances[j] = chainDist;
								predNodes[j] = (String) temp[i];

								while(pathSeq[j].size() > 1){
									//if new chained path is found, pop all trash nodes off
									pathSeq[j].pop();
								}
									//push new path sequence on
								pathSeq[j].push((String) predNodes[i]);
								pathSeq[j].push((String) predNodes[j]);
								pathSeq[j].push(pred);
						
							}//if chain distance beats old distance
						}
					}
				}//if valid weight	
			}//for col
		}//for row
		
		//print distances from the source node to destinations
		System.out.println("***RESULTS ARE TO BE READ FROM RIGHT TO LEFT. ***\n"
						+  "***RIGHTMOST NODE IS THE SOURCE NODE. ***");
		String format = "\n%14s%4s%-14s%-6s%-1d%19s%-14s\n";
		int hopCt = 0;
		int pathInd = 0;
		for (int i = 0; i < distances.length; i++){
			hopCt = 0;
			pathInd = 0;
			if(distances[i] != inf){
				System.out.format(format, temp[source], " -> ", temp[i], ":...: ", 
						distances[i]," :...: predecessor: ", predNodes[i]);	
			
				//print path
				String tmp = "";
				System.out.print("\t");
				while(!pathSeq[i].isEmpty()){
					tmp = (String) pathSeq[i].peek();
					if(tmp.compareTo((String) temp[i]) == 0){
						//if we peek at the stack and see the last node (the source node)
						//We know we are at the end of the path.
						bfpath[pathInd] = pathSeq[i].pop();
						System.out.print(bfpath[pathInd]);
						 
						
					}
					else{
					System.out.print(" < ");
					bfpath[pathInd] = pathSeq[i].pop();
				System.out.print(bfpath[pathInd]);
					}
					
					pathInd++;
					if(!pathSeq[i].isEmpty())
						hopCt++;
					
					if(!pathSeq[i].isEmpty() && tmp.compareTo((String) pathSeq[i].peek()) == 0){
						//if there is a dupe node in the path, disregard
						pathSeq[i].pop();
					}
				}
					System.out.println("(source) \n\t**hops: " + hopCt);
			}
		}

		System.out.println("\n**********************************************************");
	}//end bellman ford
	
	/*	String to int.
	 * 	input:	String
	 * 	output:	int (a number value)
	 */
	public int strToInt(String s){
		//string to int
		int num = 0;
		int power;
		char pos[] = s.toCharArray();
		if(s.charAt(0) != '-'){

			if (pos.length == 1){
				//if input is a positive
				return pos[0] - '0';
			}
			else
				for (int i = 0; i < pos.length; i++){
					power = pos.length - i - 1;
					num += (pos[i] - '0') * Math.pow(10, power);
				}
			
		}
		else{
			//if input is a negative
			if (pos.length == 2){
				return (pos[1] - '0') * -1;
			}
			else
				for (int i = 1; i < pos.length; i++){
					power = pos.length - i - 1;
					num -= (pos[i] - '0') * Math.pow(10, power);
				}
		}
		return num;
	}

	/*	Input:	filename, referencing the adjacency matrix file
	 * 	output:	none, sends adjacency & Warshalls
	 *  matrices to static variables
	 * 
	 * 	implemented to build from a file that is formatted as
	 * 	an adjacency matrix.
	 */
	public static GeneralGraph buildAdjMat(String f, String nm) throws FileNotFoundException{
		
		String line;		
		Scanner in = new Scanner(new File(f));
		size = 0;
		line = in.nextLine();
		for (int i = 0; i < line.length(); i++){
			//walks through first line, increments the size
			//to be used in creating matrix
			size++;
		}
		in.close();
		
		//close and reopen the file
		in = new Scanner(new File(f));
			
		//create adjacency & transitive matrix
		adjMatrix = new int[size][size];
		warshallsMatrix = new int[size][size];
		GeneralGraph graph = new GeneralGraph(size, nm);
		line = in.nextLine();
		
		for (int i = 0; i < size; i++){
			for (int j = 0; j < size; j++){
				//traverses through file
				//load values to matrix(s)				
				graph.tca(i,j,line.charAt(j) - '0');
			}
			//check if another line exists, move to next line if so
			if(in.hasNextLine())
			line = in.nextLine();
		}//build adjacency matrix
		
		in.close();
		return graph;
	}//build from adjacency matrix
	
	/*	Input: none, uses static adjacency matrix
	 *  output: none, uses static matrix to append information about the 
	 * 			transitive closure
	 * 
	 * 	Implements Warshall's transitive closure algorithm 
	 */
	public void performWarshalls(){
		
		int edge1 = 0;
		int edge2 = 0;
		int edge3 = 0;

		for (int i = 0; i < size; i++){

			for (int j = 0; j < size; j++){

				for (int k = 0; k < size; k++){

					edge1 = warshallsMatrix[j][k];
					edge2 = warshallsMatrix[j][i];
					edge3 = warshallsMatrix[i][k];
					if (edge1 == 1 || (edge2 == 1 && edge3 == 1))
						warshallsMatrix[j][k] = 1;
					else
						warshallsMatrix[j][k] = 0;
				}
			}	
		}
	}//Warshall's algorithm
	
	/*	Performs Traveling Salesman Problem for a graph
	 * 	Input:	a graph
	 * 	output:	void, prints out shortest hamiltonian cycle
	 * 
	 */
	public void performTSP(GeneralGraph graph){
		
		  Object[] nodeSet = graph.getNodeSet(); 	//clone of the set of nodes. String names
		  boolean exhausted = false; 	// flag to stop the algorithm
		  int to, from; 				// values for source and destination nodes
		  Object startNode; 			// the source node for all possible TSP paths
		  Object[] optimalPath = null; 	// an array of node names that make up the optimal path
		  int currWeight = 0; 			// weight value of current edge
 		  int pathWeight = 0;			// weight of the path being analyzed
		  int bestWeight = 0;			// optimal weight for each decision made at a node
		  boolean pathExists = true;
		  boolean incGraph = false;		//flag for incomplete graph. For displaying info
		  bestWeight = Integer.MAX_VALUE;	//best weight is initially max integer
		  optimalPath = nodeSet.clone();
		  
		  //Runs each permutation and tracks the number of permutations done
		  //Exhausted should be true when the 2nd node becomes the 1st
		  while (!exhausted) {
			    //Resets for each permutation
			   startNode = nodeSet[0];
			   pathWeight = 0;
			   pathExists = true;
			   
			   //Goes through each permutation of an array, save the path with the 
			   //lowest weight.
			   for (int i = 0; i < nodeSet.length; i++) {
				   
				    from = Character.valueOf((char)nodeSet[i]) - 97;
				    if (i == nodeSet.length - 1)
				    	to = (char) 0;
				    else
				    	to = Character.valueOf((char)nodeSet[i + 1]) - 97;
				    // if nodeSet index is at the end of the path, the "to"
				    // variable is set to the starting node, to complete the cycle
				
				    
				    currWeight = graph.getEdgeWeight(from, to, 0);
				    //next comparison solely exists to account for incomplete
				    //graphs. if graph is complete, this compare is extraneous
				    if (currWeight == -1 || currWeight == inf) {
				     // if returned with invalid value, the edge DNE, therefore the path DNE
				    	pathExists = false;
				    	incGraph = true;
				    	pathWeight = 0;
				    	//System.out.println("non-complete graph");
				    }
				
				    else {
				     // else, the path exists, and weight (distance) is
				     // accumulated
				    	pathWeight += currWeight;
				    }
			
			   }// for the path array, (the node set)
			
			   //System.out.println("comparing: " + pathWeight + " <> best: " + bestWeight );
			   if (pathWeight < bestWeight && pathExists) {
				   bestWeight = pathWeight;
				   optimalPath = nodeSet.clone();
			   }
			
			   // loads next node permutation into nodeSet array
			   nodeSet = nextPermutation(nodeSet);

			   if (nodeSet[0] != startNode)
				   exhausted = true;
		   // end loop condition
		  }// while TSP not exhausted
		  
		  //*******************************
		  // Print Results
		  //*******************************
		  
		  if(incGraph)
			  System.out.println("Graph \"" + name +"\"" + " is an incomplete graph.");
		  
		  for (int i = 0; i < optimalPath.length; i++) {
			  System.out.print(nodes[(char)optimalPath[i] - 97] + " > ");
		  }
		  if (pathExists)
			  System.out.println(nodes[(char)optimalPath[0] - 97] + " : " + bestWeight);
		  System.out.println("optimal distance: " + bestWeight);
	}
	
	/*	Builds a set of nodes represented as a, b, c, d, e, f... regardless
		of official node names. For TSP algorithm
	*/
	public static void buildTSPnodeSet(){
		
		for(int i = 0; i < size; i++){
			nodesTSP[i] = Character.toChars(i+97)[0]; //returns array. pick first index
		}
	}
	
	/*	generates the next permutation for a set of values
	 * 	Ex.		input: 	[b, a, c, d]
	 * 			output:	[b, a, d, c] 
	 * 
	 */
	public static Object[] nextPermutation(Object[] nodeSet) {
		  // input: array of characters
		  // output: next permutation of those characters
		  // Ex. >bacd -> badc
		  int incLoc;
		  int length = nodeSet.length - 1;
		  int i = length;
		  char temp;
		  Object end[] = { '\0' };

		  while (i > 0 && (char)nodeSet[i - 1] > (char)nodeSet[i]) {
			  i--;
		  	}
		  incLoc = i - 1;
		  if (incLoc < 0)
		   return end;

		  int j = length;
		  while ((char)nodeSet[j] < (char)nodeSet[incLoc]) {
			  j--;
		  	}
		  
		  temp = (char)nodeSet[j];
		  nodeSet[j] = nodeSet[incLoc];
		  nodeSet[incLoc] = temp;

		  int k = length;
		  while (i < k) {
			   temp = (char)nodeSet[i];
			   nodeSet[i] = nodeSet[k];
			   nodeSet[k] = temp;
			   i++;
			   k--;
		  	}

		  return nodeSet;

		 }// nextPermutation
	
	//for debugging purposes. Prints n x n matrix
	public void printMatrix(){
		
		
		for(int k = 0; k < nodes.length; k++){
			System.out.print(nodes[k] + " ");
		}
		
		System.out.println("\n***********Displaying Matrix*************");
		for (int i = 0; i<weightMatrix.length; i++){
			for (int j = 0; j < weightMatrix[i].length; j++){
				
				if(weightMatrixun[i][j] == inf){
					System.out.print("inf ");
				}
				else
					System.out.print(weightMatrixun[i][j] + " ");
			}
			System.out.println();
		}
		
	}//printMatrix
	
	//for debugging purposes. Prints the inputted n x n matrix
	//representing the given adjacency matrix
	public void printAdjMatrix(){
		
		System.out.println("Input - Adjacency Matrix: n = " + size + ", name = \"" + name + "\"");

		for (int i = 0; i<adjMatrix.length; i++){
			for (int j = 0; j < adjMatrix[i].length; j++){
				System.out.print(adjMatrix[i][j] + " ");
			}
			System.out.println();
		}
	}//printAdjMatrix
	
	/*	
	Print current n x n matrix for Warshall's algorithm
	representing the current progress in determining 
	transitive closure. If called at end of algorithm, will print the
	final transitive closure matrix
	 */	
	public void printWarshallsMatrix(){
			
		System.out.println("Transitive Matrix (Warshall's):");
		
		for (int i = 0; i<warshallsMatrix.length; i++){
			for (int j = 0; j < warshallsMatrix[i].length; j++){
				System.out.print(warshallsMatrix[i][j] + " ");
			}
			System.out.println();
		}
	}//printMatrix
	
	/*	Sequential search of the requested node via string comparison
	 * 	input:	A requested string to find 
	 * 	output:	The index of that node
	 */
	public int indexOf(String s){
		
		
		for(int i = 0; i < size; i ++){
			if (s.toUpperCase().compareTo(((String) nodes[i]).toUpperCase()) == 0){
				return i;
			}
		}
		
		System.err.println("requested string (\"" + s + "\") not found in list of node names.");
		displayNodes();
		
		System.exit(0);
		return -1;
	}
	
}
