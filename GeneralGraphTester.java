/*	Author:	Benjamin Pipes
 * 
 * 	This file contains test cases and room for sample executions
 * 	to demonstrate the GeneralGraph structure.
 *
 *  Files needed: GeneralGraphTester.java, GeneralGraph.java
 *  
 */

import java.io.FileNotFoundException;


public class GeneralGraphTester {

	public static void main(String[] args) throws FileNotFoundException {
		
		System.out.println("Beginning execution...");
		
		GeneralGraph graph = GeneralGraph.buildTsp("TSPTESTCASE.dat", "TSP Test Case");
		graph.displayNodes();
		graph.performTSP(graph);
		graph.performBellmanFord(graph.indexOf("C"));
			
		GeneralGraph graph1 = GeneralGraph.buildTc("TCTESTCASE.dat", "Transitive Closure Test Case");
		graph1.displayNodes(); //should not have names for node values
		graph1.performWarshalls();
		graph1.printWarshallsMatrix();
		
		GeneralGraph graph2 = GeneralGraph.buildBf("BFTESTCASE.DAT", "Bellman Ford Test Case");
		graph2.displayNodes();
		
		int source = graph2.indexOf("Appl"); //find index of town
		graph2.performBellmanFord(source); //perform BF with that node as source
		
		graph2.printAdjMatrix();
		graph2.performWarshalls();
		graph2.printWarshallsMatrix();
		
		graph2.displayNodes();
		
		graph2.printMatrix();
		graph2.performTSP(graph2);

		System.out.println("Program terminated");
	}

}
