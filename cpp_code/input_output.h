#ifndef INPUT_OUTPUT_H
#define INPUT_OUTPUT_H

// This is a file which deals with the input and output of the 
// code written, that is it can generate and read input,
// as well as write and display output

#include <random>
#include "line_geometry.h"
#include "simple_svg_1.0.0.hpp"
#include "cgal_interface.h"

#define TOLERANCE 0.000001



// A wrapper function which checks whether the graph defined by 
// the adjacency matrix edges is connected
bool isConnected(vector<vector<bool> > &edges);

// Recursive implementation of a depth first search; it stores the blocks inside of visited
void depthFirstSearch(vector<vector<bool> > &edges, vector<int> &visited, int vertexId, int group);

// A new function which creates a graph from the top bottom; that is 
// if creates a random set of n vertices, creates the complete graph, and then proceeds
// to create a random connected planar graph. It first removes all intersecting edges.
// It then adds edges to make the graph connected. Finally it randomly chooses how many edges
// to remove and removes them whilst keeping the graph connected.
void generate_graph_top_bottom(Graph* graph, int n, double X, double Y, int seed);

// Function used to generate a random graph
// It is not ideal and maybe a new function should be looked at
// Specifies number of vertices (n), edges (e), max x and y coordinates (X,Y respectively), and the seed of the 
// random function to reproduce results. Everything is stored inside of graph.
// An extra option specifies how long the edges can be (to not get "weird" graphs with long edges)
void generate_graph(Graph* graph, int n, int e, double X, double Y, int seed, double maxEdgeLength);

// Reads graph from a file input
int read_graph(Graph* graph, string filename);

// Writes graph into a file
void write_graph(Graph* graph, string filename);

// Writes out graph into an SVG file. Need to specify the size of the disaster and the dimensions of the graph
// A bunch of options are give:
//				- showLeaves: Shows nodes which are part of a leaf as a red circle rather than black
// 				- showCuts: Draws the edges that are in a cut in red
// 				- showMinCuts: Draws the edges that are in a minimal cut in red
//				- showVulnerableRegions: Draws the parallelogram of the vulnerable region of every given
//										 cut (if showCuts = true) and minimal cut (if showMinCuts = true)
//				- ShowDestructionRegions: Same as showVulnerableRegions, but with destruction regions
//				- showAuxEdges: If false, does not label the vertices that are part of a new edge
//				- showDisaster: Prints the disaster at the bottom right of the graph 

void write_graph_to_svg(Graph* graph, double length, string filename, double dimX, double dimY, 
								bool showLeaves, bool showCuts, bool showMinCuts,
								bool showVulnerableRegions,
								bool showDestructionRegions,
								bool showAuxEdges, bool showDisaster);

void visibilityToSVG(string filename, Face_Handle &face, Poly_Arr visibility);

#endif