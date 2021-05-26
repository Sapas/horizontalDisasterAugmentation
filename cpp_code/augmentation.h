#ifndef AUGMENTAION_H
#define AUGMENTAION_H

#include <vector>
#include <iostream>
#include <fstream>
#include <random>
#include <stdio.h>
#include <time.h>
#include <utility>
#include <cstring>
#include <ctime>
#include <queue>
#include <math.h> 
#include <limits>
#include "input_output.h"
#include "line_geometry.h"
#include "sweep_binary_search_tree.h"
#include "cgal_interface.h"

using namespace std;

#define TOLERANCE 0.000001

// This file is the main one; it is used to add edges to a graph to make it 
// l-resilient, where l is the length of a disaster, and a disaster is a 
// horizontal line. Each approach is explained in the function corresponding
// to it, but note that all approaches are heuristics. There are two main 
// components in this file: Finding the l-cuts of the graph, and then
// adding edges to make the graph l-resilient. For questions and support 
// of the code, please contact Nico at nandresthio@gmail.com


// A sweep algorithm is implemented to find all of the cuts; these 
// events are the possible ones during the sweep algorithm
struct Event{
	int id;		// Id of event, used to find events in queue
	char type;	// Options: 'a' activaton
			 	// 		    'd' deactivation
				//			'c' close enough
	double height; 	// Where the event happens
	Line* line1;    // Line being activated/deactivated or one of near lines
	Line* line2;	// Other line (in case of close enough)
	int score;		// Importance score of the even, lowest the higher the priority
	
	bool operator==(const Event& otherEvent) const{
		return (type == otherEvent.type) &&
				(std::abs(height - otherEvent.height) < TOLERANCE) &&
				((line1->id == otherEvent.line1->id && line2->id == otherEvent.line2->id) || 
					(line1->id == otherEvent.line2->id && line2->id == otherEvent.line1->id));
	}

	bool operator<(const Event& otherEvent) const{
		// If not in line vertically, smaller is whichever is higher
		if(std::abs(height - otherEvent.height) > TOLERANCE){return height < otherEvent.height;}
		// Same height, check score of event. If same, either one can take preference,
		// break tie by id of event for consistency.
		if(score == otherEvent.score){return id < otherEvent.id;}
		return score < otherEvent.score;
	}

	Event(int id, char type, double height, Line* line1, Line* line2):
	id(id),
	type(type),
	height(height),
	line1(line1),
	line2(line2)
	{
		if((std::abs(line1->start->y - line1->end->y) > TOLERANCE) && type == 'd'){score = 5;}
		else if((std::abs(line1->start->y - line1->end->y) < TOLERANCE) && type == 'a'){score = 4;}
		else if((std::abs(line1->start->y - line1->end->y) < TOLERANCE) && type == 'd'){score = 3;}
		else if((std::abs(line1->start->y - line1->end->y) > TOLERANCE) && type == 'a'){score = 2;}
		else if(type == 'c' && line2 != NULL && 
									std::abs(line1->start->y - line1->end->y) > TOLERANCE &&
									std::abs(line2->start->y - line2->end->y) > TOLERANCE){score = 1;}
		else{
			score = 0;
			printf("This score should not be assigned!\n");
		}
	}
};

// Next two functions are used to find the
// number of ones in binary representation
// of a number. Has O(1) lookup, O(n) done 
// when code starts. Taken from
// https://www.geeksforgeeks.org/count-set-bits-in-an-integer/

// Initialises a lookup table used to know what a certain binary bit is 
// for a given int. Used when looking at subsets
void initialize();

// Function to return the count of set bits in n
int countSetBits(int n);

// Implementation of the cantor function, which provides a unique int for
// any pair of ints
int cantorPair(int n1, int n2);

// A function which takes two vectors with ints, and returns a vector with 
// ints such that for any two indeces, the returned vector has the same value
// only if the input vectors have the same value at the indeces. In order for 
// the ints to not blow up, after using the cantor coder, the ints are 
// relabelled.
void vectorIntEncoder(vector<int> &currentVector, vector<int> &addedVector);

// Simple function which checks if a graph is connected. Runs a depth first search
// and checks whether every node can be reached.
bool isConnected(Graph* graph);

// Function which checks whether a vector is a subset of another.
// Assumes the integers in the two vectors are unique and are increasing
bool isSubset(vector<int> &smallVector, vector<int> &bigVector);

// Given a new set of nodes (potentialLeaf) which could be leaf, it gets compared
// to the current leaves of the graph. The current leaves of the graph are updated
// if no existing leaf is a subset, or if this new leaf is a subset of existing leaves
void checkPotentialLeaf(Graph* graph, vector<int> &potentialLeaf, vector<int> &minCut);

// Wrapper function used to find the leaves, simply cycles through the minimal cuts,
// finds the shores of the minimal cut and checks whether either one (or both) are
// leaves
void findLeaves(Graph* graph);

// Does a depth first search of the graph, ignoring the set of edges forbidden.
// It labels all the vertices in finds in visited with the label group
// This is a recursive function, the input vertexId specifies the current location
// in the graph of the algorithm
void depthFirstSearch(Graph* graph, vector< pair <int,int> > &forbidden, vector<int> &visited, int vertexId, int group);

// Function which checks whether the input graph (in the form of an adjacency matrix) is bipartite.
// This is not used by the current code but I think it was a pretty smart idea so I am leaving it here.
// Originally it was used to check whether a group of edges whose removal disconnected the graph made up
// an edge cut; this was only the case if the condensation of the graph was bipartite.
bool checkBipartite(vector<vector<int> > &adjacencyMatrix, vector<int> &group, int vertexId, int groupId);

// A function which checks every edge in graph to check whether it is an l-cut or not.
// It simply performs a DFS (ignoring an edge) and checks if the graph is connected!
void checkSingleEdgesCut(Graph* graph);

// A function which cycles through all of the cuts of the graph (which have already been found)
// and finds the blocks obtained when removing the cut. Also takes the chance to find which cuts are
// minimal l-cuts; these are the ones that produce exactly two blocks.
void findLabelsAndMinCuts(Graph* graph);

// Function called to analyse a set of edges which can be destroyed at the same time.
// Function finds all the subsets of the set of edges which include left, right and middle
// which are l-cuts. It starts by finding the connected components of the graph after removing 
// the edges in vectorEdgesId. It then looks at every partition of the blocks into two sets,
// knowing every set of edges going between the two is a subset of the set of edges and can be 
// destroyed by a disaster, thus being an l-cut.
void analyseEdgeSet(Graph* graph, vector< pair <int,int> > &vectorEdgesId, int left, int right, int middle);


// Function which deals with the activate event in the sweep algorithm which finds the l-cuts of the graph
// When an edge e is activated, the function finds any pair of edges e_l and e_r (to the left and right
// of e respectively) which are closer than length, the disaster length. For each such pair, the set of edges
// E between the two (including e_l and e_r) is analysed by analyseEdgeSet, looking for subsets of E which include
// e_l, e_r and e which are cuts (and thus l-cuts) of G. The function also finds every edge e' which is closer than length
// to e and calls analyseEdgeSet on the set of edges between and including e and e', and calls analyseEdgeLength.
// If e and e' are further than "length" appart, it also checks whether at some stage they become close enough, and if 
// so, the event "close enough" is added to the queue.
void event_activate(Graph* graph, double length, SweepBST* activeEdges, std::priority_queue <Event> &q, int* queue_node_id);

// Function which deals with the deactivate event in the sweep algorithm which finds the l-cuts of the graph
// This function simply removes the edge from the binary search tree of active edges
void event_deactivate(Graph* graph, double length, SweepBST* activeEdges, std::priority_queue <Event> &q);

// Function which deals with the close enough event in the sweep algorithm which finds the l-cuts of the graph
// It only finds all the edges between and including the two edges that have just become close enough,
// and calls the analyseEdgeSet function on this set of edges
void event_close_enough(Graph* graph, double length, SweepBST* activeEdges, std::priority_queue <Event> &q);

// Parent function which finds all the l-cuts of the inputed graph, where length is the length of the disaster
// It works as a sweep algorithm from top to bottom, with three subfunctions doing most of the work
// This function creates activation and deactivation events for all of the edges in the graph and puts them 
// in an even priority queue, organised by height of event. For each event, it calls the relevant function and
// and pops the event. At the end it checks each individual edge to see if it produces an l-cut
void l_cut_finder(Graph* graph, double length);

// Parent function used when adding a new edge to graph. It should be called after l_cut_finder, as it relies
// on the structure graph knowing what the l-cuts are. It finds the best edge based on the search type specified;
// it does some precomputation then calls the relevant sub function. Before doing so, it creates a CGAL arrangement
// (a way to store the graph so that it is easy to find faces, calculate visibility, etc...) of the graph, as well as
// adding a "bounding box" to the graph: Four straight edges which do not intersect the edges of the graph. This is done
// so that the true unbounded face of the graph becomes a bounded face with holes in the arrangement; this is a requirement
// for the visibility functions of CGAL to work. Also if the search type is realted to leaves, it finds the leaf id of 
// each node in G. The search options are the following:
//
//                  - searchType = 'n': Finds the shortest edge between a pair of nodes which are separated
//										by at least one l-cut
//
//                  - searchType = 'h': NEED TO SPECIFY THIS ONCE FINALISED 
//
//                  - searchType = 'l': Finds the shortest edge between a node in a leaf and a node outside
//										the leaf
//
//                  - searchType = 's': Finds L, the leaf containing the lowest indexed node inside a leaf,
//										then finds the shortes edge with an endpoint in L and an endpoint
//										outside of L
//
//					- searchType = 'p': Finds the shortest edge with endpoints in different leaves.
//										Note that this is not only possible even if there are leaves in the graph, 
//										so if such an edge does not exist the search is done like in searchType 'l'
//
//					- searchType = 'c': Combines 'l' and 'p' and chooses the solution based on a weight; for instance
//										it might prefer an edge between leaves if it is no longer than twice an edge 
//										between a leaf and something else
//
//					- searchType = 'x': Custom search, SPECIFY ONCE FINISHED!
//
pair<list<Vertex>, double> find_best_edge(Graph* graph, double length, char searchType, double weight);

bool doIntersect(vector<int> &v1, vector<int> &v2);

pair<list<Vertex>, double> blockSearch(Graph* graph, double length, vector<Vertex*> &sortedVertices, Poly_Arr &graphArr,
	int startType, int endType, double weight);

// Function called when using the searchType 'n'. Requires the graph, the length of the disaster, the vertices of the graph
// (in order to label the arrangement once it is done), the arrangement of the graph and the block id of each node
pair<list<Vertex>, double> node_pair_search(Graph* graph, double length, vector<Vertex*> &sortedVertices, Poly_Arr &graphArr);

// Function called when using the searchType 'h'. Requires the graph, the length of the disaster, the vertices of the graph
// (in order to label the arrangement once it is done), the arrangement of the graph and the block id of each node.
// For each pair of nodes for which an edge exists, a heuristic score is given to the solution based on how many cuts
// are "fixed" as well as well as the edge cost. The heuristic score of the solution is given by
// (1 + K (C - C') / C) * edgeLength, where K is the input weight, C is the number of cuts of G, and C' is the number
// of cuts "fixed" by the new edge 
pair<list<Vertex>, double> node_pair_heuristic_search(Graph* graph, double length, vector<Vertex*> &sortedVertices, Poly_Arr &graphArr, double weight);

// Function called when using the searchType 'l'. Requires the graph, the length of the disaster, the vertices of the graph
// (in order to label the arrangement once it is done), the arrangement of the graph and the leaf id of each node
pair<list<Vertex>, double> leaf_search(Graph* graph, double length, vector<Vertex*> &sortedVertices, Poly_Arr &graphArr);

// Function called when using the searchType 's'. Requires the graph, the length of the disaster, the vertices of the graph
// (in order to label the arrangement once it is done), the arrangement of the graph and the leaf id of each node
pair<list<Vertex>, double> single_leaf_search(Graph* graph, double length, vector<Vertex*> &sortedVertices, Poly_Arr &graphArr);

pair<list<Vertex>, double> random_leaf_search(Graph* graph, double length, vector<Vertex*> &sortedVertices, Poly_Arr &graphArr, int seed);

// Function called when using the searchType 'p'. Requires the graph, the length of the disaster, the vertices of the graph
// (in order to label the arrangement once it is done), the arrangement of the graph and the leaf id of each node
pair<list<Vertex>, double> pair_leaf_search(Graph* graph, double length, vector<Vertex*> &sortedVertices, Poly_Arr &graphArr);

// Function called when using the searchType 'c'. Requires the graph, the length of the disaster, the vertices of the graph
// (in order to label the arrangement once it is done), the arrangement of the graph and the leaf id of each node
// The input value weight specifies how much better the single leaf solution must be over the double leaf solution,
// that is the single leaf solution is picked if costDouble > weight * costSinge
pair<list<Vertex>, double> combined_leaf_search(Graph* graph, double length, vector<Vertex*> &sortedVertices, Poly_Arr &graphArr, double weight);

pair<list<Vertex>, double> random_combined_leaf_search(Graph* graph, double length, vector<Vertex*> &sortedVertices, Poly_Arr &graphArr, double weight, int seed);

// Function called when using the searchType 'x'. Requires the graph, the length of the disaster, the vertices of the graph
// (in order to label the arrangement once it is done), the arrangement of the graph and the leaf id of each node
// The input value weight specifies how much better the single leaf solution must be over the double leaf solution,
// that is the single leaf solution is picked if costDouble > weight * costSinge
pair<list<Vertex>, double> custom_leaf_search(Graph* graph, double length, vector<Vertex*> &sortedVertices, Poly_Arr &graphArr, double weight);

pair<list<Vertex>, double> single_leaf_pair_search(Graph* graph, double length, vector<Vertex*> &sortedVertices, Poly_Arr &graphArr);

pair<list<Vertex>, double> random_leaf_pair_search(Graph* graph, double length, vector<Vertex*> &sortedVertices, Poly_Arr &graphArr, int seed);

double randomAugment(Graph* graph, double length, double weight);

// Parent function which adds edges to graph until it is l-resilient, where l = length. The input searchType specifies
// the search strategy when adding new edges. The function returns the cost of the solution; if the input graph 
// is disconnected the algorithm stops (as it cannot deal with this) and returns NAN. The function calls the functions
// relevant to finding the next best edge in turn: It resets all the helper data, finds the l-cuts of the current graph,
// finds the minimal cuts and leaves, finds the best edge to add and finally it adds it to the graph. It only stops 
// when no edge can be found; if l-cuts remain it returns a cost of NAN.
double augment(Graph* graph, double length, char searchType, double weight);


// Little side function which adds edges to the graph until it is connected.
// It adds the shortest edge which does not intersect any existing edge.
// If the input graph has no edges, the output will be a minimum spanning tree
void connectGraph(Graph* graph);




#endif
