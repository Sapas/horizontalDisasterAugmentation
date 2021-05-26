// These are the sets of functions used in plane geometry,
// stuff such as whether two lines intersect, definition of 
// vertices and graphs, length of edge etc...
// Initially functions were taken from 
// https://www.geeksforgeeks.org/check-if-two-given-line-segments-intersect/

#ifndef LINE_GEOMETRY_H
#define LINE_GEOMETRY_H

#include <iostream>
#include <stdio.h>
#include <vector>
#include <list>
#include <cmath>
#include <float.h>
#include <cstdlib>
#include <algorithm>
using namespace std;

#define TOLERANCE 0.000001


// Start with definitions of vertices, lines and graphs
// to be used as underlying data structures in the code
struct Vertex {
	int id;
	double x, y;
	bool usable;

	bool operator==(const Vertex& otherPoint) const{
		return abs(otherPoint.x - x) < TOLERANCE && abs(otherPoint.y - y) < TOLERANCE;
	}

	bool operator<(const Vertex& otherPoint) const{
		if(abs(otherPoint.y - y) < TOLERANCE){
			if(abs(otherPoint.x - x) < TOLERANCE){
				return false;
			}
			return x < otherPoint.x;
		}
		return y < otherPoint.y;
	}
	Vertex():
	usable(true)
	{}

	Vertex(int id, double x, double y):
	id(id),
	x(x),
	y(y),
	usable(true)
	{}

	Vertex(int id, double x, double y, bool usable):
	id(id),
	x(x),
	y(y),
	usable(usable)
	{}
};

double edgeLength(Vertex* v1, Vertex* v2);

struct Line {
	int id;
	Vertex *start;
	Vertex *end;

	bool operator<(const Line& otherLine) const{
		return edgeLength(start, end) < edgeLength(otherLine.start, otherLine.end);
	}

	Line(int id, Vertex* start, Vertex* end):
	id(id),
	start(start),
	end(end)
	{}
};

struct Graph{
	vector<Vertex*> vertices;
	vector<Line*> edges;
	list< vector<int> > cuts;
	list< vector<int> > minCuts;
	list< vector<int> > cutLabels;
	list< vector<int> > minCutLabels;
	vector< vector<int> > leaves;
	vector< vector<int> > leafMinCut;
	vector< vector<int> > adjacencyMatrix;
	vector<int> blockIds;
	vector< vector<int> > nodeLeaves; 

	Graph(){}
	Graph(int n):
	vertices(n)
	{};

	~Graph(){
		for(int i = 0; i < vertices.size(); i++){
			delete vertices[i];
			vertices[i] = NULL;
		}
		for(int i = 0; i < edges.size(); i++){
			delete edges[i];
			edges[i] = NULL;
		}
		vertices.erase(vertices.begin(), vertices.end());
		edges.erase(edges.begin(), edges.end());

		for(auto it = cuts.begin(); it != cuts.end(); ++it){
			it->erase(it->begin(), it->end());
		}
		cuts.clear();
		
		for(auto it = minCuts.begin(); it != minCuts.end(); ++it){
			it->erase(it->begin(), it->end());
		}
		minCuts.clear();

		for(auto it = cutLabels.begin(); it != cutLabels.end(); ++it){
			it->erase(it->begin(), it->end());
		}
		cutLabels.clear();

		for(auto it = minCutLabels.begin(); it != minCutLabels.end(); ++it){
			it->erase(it->begin(), it->end());
		}
		minCutLabels.clear();

		for(int i = 0; i < leaves.size(); i++){
			leaves[i].erase(leaves[i].begin(), leaves[i].end());
		}
		leaves.erase(leaves.begin(), leaves.end());

		for(int i = 0; i < leafMinCut.size(); i++){
			leafMinCut[i].erase(leafMinCut[i].begin(), leafMinCut[i].end());
		}
		leafMinCut.erase(leafMinCut.begin(), leafMinCut.end());

		for(int i = 0; i < adjacencyMatrix.size(); i++){
			adjacencyMatrix[i].erase(adjacencyMatrix[i].begin(), adjacencyMatrix[i].end());
		}
		adjacencyMatrix.erase(adjacencyMatrix.begin(), adjacencyMatrix.end());

		blockIds.erase(blockIds.begin(), blockIds.end());

		for(int i = 0; i < nodeLeaves.size(); i++){
			nodeLeaves[i].erase(nodeLeaves[i].begin(), nodeLeaves[i].end());
		}
		nodeLeaves.erase(nodeLeaves.begin(), nodeLeaves.end());
	}

};

// Add a segment to the graphs. Since the edge being added
// might have multiple straight lines, might need to define new 
// vertices. It is also possible that overlapping edges are added.
// Need to make sure these are NOT added to the structure of the graph
void addSegment(Graph* graph, list<Vertex> &edge);

// Given a set of sorted vertices and a new point, find the id of the vertex
// on the plane. This has complexity O(log n) on the size of the vertices
int findVertexId(vector<Vertex*> &sortedVertices, Vertex* point);


// Compare two vertices, used as a wrapper for pointer.
// Calls comparison of two vertices which compares their position
bool compareVertices(Vertex* v1, Vertex* v2);

// Given a graph with nodes and edges, and ajacendy matrix is created
void makeAdjacencyMatrix(Graph* graph);

// Wipes all the data from the graph appart from nodes and edges.
void resetGraph(Graph* graph);

// Distance between two lines at a certain height
// Note this function assumes a correct input; that is 
// if lines are horizontal it assumes the height is where the
// line is at
double distanceAtHeight(Line* line1, Line* line2, double height);

// Wrapper function made for pointers for the function below
bool onSegment(Vertex* p, Vertex* q, Vertex* r);

// Given three colinear points p, q, r, the function checks if 
// point q lies on line segment 'pr' 
bool onSegment(Vertex p, Vertex q, Vertex r);
  
// Wrapper function made for pointers for the function below
int orientation(Vertex* p, Vertex* q, Vertex* r);

// To find orientation of ordered triplet (p, q, r). 
// The function returns following values 
// 0 --> p, q and r are colinear 
// 1 --> Clockwise 
// 2 --> Counterclockwise
// See https://www.geeksforgeeks.org/orientation-3-ordered-points/ 
// for details of below formula.  
int orientation(Vertex p, Vertex q, Vertex r);
 
// Wrapper function made for pointers for the function below
bool doIntersect(Vertex* p1, Vertex* q1, Vertex* p2, Vertex* q2);

// The main function that returns true if line segment 'p1q1' 
// and 'p2q2' intersect. 
bool doIntersect(Vertex p1, Vertex q1, Vertex p2, Vertex q2);

// Function which finds at which height two lines have a
// horizontal disatnce of d. Returns NAN if the height is 
// outside of the location of the lnes, or exactly at the endpoints
double when_d_appart(Line* line1, Line* line2, double d);

// Function which finds the length of an edge degined by to vertices
double edgeLength(Vertex* v1, Vertex* v2);

// Function which finds the length of an edge degined by to vertices
double edgeLength(Vertex v1, Vertex v2);

// A function that, given an edge, finds a point collinear to the edge slightly
// further away in the direction of target
Vertex extendEdge(Vertex source, Vertex target);

// A function which finds a point slightly inside the face 
// Note source and target must be different
Vertex moveInsideFace(Vertex source, Vertex middle, Vertex target);

// Helper function, checks if an int is inside a vector
bool intInVect(vector<int> &vector, int id);

#endif