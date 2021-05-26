#ifndef CGAL_INTERFACE_H
#define CGAL_INTERFACE_H

#include <CGAL/Cartesian.h>
#include <CGAL/CORE_algebraic_number_traits.h>
#include <CGAL/Arr_conic_traits_2.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Exact_rational.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/General_polygon_2.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arr_extended_dcel.h>
#include <CGAL/Arr_overlay_2.h>
#include <CGAL/Arr_default_overlay_traits.h>
#include <CGAL/minkowski_sum_2.h>
#include <CGAL/Simple_polygon_visibility_2.h>
#include <CGAL/Triangular_expansion_visibility_2.h>

using namespace std;

#define TOLERANCE 0.000001

// This is a file which deals with the functions provided by CGAL. 
// It is meant to do all of the visibility and minkowski sum stuff

// Start by defining the types that are needed to use CGAL
// Using Exact_predicates_exact_constructions_kernel 
// This Kernel is quite fast and can be used as all values throughout the code
// should be rationals (as there are no disks for now)
typedef CGAL::Exact_predicates_exact_constructions_kernel 							  Rat_Kernel;
typedef Rat_Kernel::Point_2                     									        Rat_Point;
typedef Rat_Kernel::Segment_2                   									        Rat_Segment;
typedef CGAL::Polygon_2<Rat_Kernel>             									        Rat_Polygon;
typedef CGAL::Polygon_with_holes_2<Rat_Kernel>  									        Rat_Polygon_Holes;
typedef CGAL::Arr_segment_traits_2<Rat_Kernel>  									        Poly_Traits;

// Use a "DCEL" which adds data values to the node, edge and face
// classes of an arrangement. Adding a pair of ints to nodes, and int 
// to edges, and a bool to faces. Also define an overlay which deals
// with this addition of types
typedef CGAL::Arr_extended_dcel<Poly_Traits, pair<int,int>, int, bool >   Dcel;
typedef CGAL::Arrangement_2<Poly_Traits, Dcel>        								    Poly_Arr;
typedef CGAL::Arr_face_overlay_traits<Poly_Arr,
                                      Poly_Arr,
                                      Poly_Arr,
                                      logical_and<bool> >                 Overlay_traits;

// Defining a couple of classes that are used in the code
typedef Poly_Arr::Face_handle                              							  Face_Handle;
typedef CGAL::Arrangement_2<Poly_Traits, Dcel >::Ccb_halfedge_circulator	Arr_Edge_Circulator;

// Define the necessary classes for visibility functions
typedef CGAL::Triangular_expansion_visibility_2<Poly_Arr>                 TEV;
typedef CGAL::Arr_naive_point_location<Poly_Arr> 									        Point_Searcher;
typedef CGAL::Arr_point_location_result<Poly_Arr>::Type    							  Point_location_result;

#include "line_geometry.h"
#include "input_output.h"


// It turns out CGAL does not deal well when doing the minkowski sum of two edges
// Luckily this is easy (get either a line or a parallelogram) so doing
// and implementation here. The result is returned as a two vertex polygon
Rat_Polygon edge_minkowski_sum(Line* edge, double length);

// This function takes an edge cut and produces the vulnerable zone based on the length
// of the horizontal disaster. It finds the minkowski sum of every edge in the cut and then
// finds the intersection. The result is returned as a polygon
Rat_Polygon cut_to_vulnerable(Graph* graph, vector<int> &cut, double length);

// This is an extra function needed when the cut includes a horizontal edge (as the minkowski
// sum does not work for two edges in CGAL). It specifies the index of a horizontal edge in 
// cut so the function itself does not have to search for it.
Rat_Polygon cut_to_vulnerable_horizontal(Graph* graph, vector<int> &cut, double length, int flatIndex);

// This function takes the polygon representing a vulnerable zone and 
// turns it into a destruction zone by finding the minkowski sum of the two
Rat_Polygon vulnerable_to_destruction(Rat_Polygon vulnerable_zone, double length);

// This function simply combines the previous two functions cut_to_vulnerable
// and vulnearble_to_destruction to go from an edge cut to a destruction zone
Rat_Polygon cut_to_destruction(Graph* graph, vector<int> &cut, double length);

// This is rather a large function to account for different things. Essentially it works out
// what the arrangement is that a new edge needs to avoid. Note that once the arrangement
// has been worked out, any face appart from the unbounded face represents a destruction zone
// the edge must avoid, so every face which is bounded is labelled "false".
// The input length specifies the size of the disaster, and graph is the input graph.
// What drives the function is the searchType:
//
//                  - searchType = 'l': Adding an edge from a leaf to something else. Label1 specifies
//                                      what the leaf is. The function finds the destruction zone of the 
//                                      minimal cut associated with this leaf and stops
//
//                  - searchType = 'p': Adding an edge between a pair of leaves. Label1 and label2 
//                                      specify the two leaves. The function finds the destruction zones 
//                                      of the minimal cuts associated with both leaves and overlays them
//
//                  - searchType = 'n': Adding an edge between two nodes, specified by label1 and label2.
//                                      Finds the overlay of all the destruction arrangements of the 
//                                      cuts which separate the two nodes
Poly_Arr destructionArrangement(Graph* graph, double length, char searchType, int label1, int label2, vector<pair<int, int> > horizontal_forbidden);

// The way visibility works in CGAL is a polygon is represented by an arrangement face
// This function creates a adjacency matrix for each vertex on the face (and polygon) inputed
// First of all if goes through the face and stores all the nodes on the facial walk in faceVertices
// Then for each vertex v it creates a visibility arrangement from the point of view of v,
// and finds out vertices can be seen from v. If two vertices can see each other, then a straight "edge" is added
// to the adjecency matrix. Note tev is what is used to analse visibility
void makeVisibilityGraph(Graph* graph, Face_Handle &face, TEV &tev, vector<Poly_Arr::Vertex_handle> &faceVertices, vector<list< pair<int,double> > > &adjacentMatrix, vector<pair<int, int> > horizontal_forbidden);

// Defined a custom operator needed for the queue used in this implementation of Disjktra's
class dijkstraCustomComparator { 
public: 
    int operator() (const pair<int, double> &data1, const pair<int, double> &data2){ 
        return data1.second > data2.second;
    } 
};
// This does not use any of the functionality of CGAL. It is an impolementation of Disjktra's shortest path algorithm.
// It is meant to be used after calling the makeVisibilityGraph function. Given an index adjacency list (i.e. the nodes
// are not specified, only their ids on the facial walk) as well as a specification of which nodes are starts (1), goals (2)
// and mids (0), the algorithm finds the shortest path between a start and a goal, and it stores that in path
double dijkstraMultipleSourceGoal(vector<int> &startGoalMid, vector<list< pair<int,double> > > &adjacentMatrix, vector<int> &path);

// This is the main function which finds the shortest edge based on a search type. It takes as inputs the graph and length of disaster, 
// The vertices of the graph sorted by location (in sortedVertices), the arrangement representing the graph (as it is created just once
// so that does not need to be created each time this function is called, as well as label1 and label2 which are used deppending 
// on the value of searchType:
//                
//                  - searchType = 'l': Adding an edge from a leaf to something else. Label1 specifies
//                                      what the leaf is.
//
//                  - searchType = 'p': Adding an edge between a pair of leaves. Label1 and label2 
//                                      specify the two leaves.
//
//                  - searchType = 'n': Adding an edge between two nodes, specified by label1 and label2.
//                                      
pair<list<Vertex>, double> findShortestEdge(Graph* graph, double length, char searchType, int label1, int label2, vector<Vertex*> &sortedVertices, Poly_Arr &graphArr);



#endif