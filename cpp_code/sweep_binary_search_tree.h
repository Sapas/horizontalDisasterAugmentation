#ifndef SWEEP_BINARY_SEARCH_TREE_H
#define SWEEP_BINARY_SEARCH_TREE_H

#include <iostream>
#include <stdio.h>
#include <vector>
#include <math.h>
#include "line_geometry.h"
using namespace std;

#define TOLERANCE 0.000001

// Special implementation of a binary search tree 
// Code adapted mainly from (with other additions)
// https://gist.github.com/harish-r/a7df7ce576dda35c9660
class SweepBST{
    struct Node {
        Line* line;
        Node* left;
        Node* right;
    };
    
    Node* root;

    // Clear up memory from tree
    Node* makeEmpty(Node* t);
    // Find smallest node (uses tree geomtry)
    Node* findMin(Node* t);
	// Find largets node (uses tree geomtry)
    Node* findMax(Node* t);
	// Inserts line, sweep specifies y coordinate of sweep line,
	// used to compare segment's positions
    Node* insert(Line* newLine, Node* t, double height);
    // Remove given line
    Node* remove(Line* oldLine, Node* t, double height);
    // Prints all entries in tree in order
    void inorder(Node* t, int order);
    // Puts all entries inside vector in order
    void inOrderIdVector(Node* t, vector<int> &idVector, int *current);
    // Find smaller neighbour, procedure taken and adapted from
    // https://codereview.stackexchange.com/questions/204105/find-largest-smaller-key-in-binary-search-tree
    Node* find_largest_smaller(Node* root, Line* line, double height);
    // Find larger neighbour, adapted from previous link
    Node* find_smallest_larger(Node* root, Line* line, double height);
    // Compares two lines based on state of sweep line
    // 0 = same
    // to the left = -1
    // to the right = 1
    int compare(Line* newLine, Line* line, double height);

public:
    SweepBST();

    ~SweepBST();

    // Wrapper inset function
    void insert(Line* line, double height);

    // Wrapper remove function
    void remove(Line* line, double height);

    // Prints all the nodes of the tree in order
    void display();

    // Puts all entry ids inside vector in order, fills the entry after last with -1
    void makeOrderedVector(vector<int> &idVector);

    // Wrapper left neighbour function
    Line* left_neighbour(Line* line, double height);

    // Wrapper left neighbour function
    Line* right_neighbour(Line* line, double height);

    // Checks if the BST is empty
    bool isEmpty();

};

#endif

