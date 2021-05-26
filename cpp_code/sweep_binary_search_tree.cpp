#ifndef SWEEP_BINARY_SEARCH_TREE_CPP
#define SWEEP_BINARY_SEARCH_TREE_CPP

#include "sweep_binary_search_tree.h"

SweepBST::Node* SweepBST::makeEmpty(Node* t){
	if(t == NULL){return NULL;}
	makeEmpty(t->left);
    makeEmpty(t->right);
    delete t;
    return NULL;
}

SweepBST::Node* SweepBST::findMin(Node* t){
    if(t == NULL){
        return NULL;
    }
    else if(t->left == NULL){
        return t;
    }
    else{
        return findMin(t->left);
    }
}

SweepBST::Node* SweepBST::findMax(Node* t){
	if(t == NULL){
		return NULL;
	}
	else if(t->right == NULL){
		return t;
	}
	else {
		return findMax(t->right);
	}
}

SweepBST::Node* SweepBST::insert(Line* newLine, Node* t, double height){
    // Tree is empty
    if(t == NULL){
        t = new Node;
        t->line = newLine;
        t->left = t->right = NULL;
        return t;
    }
    int comparison = compare(newLine, t->line, height);
    if(comparison == 0){
    	comparison = compare(newLine, t->line, height - 1);
    }
    if(comparison == -1){
        t->left = insert(newLine, t->left, height);
    }else if(comparison == 1){
    	t->right = insert(newLine, t->right, height);
    }else{
    	printf("Shouldn't happen");
    }
    return t;
}

SweepBST::Node* SweepBST::remove(Line* oldLine, Node* t, double height){
    Node* temp;
    if(t == NULL){return NULL;}
    int comparison = compare(oldLine, t->line, height);
    // Check if current node is the one we want to delete
    if(t->line->id == oldLine->id){
    	if(t->left && t->right){
    		temp = findMin(t->right);
    		t->line = temp->line;
    		t->right = remove(t->line, t->right, height);
    	}else{
    		temp = t;
    		if(t->left == NULL){
            	t = t->right;
    		}
        	else if(t->right == NULL){
            	t = t->left;
            }
        	delete temp;
    	}
    	return t;
    }
    if(comparison == 0){
    	comparison = compare(oldLine, t->line, height + 1);
    }
    if(comparison == -1){
        t->left = remove(oldLine, t->left, height);
    }else if(comparison == 1){
        t->right = remove(oldLine, t->right, height);
    }else{
    	printf("Should not happen!\n");
    }

    return t;
}

void SweepBST::inorder(Node* t, int order){
    if(t == NULL){return;}
    inorder(t->left, order + 1);
    cout << "- " << t->line->id << " (" << t->line->start->id << ", " << t->line->end->id << ") (" << order << ") -";
    inorder(t->right, order + 1);
}

void SweepBST::inOrderIdVector(Node* t, vector<int> &idVector, int *current){
	if(t == NULL){return;}
	inOrderIdVector(t->left, idVector, current);
	idVector[*current] = t->line->id;
	*current = *current + 1;
	inOrderIdVector(t->right, idVector, current);
}

SweepBST::Node* SweepBST::find_largest_smaller(Node* root, Line* line, double height){

    Node* temp;
    if(root == NULL){
        return NULL;
    }
    if(root->line->id == line->id){
    	return find_largest_smaller(root->left, line, height);
    }

    int comparison = compare(line, root->line, height);
    if(comparison == 0){
    	comparison = compare(line, root->line, height - 1);
    }
    if(comparison == -1){
    	return find_largest_smaller(root->left, line, height);
    }
    // comparison is 1, need to explore and remember this
    temp =  find_largest_smaller(root->right, line, height);

    if(temp == NULL){return root;}
    return temp;
}

SweepBST::Node* SweepBST::find_smallest_larger(Node* root, Line* line, double height){

    Node* temp;
    if(root == NULL){
        return NULL;
    }
    if(root->line->id == line->id){
    	return find_smallest_larger(root->right, line, height);
    }

    int comparison = compare(line, root->line, height);
    if(comparison == 0){
    	comparison = compare(line, root->line, height - 1);
    }
    if(comparison == 1){
    	return find_smallest_larger(root->right, line, height);
    }
    // comparison is -1, need to explore and remember this
    temp = find_smallest_larger(root->left, line, height);
    if(temp == NULL){return root;}
    return temp;
}

int SweepBST::compare(Line* line1, Line* line2, double height){
	double changeHeight1 = line1->start->y - line1->end->y;
	double changeHeight2 = line2->start->y - line2->end->y;
	double changeWidth1 = line1->start->x - line1->end->x;
	double changeWidth2 = line2->start->x - line2->end->x;
	
	// Start by dealing with horizontal case
	if(std::abs(changeHeight1) < TOLERANCE && std::abs(changeHeight2) < TOLERANCE){
		// Check leftmost and rightmost for both
		double left1 = std::min(line1->start->x,line1->end->x);
		double left2 = std::min(line2->start->x,line2->end->x);
		// Account for sharing endpoints
		if(left1 < left2){return -1;}
		if(left1 > left2){return 1;}
		return 0;
	}else if(std::abs(changeHeight1) < TOLERANCE){
		double x2 = line2->end->x + (height - line2->end->y)*changeWidth2/changeHeight2;
		double left1 = std::min(line1->start->x,line1->end->x);
		// Account for shared endpoints
		if(std::abs(left1 - x2) > TOLERANCE && left1 < x2){
			return -1;
		}else{
			return 1;
		}
	}else if(std::abs(changeHeight2) < TOLERANCE){
		double x1 = line1->end->x + (height - line1->end->y)*changeWidth1/changeHeight1;
		double left2 = std::min(line2->start->x,line2->end->x);
		// Account for shared endpoints
		if(std::abs(left2 - x1) > TOLERANCE && left2 < x1){
			return 1;
		}else{
			return -1;
		}
	}else{
		double x1 = line1->end->x + (height - line1->end->y)*changeWidth1/changeHeight1;
		double x2 = line2->end->x + (height - line2->end->y)*changeWidth2/changeHeight2;
		// Need to deal with shared points; if same, "lower" sweep line and analyse again
		if(std::abs(x1-x2) < TOLERANCE){
			return 0;
		}
		else if(x1 < x2){
			return -1;
		}else{
			return 1;
		}
	}
}

SweepBST::SweepBST(){
	root = NULL;
}

SweepBST::~SweepBST(){
	root = makeEmpty(root);
    }

void SweepBST::insert(Line* line, double height){
    root = insert(line, root, height);
}

void SweepBST::remove(Line* line, double height){
    root = remove(line, root, height);
}

void SweepBST::display()
{
	inorder(root, 0);
    cout << endl;
}

void SweepBST::makeOrderedVector(vector<int> &idVector){
	int current = 0;
	inOrderIdVector(root, idVector, &current);
	idVector[current] = -1;
}

Line* SweepBST::left_neighbour(Line* line, double height){
	Node* temp = find_largest_smaller(root, line, height);
	if(temp){return temp->line;}
	else{return NULL;}
}

Line* SweepBST::right_neighbour(Line* line, double height){
	Node* temp = find_smallest_larger(root, line, height);
	if(temp){return temp->line;}
	else{return NULL;}
}

bool SweepBST::isEmpty(){
    return root == NULL;
}

#endif