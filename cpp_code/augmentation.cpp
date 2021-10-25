#ifndef AUGMENTAION_CPP
#define AUGMENTAION_CPP

#include "augmentation.h"



int BitsSetTable256[256]; 

void initialize(){  
    // To initially generate the  
    // table algorithmically  
    BitsSetTable256[0] = 0;  
    for (int i = 0; i < 256; i++) {  
        BitsSetTable256[i] = (i & 1) +  
        BitsSetTable256[i / 2];  
    }  
}    

int countSetBits(int n){  
    return (BitsSetTable256[n & 0xff] +  
            BitsSetTable256[(n >> 8) & 0xff] +  
            BitsSetTable256[(n >> 16) & 0xff] +  
            BitsSetTable256[n >> 24]);  
}

int cantorPair(int n1, int n2){
	return (n1 + n2)*(n1 + n2 + 1)/2 + n2;
}

void vectorIntEncoder(vector<int> &currentVector, vector<int> &addedVector){
	// First find max of both vectors
	int max1 = 0;
	int max2 = 0;
	for(int i = 0; i < currentVector.size(); i++){
		if(currentVector[i] > max1){max1 = currentVector[i];}
		if(addedVector[i] > max2){max2 = addedVector[i];}
	}
	// Get max possible id
	int maxId = cantorPair(max1, max2);
	// Get vector with relabelling, assign size
	vector<int> newIds(maxId + 1,-1);
	int currentId = 0;
	// Combine vectors and assign new label
	for(int i = 0; i < currentVector.size(); i++){
		int midId = cantorPair(currentVector[i], addedVector[i]);
		if(newIds[midId] == -1){
			// Have not assigned id to this yet
			newIds[midId] = currentId;
			currentId++;
		}
		currentVector[i] = newIds[midId];
	}
}

bool isConnected(Graph* graph){
	vector<int> visited(graph->vertices.size(), 0);
	vector< pair <int,int> > forbidden(0);
	int group = 1;
	for(int i = 0; i < graph->vertices.size(); i++){
		if(visited[i] != 0){continue;}
		depthFirstSearch(graph, forbidden, visited, i, group);
		group++;
	}
	// if(group != 2){
	// 	for(int i = 0; i < visited.size(); i++){
	// 		printf("%d ", visited[i]);
	// 	}
	// 	printf("\n");
	// 	if(group == 1){
	// 		printf("EMPTY!!\n");
	// 	}
	// }
	
	return group == 2;
}

bool isSubset(vector<int> &smallVector, vector<int> &bigVector){
	int extra = 0, i;
	for(i = 0; i < smallVector.size(); i++){
		while(i + extra < bigVector.size() && smallVector[i] != bigVector[i+extra]){extra++;}
		if(i + extra >= (int)bigVector.size()){break;}
	}
	return i == (int)smallVector.size();
}

void checkPotentialLeaf(Graph* graph, vector<int> &potentialLeaf, vector<int> &minCut){
	bool isLeaf = true;
	bool added = false;
	for(int i = 0; i < graph->leaves.size() && isLeaf; i++){
		// Start with leaf 1
		if(isLeaf && (int)potentialLeaf.size() < (int)graph->leaves[i].size() && isSubset(potentialLeaf, graph->leaves[i])){
			// Clean up vector
			graph->leaves[i].erase(graph->leaves[i].begin(), graph->leaves[i].end());
			graph->leafMinCut[i].erase(graph->leafMinCut[i].begin(), graph->leafMinCut[i].end());
			// If have removed at least one leaf already, so only need to remove current leaf of graph
			if(added){
				graph->leaves.erase(graph->leaves.begin()+i);
				graph->leafMinCut.erase(graph->leafMinCut.begin()+i);
			}
			// First time removing, so can simply replace
			else{
				added = true;
				graph->leaves[i] = potentialLeaf;
				graph->leafMinCut[i] = minCut;
			}

		}else if(isLeaf && (int)potentialLeaf.size() >= (int)graph->leaves[i].size() && isSubset(graph->leaves[i], potentialLeaf)){
			// A subset exists, so potential leaf 1 is not a leaf. No need to keep checking as it cannot
			// be a subset of anything in graph leaves (otherwise the current graph leaf would also be a subset)
			isLeaf = false;
		}
	}
	// If it is leaf and have not added yet, do so now
	if(isLeaf && !added){
		graph->leaves.push_back(potentialLeaf);
		graph->leafMinCut.push_back(minCut);
	}
}

void storeNodeLeaves(Graph* graph){
	// Store the leaves of each node. Note there might be more than one!
	// Start by allocating memory
	vector<vector<int> > nodeLeaves;
	nodeLeaves.reserve(graph->vertices.size());
	for(int i = 0; i < graph->vertices.size(); i++){
		vector<int> localIds;
		localIds.reserve(graph->leaves.size());
		nodeLeaves.push_back(localIds);
	}
	// Go through the leaves and "add" them to this
	for(int leaf = 0; leaf < graph->leaves.size(); leaf++){
		bool leafHasUsableNode = false;
		for(int nodeIndex = 0; nodeIndex < graph->leaves[leaf].size(); nodeIndex++){
			// First check the node is usable
			if(!(graph->vertices[graph->leaves[leaf][nodeIndex]]->usable)){continue;}
			leafHasUsableNode = true;
			nodeLeaves[graph->leaves[leaf][nodeIndex]].push_back(leaf);
		}
		if(!leafHasUsableNode){
			printf("Leaf %d has only unusable nodes, this should not happen!! Nodes are ", leaf);
			for(int nodeIndex = 0; nodeIndex < graph->leaves[leaf].size(); nodeIndex++){
				printf("%d ", graph->leaves[leaf][nodeIndex]);
			}
			printf("\n");
		}
	}
	graph->nodeLeaves = nodeLeaves;
}

void findLeaves(Graph* graph){
	// Cycle through the min cuts and their labels, then call another function
	// which looks at whether the potential leaf is a leaf (at least currently)
	auto cutIt = graph->minCuts.begin();
	auto labelIt = graph->minCutLabels.begin();
	int count = 0;
	while(cutIt != graph->minCuts.end() && labelIt != graph->minCutLabels.end()){
		count++;
		// printf("Working min cut %d/%d (%.5f%%)\n", count, (int)graph->minCuts.size(), double(100.0 * count / graph->minCuts.size()));
		vector<int> potentialLeaf1;
		vector<int> potentialLeaf2;
		potentialLeaf1.reserve((int)graph->vertices.size());
		potentialLeaf2.reserve((int)graph->vertices.size());

		for(int i = 0; i < labelIt->size(); i++){
			if((*labelIt)[i] == 1){
				potentialLeaf1.push_back(i);
			}else if((*labelIt)[i] == 2){
				potentialLeaf2.push_back(i);
			}else{
				printf("Wrong label in leaf finding, got %d\n",(*labelIt)[i]);
			}
		}
		checkPotentialLeaf(graph, potentialLeaf1, *cutIt);
		checkPotentialLeaf(graph, potentialLeaf2, *cutIt);
		++cutIt;
		++labelIt;
	}
	storeNodeLeaves(graph);
}

void depthFirstSearch(Graph* graph, vector< pair <int,int> > &forbidden, vector<int> &visited, int vertexId, int group){
	visited[vertexId] = group;
	for(int i = 0; i < graph->adjacencyMatrix[vertexId].size(); i++){
		bool ok = true;
		for(int j = 0; j < forbidden.size(); j++){
			if((vertexId == forbidden[j].first && 
				graph->adjacencyMatrix[vertexId][i] == forbidden[j].second) ||
				(vertexId == forbidden[j].second && 
				graph->adjacencyMatrix[vertexId][i] == forbidden[j].first)){
				
				ok = false;
				break;
			}
		}
		if(ok && visited[graph->adjacencyMatrix[vertexId][i]] == 0){
			depthFirstSearch(graph, forbidden, visited, graph->adjacencyMatrix[vertexId][i], group);
		}
	}
}

bool checkBipartite(vector<vector<int> > &adjacencyMatrix, vector<int> &group, int vertexId, int groupId){
	group[vertexId] = groupId;
	int otherId = 1;
	if(groupId == 1){otherId = 2;}
	for(int i = 0; i < adjacencyMatrix[vertexId].size(); i++){
		int neighbour = adjacencyMatrix[vertexId][i];
		if(group[neighbour] == groupId){
			//Looking at node with same id, false
			return false;
		}
		if(group[neighbour] != groupId && 
			group[neighbour] != 0){
			//Looking at node with different id, continue
			continue;
		}
		if(!checkBipartite(adjacencyMatrix, group, neighbour, otherId)){
			//Previous function was false, so return false
			return false;
		}
	}
	//All good return true!
	return true;
}

void checkSingleEdgesForLeaves(Graph* graph){
	for(int i = 0; i < graph->edges.size(); i++){
		vector<int> cut;
		cut.push_back(graph->edges[i]->id);
		vector< pair <int,int> > forbidden;
		forbidden.push_back(make_pair(graph->edges[i]->start->id,graph->edges[i]->end->id));
		// Simply check if graph is disconnected, if so this is a cut!
		vector<int> visited(graph->vertices.size(), 0);
		int group = 1;
		for(int i = 0; i < graph->vertices.size(); i++){
			if(visited[i] != 0){continue;}
			depthFirstSearch(graph, forbidden, visited, i, group);
			group++;
		}
		// If the graph is connected; if it is not an l-cut!
		if(group == 2){continue;}

		// Since it is not connected, either side could be a leaf, process it
		// First store potential leaf (i.e. all nodes in the block)
		vector<int> potentialLeaf1;
		vector<int> potentialLeaf2;
		potentialLeaf1.reserve(graph->vertices.size());
		potentialLeaf2.reserve(graph->vertices.size());
		for(int i = 0; i < graph->vertices.size(); i++){
			if(visited[i] == 1){
				potentialLeaf1.push_back(i);
			}else if(visited[i] == 2){
				potentialLeaf2.push_back(i);
			}else{
				printf("Something weird going on, more than 2 blocks when removing a single edge!\n");
			}
		}
		// Now store potential minimal l-cut (which is a single edge)
		vector<int> potentialLeafCut;
		potentialLeafCut.push_back(i);
		// Have two potential leaves, process them and store them if it could be!
		checkPotentialLeaf(graph, potentialLeaf1, potentialLeafCut);
		checkPotentialLeaf(graph, potentialLeaf2, potentialLeafCut);
	}
}


void checkSingleEdgesCut(Graph* graph, bool onlyCheckForLeaves){
	if(onlyCheckForLeaves){return checkSingleEdgesForLeaves(graph);}
	for(int i = 0; i < graph->edges.size(); i++){
		vector<int> cut;
		cut.push_back(graph->edges[i]->id);
		vector< pair <int,int> > forbidden;
		forbidden.push_back(make_pair(graph->edges[i]->start->id,graph->edges[i]->end->id));
		// Simply check if graph is disconnected, if so this is a cut!
		vector<int> visited(graph->vertices.size(), 0);
		int group = 1;
		for(int i = 0; i < graph->vertices.size(); i++){
			if(visited[i] != 0){continue;}
			depthFirstSearch(graph, forbidden, visited, i, group);
			group++;
		}
		// If the graph is connected; if it is not an l-cut!
		if(group == 2){continue;}
		// graph->cuts.push_back(cut);
		graph->minCuts.push_back(cut);
		graph->minCutLabels.push_back(visited);
		// Need to update the blocks as you go, but doing it for minimal l-cuts is sufficient
		vectorIntEncoder(graph->blockIds, visited);

	}
}

void findLabelsAndMinCuts(Graph* graph){
	auto cutIt = graph->cuts.begin();
	// Initialise the block each node is in
	vector<int> ids(graph->vertices.size(), 1);
	graph->blockIds = ids;
	int count = 0; 
	while(cutIt != graph->cuts.end()){
		count++;
		printf("Working cut %d/%d (%.5f%%)\n", count, (int)graph->cuts.size(), (double)(100.0 * count / graph->cuts.size()));

		vector<pair<int,int> > forbidden;
		forbidden.reserve(cutIt->size());
		for(int i = 0; i < cutIt->size(); i++){
			forbidden.push_back(make_pair(graph->edges[(*cutIt)[i]]->start->id, graph->edges[(*cutIt)[i]]->end->id));	
		}
		// DFS to find the blocks
		vector<int> visited(graph->vertices.size(), 0);
		int group = 1;

		for(int i = 0; i < graph->vertices.size(); i++){
			if(visited[i] != 0){continue;}
			depthFirstSearch(graph, forbidden, visited, i, group);
			group++;
		}
		// Sweet, have the groups! Can save

		// I actually could comment this out as for now this is never used and it might take a lot of memory to store
		graph->cutLabels.push_back(visited); 

		// If only have 3 groups (so 2 really) this is a minimal cut
		if(group == 3){
			graph->minCuts.push_back((*cutIt));
			graph->minCutLabels.push_back(visited);
			// Need to update the blocks as you go, but doing it for minimal l-cuts is sufficient
			vectorIntEncoder(graph->blockIds, visited);
		}
		
		++cutIt;
	}
}


void analyseEdgeSetForLeaves(Graph* graph, vector<int> &vectorEdgesId, int left, int right){
	// This should be very similar to the analyseEdgeSet function
	// First extract edges that have been destroyed
	vector<int> edgeSet(right - left + 1);
	vector<pair<int,int> > forbidden(right - left + 1);
	for(int i = left; i <= right; i++){
		edgeSet[i - left] = vectorEdgesId[i];
		forbidden[i - left] = make_pair(graph->edges[vectorEdgesId[i]]->start->id, graph->edges[vectorEdgesId[i]]->end->id);
	}
	// Ok first thing I want to do, is to do a depth first search to find the blocks after removing these edges
	vector<int> visited(graph->vertices.size(), 0);
	int group = 1;
	for(int i = 0; i < graph->vertices.size(); i++){
		if(visited[i] != 0){continue;}
		depthFirstSearch(graph, forbidden, visited, i, group);
		group++;
	}
	// Check if the graph is connected; if it is no l-cuts here 
	bool connected = true;
	for(int i = 0; i < graph->vertices.size(); i++){
		if(visited[i] != 1){
			connected = false;
			break;
		}
	}
	if(connected){return;}
	// Ok so now have blocks, all of which are potential leaves. 
	// For each block B check G_{V\B} is connected
	for(int block = 1; block < group; block++){
		// First need to check that shores of cut associated 
		// with this block are connected, i.e. check the cut is minimal
		vector<int> potentialLeafCut;
		potentialLeafCut.reserve(right - left + 1);
		vector<pair<int,int> > localForbidden;
		localForbidden.reserve(right - left + 1);
		// Store l-cut at the same time
		for(int i = left; i <= right; i++){
			int startNode = graph->edges[vectorEdgesId[i]]->start->id;
			int endNode = graph->edges[vectorEdgesId[i]]->end->id;

			if(visited[startNode] == visited[endNode] ||
				(visited[startNode] != block &&
					visited[endNode] != block)){
				continue;
			}
			// This is in the cut! Store it
			potentialLeafCut.push_back(vectorEdgesId[i]);
			localForbidden.push_back(make_pair(startNode, endNode));
		}
		// Now run depth first search to see whether shores are connected
		vector<int> localVisited(graph->vertices.size(), 0);
		int localGroup = 1;
		for(int i = 0; i < graph->vertices.size(); i++){
			if(localVisited[i] != 0){continue;}
			depthFirstSearch(graph, localForbidden, localVisited, i, localGroup);
			localGroup++;
		}
		// If more than 2 groups, shores not connected
		if(localGroup > 3){continue;}

		// Now can store potential leaf (i.e. all nodes in the block)
		vector<int> potentialLeaf;
		potentialLeaf.reserve(graph->vertices.size());
		for(int i = 0; i < graph->vertices.size(); i++){
			if(visited[i] == block){
				potentialLeaf.push_back(i);
			}
		}

		// Have a potential leaf, process it and store it if it could be!
		checkPotentialLeaf(graph, potentialLeaf, potentialLeafCut);
	}

}




void analyseEdgeSet(Graph* graph, vector<int> &vectorEdgesId, int left, int right, int middle, bool onlyCheckForLeaves){
	if(onlyCheckForLeaves){return analyseEdgeSetForLeaves(graph, vectorEdgesId, left, right);}
	// First of all need to extract the edges from the active edges vector
	vector<int> edgeSet(right - left + 1);
	vector<pair<int,int> > forbidden(right - left + 1);
	for(int i = left; i <= right; i++){
		edgeSet[i - left] = vectorEdgesId[i];
		forbidden[i - left] = make_pair(graph->edges[vectorEdgesId[i]]->start->id, graph->edges[vectorEdgesId[i]]->end->id);
	}
	int newMiddle = middle - left;
	// Ok first thing I want to do, is to do a depth first search to find the blocks after removing these edges
	vector<int> visited(graph->vertices.size(), 0);
	int group = 1;
	for(int i = 0; i < graph->vertices.size(); i++){
		if(visited[i] != 0){continue;}
		depthFirstSearch(graph, forbidden, visited, i, group);
		group++;
	}
	// Check if the graph is connected; if it is no l-cuts here 
	bool connected = true;
	for(int i = 0; i < graph->vertices.size(); i++){
		if(visited[i] != 1){
			connected = false;
			break;
		}
	}
	if(connected){return;}
	// Ok so now have the blocks. First need to check left, right and middle (if applicable)
	// can be part of a cut. The only reasons this might not be the case is if either of the three
	// have the same endpoints in the same block, or if the three form a cycle of size three.
	// First check none of these edges have both endpoints in the same block
	if(visited[forbidden[0].first] == visited[forbidden[0].second]){return;}
	if(visited[forbidden.back().first] == visited[forbidden.back().second]){return;}
	if(middle != -1 && visited[forbidden[newMiddle].first] == visited[forbidden[newMiddle].second]){return;}
	// First create the divisions related to the left, right and middle edges. At most 4 such divisions
	vector< vector<int> > shoresA(4), shoresB(4);
	shoresA.reserve(4);
	shoresB.reserve(4);
	// Start with left
	vector<int> shoreA(1, visited[forbidden[0].first]);
	vector<int> shoreB(1, visited[forbidden[0].second]);
	shoresA[0] = shoreA;
	shoresB[0] = shoreB;
	// Ok now want to add in right. Might not add anything, will find out
	// Just going to spell out all options and then condense
	bool firstInA = intInVect(shoresA[0], visited[forbidden.back().first]);
	bool firstInB = intInVect(shoresB[0], visited[forbidden.back().first]);
	bool secondInA = intInVect(shoresA[0], visited[forbidden.back().second]);
	bool secondInB = intInVect(shoresB[0], visited[forbidden.back().second]);
	
	if(firstInA && secondInB){/*nothing*/}
	else if(firstInB && secondInA){/*nothing*/}
	else if(firstInA && !secondInB){shoresB[0].push_back(visited[forbidden.back().second]);}
	else if(secondInA && !firstInB){shoresB[0].push_back(visited[forbidden.back().first]);}
	else if(!firstInA && secondInB){shoresA[0].push_back(visited[forbidden.back().first]);}
	else if(!secondInA && firstInB){shoresA[0].push_back(visited[forbidden.back().second]);}
	else{
		shoresA[1] = shoresA[0];
		shoresB[1] = shoresB[0];
		shoresA[0].push_back(visited[forbidden.back().first]);
		shoresB[0].push_back(visited[forbidden.back().second]);
		shoresA[1].push_back(visited[forbidden.back().second]);
		shoresB[1].push_back(visited[forbidden.back().first]);
	}
	// Perfect, time to add in middle if needed
	// If middle not used, can skip
	if(middle == -1){/*nothing*/}
	else{
		// Again spell out options and maybe condensate later
		bool firstInA = intInVect(shoresA[0], visited[forbidden[newMiddle].first]);
		bool firstInB = intInVect(shoresB[0], visited[forbidden[newMiddle].first]);
		bool secondInA = intInVect(shoresA[0], visited[forbidden[newMiddle].second]);
		bool secondInB = intInVect(shoresB[0], visited[forbidden[newMiddle].second]);

		// Need to check a potential cycle does not exist
		if((firstInA && secondInA) || (firstInB && secondInB)){shoresA[0].clear(); shoresB[0].clear();}
		else if(firstInA && secondInB){/*nothing*/}
		else if(firstInB && secondInA){/*nothing*/}
		else if(firstInA && !secondInB){shoresB[0].push_back(visited[forbidden[newMiddle].second]);}
		else if(secondInA && !firstInB){shoresB[0].push_back(visited[forbidden[newMiddle].first]);}
		else if(!firstInA && secondInB){shoresA[0].push_back(visited[forbidden[newMiddle].first]);}
		else if(!secondInA && firstInB){shoresA[0].push_back(visited[forbidden[newMiddle].second]);}
		else{
			// Want new row
			shoresA[2] = shoresA[0];
			shoresB[2] = shoresB[0];
			shoresA[0].push_back(visited[forbidden.back().first]);
			shoresB[0].push_back(visited[forbidden.back().second]);
			shoresA[2].push_back(visited[forbidden.back().second]);
			shoresB[2].push_back(visited[forbidden.back().first]);
		}

		if(shoresA[1].size() > 0){
			// Again spell out options and maybe condensate later
			bool firstInA = intInVect(shoresA[1], visited[forbidden[newMiddle].first]);
			bool firstInB = intInVect(shoresB[1], visited[forbidden[newMiddle].first]);
			bool secondInA = intInVect(shoresA[1], visited[forbidden[newMiddle].second]);
			bool secondInB = intInVect(shoresB[1], visited[forbidden[newMiddle].second]);

			// Need to check a potential cycle does not exist
			if((firstInA && secondInA) || (firstInB && secondInB)){shoresA[1].clear(); shoresB[1].clear();}
			else if(firstInA && secondInB){/*nothing*/}
			else if(firstInB && secondInA){/*nothing*/}
			else if(firstInA && !secondInB){shoresB[1].push_back(visited[forbidden[newMiddle].second]);}
			else if(secondInA && !firstInB){shoresB[1].push_back(visited[forbidden[newMiddle].first]);}
			else if(!firstInA && secondInB){shoresA[1].push_back(visited[forbidden[newMiddle].first]);}
			else if(!secondInA && firstInB){shoresA[1].push_back(visited[forbidden[newMiddle].second]);}
			else{
				// Want new row
				shoresA[3] = shoresA[1];
				shoresB[3] = shoresB[1];
				shoresA[1].push_back(visited[forbidden.back().first]);
				shoresB[1].push_back(visited[forbidden.back().second]);
				shoresA[3].push_back(visited[forbidden.back().second]);
				shoresB[3].push_back(visited[forbidden.back().first]);
			}
		}

	}
	// Sweet, so now have the "back bones" of the division of the blocks. 
	// Now want to go through the remaining blocks and assign them.
	for(int j = 0; j < shoresA.size(); j++){
		if(shoresA[j].size() == 0){continue;}
		int options = group - 1 - shoresA[j].size() - shoresB[j].size();
		int n_subsets = pow(2,options);
		for(int subset = 0; subset < n_subsets; subset++){
			vector<int> division(group, 0);
			// Assign already decided groups
			for(int k = 0; k < shoresA[j].size(); k++){
				division[shoresA[j][k]] = 1;
			}
			for(int k = 0; k < shoresB[j].size(); k++){
				division[shoresB[j][k]] = 2;
			}
			int offset = 1;
			for(int digit = 1; digit < group; digit++){
				if(division[digit] != 0){offset++; continue;}
				if(((subset >> (digit - offset)) & 1) == 1){division[digit] = 1;}
				else{division[digit] = 2;}
			}
			// Finally, cycle through each edge and check whether it is part of this partition!
			vector<int> newCut;
			vector<pair<int,int> > newCutPoints;
			newCut.reserve(right - left + 1);
			for(int k = 0; k < edgeSet.size(); k++){
				if(division[visited[forbidden[k].first]] != division[visited[forbidden[k].second]]){
					newCut.push_back(edgeSet[k]);
				}
			}
			// Got new cut, should store it
			// graph->cuts.push_back(newCut);
			// Instead, check if this is a minimal l-cut and store it if it is
			vector<pair<int,int> > newForbidden((int)newCut.size());
			for(int i = 0; i < (int)newCut.size(); i++){
				newForbidden[i] = make_pair(graph->edges[newCut[i]]->start->id, graph->edges[newCut[i]]->end->id);
			}
			vector<int> newVisited(graph->vertices.size(), 0);
			int newGroup = 1;
			for(int i = 0; i < graph->vertices.size(); i++){
				if(newVisited[i] != 0){continue;}
				depthFirstSearch(graph, newForbidden, newVisited, i, newGroup);
				newGroup++;
			}
			if(newGroup == 3){
				graph->minCuts.push_back(newCut);
				graph->minCutLabels.push_back(newVisited);
				// Need to update the blocks as you go, but doing it for minimal l-cuts is sufficient
				vectorIntEncoder(graph->blockIds, newVisited);
			}


		}
	}
}

void event_activate(Graph* graph, double length, SweepBST* activeEdges, std::priority_queue <Event> &q, int* queue_node_id, bool onlyCheckForLeaves){
	// printf("Height %.2f edge %d (%d, %d) start\n", q.top().height, q.top().line1->id, q.top().line1->start->id, q.top().line1->end->id);
	std::vector<int> vectorEdgesId(graph->edges.size() + 1);
	activeEdges->insert(q.top().line1, q.top().height);
	activeEdges->makeOrderedVector(vectorEdgesId);
	int middle, left, right;
	//Find position of current edge
	middle = 0;
	while(vectorEdgesId[middle] != q.top().line1->id){middle++;}
	left = middle - 1;
	right = middle + 1;
	// For existing cuts, need to add this edge if in between close enough edges. Note must be strictly less than disaster,
	// If exactly the same as the disaster, it will be added in the even "close enough"
	while(left >= 0){
		right = middle + 1;
		while(vectorEdgesId[right] != -1){
			double dist = distanceAtHeight(graph->edges[vectorEdgesId[left]], graph->edges[vectorEdgesId[right]], q.top().height);
			if(abs(dist - length) < TOLERANCE || dist > length){break;}
			analyseEdgeSet(graph, vectorEdgesId, left, right, middle, onlyCheckForLeaves);
			right++;
		}
		left--;
	}
	// Will need a check on whether the current line is horizontal
	bool isHorizontal = std::abs(q.top().line1->start->y - q.top().line1->end->y) < TOLERANCE;
	// For every other edge, if close enough need to register cuts, 
	// otherwise need to find when and if will be close enough and register event
	// Start with left
	left = middle - 1;
	while(left >= 0){
		// If close enough, add subsets now
		double dist = distanceAtHeight(graph->edges[vectorEdgesId[left]], graph->edges[vectorEdgesId[middle]], q.top().height);
		if(dist < length && abs(dist - length) > TOLERANCE){
			analyseEdgeSet(graph, vectorEdgesId, left, middle, -1, onlyCheckForLeaves);
		}

		// Otherwise, check if eventually will be close enough. Only do this if both lines are not horizontal
		// Note that due to order of events, cannot activate a non horizontal line while a horizontal line is active
		// So only need to check if "middle" line is horizontal
		else if(isHorizontal){break;}
		else{
			double height = when_d_appart(graph->edges[vectorEdgesId[left]], graph->edges[vectorEdgesId[middle]], length);
			if(abs(dist - length) < TOLERANCE){
				//Need to check if growing closer or further apart. If staying the same ignore too!
				double futureDist = distanceAtHeight(graph->edges[vectorEdgesId[left]], graph->edges[vectorEdgesId[middle]], q.top().height - 1);
				if(futureDist < length && abs(futureDist - length) > TOLERANCE){
					q.push(Event(*queue_node_id, 'c', q.top().height, graph->edges[vectorEdgesId[left]], graph->edges[vectorEdgesId[middle]]));
					*queue_node_id = *queue_node_id + 1;
				}
			}else if(!isnan(height) && height < q.top().height){
				q.push(Event(*queue_node_id, 'c', height, graph->edges[vectorEdgesId[left]], graph->edges[vectorEdgesId[middle]]));
				*queue_node_id = *queue_node_id + 1;
			}
		}
			
		left--;
	}
	// Repeat procedure with right
	right = middle + 1;
	while(vectorEdgesId[right] != -1){
		double dist = distanceAtHeight(graph->edges[vectorEdgesId[middle]], graph->edges[vectorEdgesId[right]], q.top().height);
		if(dist < length && abs(dist - length) > TOLERANCE){
			analyseEdgeSet(graph, vectorEdgesId, middle, right, -1, onlyCheckForLeaves);
		}
		else if(isHorizontal){break;}
		else{
			double height = when_d_appart(graph->edges[vectorEdgesId[middle]], graph->edges[vectorEdgesId[right]], length);
			if(abs(dist - length) < TOLERANCE){
				//Need to check if growing closer or further apart. If staying the same ignore too!
				double futureDist = distanceAtHeight(graph->edges[vectorEdgesId[middle]], graph->edges[vectorEdgesId[right]], q.top().height - 1);
				if(futureDist < length && abs(futureDist - length) > TOLERANCE){
					q.push(Event(*queue_node_id, 'c', q.top().height, graph->edges[vectorEdgesId[middle]], graph->edges[vectorEdgesId[right]]));
					*queue_node_id = *queue_node_id + 1;
				}
			}else if(!isnan(height) && height < q.top().height){
				q.push(Event(*queue_node_id, 'c', height, graph->edges[vectorEdgesId[middle]], graph->edges[vectorEdgesId[right]]));
				*queue_node_id = *queue_node_id + 1;
			}
		}
		right++;
	}
}

void event_deactivate(Graph* graph, double length, SweepBST* activeEdges, std::priority_queue <Event> &q){
	// printf("Height %.2f edge %d (%d, %d) end\n", q.top().height, q.top().line1->id, q.top().line1->start->id, q.top().line1->end->id);
	activeEdges->remove(q.top().line1, q.top().height);
}

void event_close_enough(Graph* graph, double length, SweepBST* activeEdges, std::priority_queue <Event> &q, bool onlyCheckForLeaves){
	// printf("Height %.2f edges %d (%d, %d) and %d (%d, %d) close enough\n", q.top().height,
		// q.top().line1->id, q.top().line1->start->id, q.top().line1->end->id,
		// q.top().line2->id, q.top().line2->start->id, q.top().line2->end->id);
	std::vector<int> vectorEdgesId(graph->edges.size() + 1);
	activeEdges->makeOrderedVector(vectorEdgesId);
	int left, right;
	left = 0;
	while(vectorEdgesId[left] != q.top().line1->id){left++;}
	right = left+1;
	while(vectorEdgesId[right] != q.top().line2->id){right++;}
	analyseEdgeSet(graph, vectorEdgesId, left, right, -1, onlyCheckForLeaves);
}

void l_cut_finder(Graph* graph, double length, bool onlyCheckForLeaves){
	// Function to find cuts of the graph
	// Start by using the binary search tree and the queue to
	// make sure all is well.
	SweepBST activeEdges;
	// If knew how, would reserve space 2*|E| + |E|^2
	std::priority_queue <Event> q;
	int queue_node_id = 0;
	std::vector<int> vectorEdgesId(graph->edges.size() + 1);

	// Put in start and end of each edge
	for(int i = 0; i < graph->edges.size(); i++){
		double start = std::max(graph->edges[i]->start->y, graph->edges[i]->end->y);
		double end = std::min(graph->edges[i]->start->y, graph->edges[i]->end->y);
		q.push(Event(queue_node_id, 'a', start, graph->edges[i], NULL));
		queue_node_id++;
		q.push(Event(queue_node_id, 'd', end, graph->edges[i], NULL));
		queue_node_id++;
	}
	while(!q.empty()){

		if(q.top().type == 'a'){event_activate(graph, length, &activeEdges, q, &queue_node_id, onlyCheckForLeaves);}
		else if(q.top().type == 'd'){event_deactivate(graph, length, &activeEdges, q);}
		else if(q.top().type == 'c'){event_close_enough(graph, length, &activeEdges, q, onlyCheckForLeaves);}
		else{printf("Weird behaviour, event is not a, d or c!\n");}		
		q.pop();
	}
	if(!activeEdges.isEmpty()){
		printf("Active edges in queue, something went wrong:\n");
		activeEdges.display();
	}
	// Also need to look at single edges (they might be cuts as well)
	checkSingleEdgesCut(graph, onlyCheckForLeaves);
}

pair<list<Vertex>, double> find_best_edge(Graph* graph, double length, char searchType, double weight){
	// Will need the vector of vertices sorted for reference
	vector<Vertex*> sortedVertices = graph->vertices;
	sort(sortedVertices.begin(), sortedVertices.end(), compareVertices);

	// Make an arrangement out of the graph for visibility purposes
	Poly_Arr graphArr;
	// Add the edges to the arrangement
	for(int i = 0; i < graph->edges.size(); i++){
		Rat_Segment newEdge = Rat_Segment(Rat_Point(graph->edges[i]->start->x, graph->edges[i]->start->y),
											Rat_Point(graph->edges[i]->end->x, graph->edges[i]->end->y));
		insert(graphArr, newEdge);
	}
	// Will need to add "boundary" box for outside face, will need min and max x,y coordinates
	double x_min = DBL_MAX, y_min = DBL_MIN, x_max = -DBL_MIN, y_max = -DBL_MIN;
	for(int i = 0; i < graph->vertices.size(); i++){
		if(graph->vertices[i]->x > x_max){x_max = graph->vertices[i]->x;}
		if(graph->vertices[i]->x < x_min){x_min = graph->vertices[i]->x;}
		if(graph->vertices[i]->y > y_max){y_max = graph->vertices[i]->y;}
		if(graph->vertices[i]->y < y_min){x_max = graph->vertices[i]->y;}
	}
	// Add four "box" edges
	insert(graphArr, Rat_Segment(Rat_Point(x_min - 2 * length, y_min - 2 * length),
														Rat_Point(x_max + 2 * length, y_min - 2 * length)));
	insert(graphArr, Rat_Segment(Rat_Point(x_max + 2 * length, y_min - 2 * length),
														Rat_Point(x_max + 2 * length, y_max + 2 * length)));
	insert(graphArr, Rat_Segment(Rat_Point(x_max + 2 * length, y_max + 2 * length),
														Rat_Point(x_min - 2 * length, y_max + 2 * length)));
	insert(graphArr, Rat_Segment(Rat_Point(x_min - 2 * length, y_max + 2 * length),
														Rat_Point(x_min - 2 * length, y_min - 2 * length)));
	// Go through each face and label it as worth analysing (i.e. might be able to add an edge here)
	// Exclude unbounded (i.e. outside of box).
	for (auto fit = graphArr.faces_begin(); fit != graphArr.faces_end(); ++fit){fit->set_data(!(fit->is_unbounded()));}
	// Here is where the search choice happens
	// Just making search types a through to r, will refine and give proper names when finalising a few

	// Modifying this so that can free memory of graphArr
	pair<list<Vertex>, double> bestSolution;
	bestSolution.second = DBL_MAX;
	if(searchType == 'a'){bestSolution = blockSearch(graph, length, sortedVertices, graphArr, 0, 0, 0);}
	else if(searchType == 'b'){bestSolution = blockSearch(graph, length, sortedVertices, graphArr, 0, 0, 10000);}
	else if(searchType == 'c'){bestSolution = blockSearch(graph, length, sortedVertices, graphArr, 1, 0, 0);}
	else if(searchType == 'd'){bestSolution = blockSearch(graph, length, sortedVertices, graphArr, 1, 0, 10000);}
	else if(searchType == 'e'){bestSolution = leaf_search(graph, length, sortedVertices, graphArr);}
	else if(searchType == 'f'){bestSolution = blockSearch(graph, length, sortedVertices, graphArr, 2, 0, 0);}
	else if(searchType == 'g'){bestSolution = blockSearch(graph, length, sortedVertices, graphArr, 2, 0, 10000);}
	else if(searchType == 'h'){bestSolution = single_leaf_search(graph, length, sortedVertices, graphArr);}
	else if(searchType == 'i'){bestSolution = pair_leaf_search(graph, length, sortedVertices, graphArr);}
	else if(searchType == 'j'){bestSolution = blockSearch(graph, length, sortedVertices, graphArr, 1, 1, 0);}
	else if(searchType == 'k'){bestSolution = blockSearch(graph, length, sortedVertices, graphArr, 1, 1, 10000);}
	else if(searchType == 'l'){bestSolution = single_leaf_pair_search(graph, length, sortedVertices, graphArr);}
	else if(searchType == 'm'){bestSolution = blockSearch(graph, length, sortedVertices, graphArr, 2, 1, 0);}
	else if(searchType == 'n'){bestSolution = blockSearch(graph, length, sortedVertices, graphArr, 2, 1, 10000);}
	else if(searchType == 'o'){bestSolution = combined_leaf_search(graph, length, sortedVertices, graphArr, weight);}
	else if(searchType == 'p'){bestSolution = blockSearch(graph, length, sortedVertices, graphArr, 1, 2, 0);}
	else if(searchType == 'q'){bestSolution = custom_leaf_search(graph, length, sortedVertices, graphArr, weight);}
	else if(searchType == 'r'){bestSolution = blockSearch(graph, length, sortedVertices, graphArr, 2, 2, 0);}
	else if(searchType == 's'){bestSolution = blockSearch(graph, length, sortedVertices, graphArr, 1, 2, 10000);}
	else if(searchType == 't'){bestSolution = blockSearch(graph, length, sortedVertices, graphArr, 2, 2, 10000);}
	else{printf("Do not understand search type!\n");}

	// delete graphArr;
	return bestSolution;
	




	// if(searchType == 'a'){return blockSearch(graph, length, sortedVertices, graphArr, 0, 0, 0);}
	// else if(searchType == 'b'){return blockSearch(graph, length, sortedVertices, graphArr, 0, 0, 10000);}
	// else if(searchType == 'c'){return blockSearch(graph, length, sortedVertices, graphArr, 1, 0, 0);}
	// else if(searchType == 'd'){return blockSearch(graph, length, sortedVertices, graphArr, 1, 0, 10000);}
	// else if(searchType == 'e'){return leaf_search(graph, length, sortedVertices, graphArr);}
	// else if(searchType == 'f'){return blockSearch(graph, length, sortedVertices, graphArr, 2, 0, 0);}
	// else if(searchType == 'g'){return blockSearch(graph, length, sortedVertices, graphArr, 2, 0, 10000);}
	// else if(searchType == 'h'){return single_leaf_search(graph, length, sortedVertices, graphArr);}
	// else if(searchType == 'i'){return pair_leaf_search(graph, length, sortedVertices, graphArr);}
	// else if(searchType == 'j'){return blockSearch(graph, length, sortedVertices, graphArr, 1, 1, 0);}
	// else if(searchType == 'k'){return blockSearch(graph, length, sortedVertices, graphArr, 1, 1, 10000);}
	// else if(searchType == 'l'){return single_leaf_pair_search(graph, length, sortedVertices, graphArr);}
	// else if(searchType == 'm'){return blockSearch(graph, length, sortedVertices, graphArr, 2, 1, 0);}
	// else if(searchType == 'n'){return blockSearch(graph, length, sortedVertices, graphArr, 2, 1, 10000);}
	// else if(searchType == 'o'){return combined_leaf_search(graph, length, sortedVertices, graphArr, weight);}
	// else if(searchType == 'p'){return blockSearch(graph, length, sortedVertices, graphArr, 1, 2, 0);}
	// else if(searchType == 'q'){return custom_leaf_search(graph, length, sortedVertices, graphArr, weight);}
	// else if(searchType == 'r'){return blockSearch(graph, length, sortedVertices, graphArr, 2, 2, 0);}
	// else if(searchType == 's'){return blockSearch(graph, length, sortedVertices, graphArr, 1, 2, 10000);}
	// else if(searchType == 't'){return blockSearch(graph, length, sortedVertices, graphArr, 2, 2, 10000);}
	// printf("Do not understand search type!\n");
	// pair<list<Vertex>, double> bestSolution;
	// bestSolution.second = DBL_MAX;
	// return bestSolution;





	// if(searchType == 'n'){return node_pair_search(graph, length, sortedVertices, graphArr);}
	// else if(searchType == 'h'){return node_pair_heuristic_search(graph, length, sortedVertices, graphArr, weight);}
	// else if(searchType == 'l'){return leaf_search(graph, length, sortedVertices, graphArr);}
	// else if(searchType == 's'){return single_leaf_search(graph, length, sortedVertices, graphArr);}
	// else if(searchType == 'p'){return pair_leaf_search(graph, length, sortedVertices, graphArr);}
	// else if(searchType == 'c'){return combined_leaf_search(graph, length, sortedVertices, graphArr, weight);}
	// else if(searchType == 'x'){return custom_leaf_search(graph, length, sortedVertices, graphArr, weight);}
	// If got here, search type specified is not defined
	
}

bool doIntersect(vector<int> &v1, vector<int> &v2){
	for(int i = 0; i < v1.size(); i++){
		for(int j = 0; j < v2.size(); j++){
			if(v1[i] == v2[j]){return true;}
		}
	}
	return false;
}



pair<list<Vertex>, double> blockSearch(Graph* graph, double length, vector<Vertex*> &sortedVertices, Poly_Arr &graphArr,
	int startType, int endType, double weight){

	pair<list<Vertex>, double> bestSolution;
	bestSolution.second = DBL_MAX;
	double bestHeurScore = DBL_MAX;

	// Deal with endType == 2
	if(endType == 2){
		pair<list<Vertex>, double> bestSolutionSingle = blockSearch(graph, length, sortedVertices, graphArr, startType, 0, weight);
		pair<list<Vertex>, double> bestSolutionDouble = blockSearch(graph, length, sortedVertices, graphArr, startType, 1, weight);
		if(bestSolutionSingle.second * 2 < bestSolutionDouble.second){
			return bestSolutionSingle;
		}else{
			return bestSolutionDouble;
		}
		
	}

	// Look for the number of l-blocks
	int numBlocks = 0;
	for(int i = 0; i < graph->blockIds.size(); i++){
		if(graph->blockIds[i] > numBlocks){numBlocks = graph->blockIds[i];}
	}
	// Find out which blocks have usable nodes in them, want a reference node
	vector<int> blockExampleNodes(numBlocks+1, -1);
	if(graph->vertices.size() != graph->blockIds.size()){printf("Different size, should not happen!\n");}
	for(int i = 0; i < graph->blockIds.size(); i++){
		if(graph->vertices[i]->usable){
			blockExampleNodes[graph->blockIds[i]] = i;
		}
	}

	// Find single leaf
	int leaf = -1;
	for(int i = 0; i < graph->nodeLeaves.size(); i++){
		if(graph->nodeLeaves[i].size() > 0){
			leaf = graph->nodeLeaves[i][0];
			break;
		}
	}
	if(leaf == -1){return bestSolution;}

	// Go through every pair of blocks
	for(int b1 = 0; b1 < blockExampleNodes.size(); b1++){
		// Check this is a usable block
		if(blockExampleNodes[b1] == -1){continue;}
		// If must be a leaf, check it is
		if(startType == 1 && graph->nodeLeaves[blockExampleNodes[b1]].size() == 0){continue;}
		// If must be a particular leaf, check it is
		if(startType == 2 && !intInVect(graph->nodeLeaves[blockExampleNodes[b1]], leaf)){continue;}

		for(int b2 = 0; b2 < blockExampleNodes.size(); b2++){
			if(b1 == b2){continue;}
			if(blockExampleNodes[b2] == -1){continue;}
			// If any type, don't want repetition
			if(startType == 0 && endType == 0 && b2 <= b1){continue;}
			// If from a leaf to something else, anything should be fine
			// as long as not the same block (as two blocks cannot be part
			// of all the leaves).
			if(startType == 1 && endType == 0){}
			// Similarly from leaf to leaf, want this to be a leaf and that is it
			if(startType == 1 && endType == 1 && graph->nodeLeaves[blockExampleNodes[b2]].size() == 0){continue;}
			// From a leaf to the outside, so this one cannot contain the same leaf
			if(startType == 2 && endType == 0 && intInVect(graph->nodeLeaves[blockExampleNodes[b2]], leaf)){continue;}
			// Want both to be a leaf and to not be a subset of the specified leaf
			if(startType == 2 && endType == 1 && (graph->nodeLeaves[blockExampleNodes[b2]].size() == 0 || intInVect(graph->nodeLeaves[blockExampleNodes[b2]], leaf))){continue;}

			// Have two blocks that can look at! Do so.
			// Show what they are
			// printf("Block %d: ", b1);
			// for(int i = 0; i < graph->blockIds.size(); i++){
			// 	if(graph->blockIds[i] == b1){printf("%d ", i);}
			// }
			// printf("\nBlock %d: ", b2);
			// for(int i = 0; i < graph->blockIds.size(); i++){
			// 	if(graph->blockIds[i] == b2){printf("%d ", i);}
			// }
			// printf("\n\n");

			// Do the heuristic if needed
			int count = 0;
			auto it2 = graph->minCuts.begin();
			for(auto it = graph->minCutLabels.begin(); it != graph->minCutLabels.end() && abs(weight) > TOLERANCE; ++it){
				if((*it)[blockExampleNodes[b1]] != (*it)[blockExampleNodes[b2]]){count++;}
				++it2;
			}

			
			double cutRatio = (double)(graph->minCuts.size() - count) / graph->minCuts.size();
			if(weight*cutRatio > bestSolution.second){continue;}
			pair<list<Vertex>, double> bestPairSolution = findShortestEdge(graph, length, 'n', b1, b2, sortedVertices, graphArr);
			double heurScore = (weight * cutRatio + 1) * bestPairSolution.second;
			if(bestHeurScore > heurScore){
				bestSolution = bestPairSolution;
				bestHeurScore = heurScore;
			}
		}
	}
	// Might not have anything as solution is not possible!
	if(bestSolution.second == DBL_MAX && (startType == 1 || startType == 2) && endType == 1){bestSolution = blockSearch(graph, length, sortedVertices, graphArr, startType, 0, weight);}
	// Might still not have a solution, do respective search!
	if(bestSolution.second == DBL_MAX){
		if(startType == 1 && endType == 0){return leaf_search(graph, length, sortedVertices, graphArr);}
		if(startType == 2 && endType == 0){return single_leaf_search(graph, length, sortedVertices, graphArr);}
		if(startType == 1 && endType == 1){return pair_leaf_search(graph, length, sortedVertices, graphArr);}
		if(startType == 2 && endType == 1){return single_leaf_pair_search(graph, length, sortedVertices, graphArr);}
		

	}
	return bestSolution;
}


pair<list<Vertex>, double> node_pair_search(Graph* graph, double length, vector<Vertex*> &sortedVertices, Poly_Arr &graphArr){
	pair<list<Vertex>, double> bestSolution;
	bestSolution.second = DBL_MAX;
	// Look for the number of l-blocks
	int numBlocks = 0;
	for(int i = 0; i < graph->blockIds.size(); i++){
		if(graph->blockIds[i] > numBlocks){numBlocks = graph->blockIds[i];}
	}
	// Find out which blocks have usable nodes in them
	vector<bool> usable(numBlocks+1, false);
	if(graph->vertices.size() != graph->blockIds.size()){printf("Different size, should not happen!\n");}
	for(int i = 0; i < graph->blockIds.size(); i++){
		usable[graph->blockIds[i]] = usable[graph->blockIds[i]] || graph->vertices[i]->usable;
	}
	// Look at every pair of blocks 
	for(int block1 = 0; block1 <= numBlocks; block1++){
		if(!usable[block1]){continue;}
		for(int block2 = block1 + 1; block1 <= numBlocks; block2++){
			if(!usable[block1]){continue;}
			// Find best edge between these two blocks
			pair<list<Vertex>, double> bestPairSolution = findShortestEdge(graph, length, 'n', block1, block2, sortedVertices, graphArr);
			if(bestPairSolution.second < bestSolution.second){
				bestSolution = bestPairSolution;
			}
		}
	}
	return bestSolution;
}

pair<list<Vertex>, double> node_pair_heuristic_search(Graph* graph, double length, vector<Vertex*> &sortedVertices, Poly_Arr &graphArr, double weight){
	pair<list<Vertex>, double> bestSolution;
	bestSolution.second = DBL_MAX;
	double bestHeurScore = DBL_MAX;
	int maxCuts = 0;
	for(int v1 = 0; v1 < graph->vertices.size(); v1++){
		if(!graph->vertices[v1]->usable){continue;}
		for(int v2 = v1 + 1; v2 < graph->vertices.size(); v2++){
			if(!graph->vertices[v2]->usable){continue;}
			if(graph->blockIds[v1] == graph->blockIds[v2]){continue;}
			int count = 0;
			for(auto it = graph->minCutLabels.begin(); it != graph->minCutLabels.end(); ++it){
				if((*it)[v1] != (*it)[v2]){count++;}
			}
			// Find the shortest edge joining the two and its associated score
			pair<list<Vertex>, double> bestPairSolution = findShortestEdge(graph, length, 'n', v1, v2, sortedVertices, graphArr);
			double cutRatio = (graph->cuts.size() - count) / graph->cuts.size();
			double heurScore = (weight * cutRatio + 1) * bestPairSolution.second;
			if(bestHeurScore > heurScore){
				bestSolution = bestPairSolution;
				bestHeurScore =  heurScore;
			}
			
		}
	}
	return bestSolution;
}

pair<list<Vertex>, double> leaf_search(Graph* graph, double length, vector<Vertex*> &sortedVertices, Poly_Arr &graphArr){	
	pair<list<Vertex>, double> bestSolution;
	bestSolution.second = DBL_MAX;
	// Simply cycle through each leaf and find the best edge out of it
	for(int leaf = 0; leaf < graph->leaves.size(); leaf++){
		printf("Singe leaf %d: ", leaf);
		for(int i = 0; i < graph->leaves[leaf].size(); i++){
			printf("%d ", graph->leaves[leaf][i]);
		}
		printf("\n");
		pair<list<Vertex>, double> bestLeafSolution = findShortestEdge(graph, length, 'l', leaf, -1, sortedVertices, graphArr);
		if(bestLeafSolution.second < bestSolution.second){
			bestSolution = bestLeafSolution;
		}
	}
	return bestSolution;
}

pair<list<Vertex>, double> single_leaf_search(Graph* graph, double length, vector<Vertex*> &sortedVertices, Poly_Arr &graphArr){
	pair<list<Vertex>, double> bestSolution;
	bestSolution.second = DBL_MAX;
	// Point is to find the smallest indexed node which is in a leaf; then use this leaf only
	int leaf = -1;
	for(int i = 0; i < graph->nodeLeaves.size(); i++){
		if(graph->nodeLeaves[i].size() > 0){
			leaf = graph->nodeLeaves[i][0];
			break;
		}
	}
	if(leaf == -1){return bestSolution;}
	return findShortestEdge(graph, length, 'l', leaf, -1, sortedVertices, graphArr);
}

pair<list<Vertex>, double> random_leaf_search(Graph* graph, double length, vector<Vertex*> &sortedVertices, Poly_Arr &graphArr, int seed){	
	pair<list<Vertex>, double> bestSolution;
	bestSolution.second = DBL_MAX;
	double sumSolutionsCost = 0;
	vector<pair<list<Vertex>, double> > foundSolutions;
	foundSolutions.reserve(graph->leaves.size());
	// Simply cycle through each leaf and find the best edge out of it
	for(int leaf = 0; leaf < graph->leaves.size(); leaf++){
		pair<list<Vertex>, double> bestLeafSolution = findShortestEdge(graph, length, 'l', leaf, -1, sortedVertices, graphArr);
		// If if it is not an impossible solution (i.e. if cost is not DBL_MAX)
		if(abs(bestLeafSolution.second - DBL_MAX) < TOLERANCE){continue;}
		// Store the solution and add the cost
		sumSolutionsCost += bestLeafSolution.second;
		foundSolutions.push_back(bestLeafSolution);
	}
	// Pick a random number and choose edge accordingly
	mt19937 mt_rand(seed);
	double chosenValue = sumSolutionsCost * ((double) mt_rand() / mt_rand.max());
	// Cycle through solutions and find the chosen one
	for(int i = 0; i < foundSolutions.size(); i++){
		chosenValue -= foundSolutions[i].second;
		if(chosenValue <= 0){return foundSolutions[i];}
	}
	return bestSolution;
}

pair<list<Vertex>, double> pair_leaf_search(Graph* graph, double length, vector<Vertex*> &sortedVertices, Poly_Arr &graphArr){
	pair<list<Vertex>, double> bestSolution;
	bestSolution.second = DBL_MAX;
	// Look at every pair of leaves
	for(int leaf1 = 0; leaf1 < graph->leaves.size(); leaf1++){
		for(int leaf2 = leaf1 + 1; leaf2 < graph->leaves.size(); leaf2++){
			printf("Pair leaves %d, %d\n", leaf1, leaf2);
			// printf("%d %d\n", graph->leaves[leaf1][0], graph->leaves[leaf2][0]);
			// if(graph->leaves[leaf1][0] != 0 || graph->leaves[leaf2][0] != 3){continue;}
			pair<list<Vertex>, double> bestLeafSolution = findShortestEdge(graph, length, 'p', leaf1, leaf2, sortedVertices, graphArr);
			if(bestLeafSolution.second < bestSolution.second){
				bestSolution = bestLeafSolution;
			}
		}
	}
	// If still no edge, might be that not two pairs of edges exist, so look for a single one
	if(bestSolution.second == DBL_MAX){return leaf_search(graph, length, sortedVertices, graphArr);}
	return bestSolution;
}

pair<list<Vertex>, double> combined_leaf_search(Graph* graph, double length, vector<Vertex*> &sortedVertices, Poly_Arr &graphArr, double weight){
	pair<list<Vertex>, double> bestSolutionSingle = leaf_search(graph, length, sortedVertices, graphArr);
	pair<list<Vertex>, double> bestSolutionDouble = pair_leaf_search(graph, length, sortedVertices, graphArr);
	if(bestSolutionDouble.second <= weight*bestSolutionSingle.second){return bestSolutionDouble;}
	return bestSolutionSingle;
}

pair<list<Vertex>, double> random_combined_leaf_search(Graph* graph, double length, vector<Vertex*> &sortedVertices, Poly_Arr &graphArr, double weight, int seed){
	pair<list<Vertex>, double> bestSolutionSingle = random_leaf_search(graph, length, sortedVertices, graphArr, seed);
	pair<list<Vertex>, double> bestSolutionDouble = random_leaf_pair_search(graph, length, sortedVertices, graphArr, seed);
	if(bestSolutionDouble.second <= weight*bestSolutionSingle.second){return bestSolutionDouble;}
	return bestSolutionSingle;	
}

pair<list<Vertex>, double> custom_leaf_search(Graph* graph, double length, vector<Vertex*> &sortedVertices, Poly_Arr &graphArr, double weight){
	pair<list<Vertex>, double> bestSolutionSingle, bestSolutionDouble;
	bestSolutionDouble.second = DBL_MAX;
	bestSolutionSingle.second = DBL_MAX;
	// First find single edge out of smallest index leaf
	int leaf = -1;
	for(int i = 0; i < graph->nodeLeaves.size(); i++){
		if(graph->nodeLeaves[i].size() > 0){
			leaf = graph->nodeLeaves[i][0];
			break;
		}
	}
	if(leaf == -1){return bestSolutionSingle;}
	bestSolutionSingle =  findShortestEdge(graph, length, 'l', leaf, -1, sortedVertices, graphArr);
	// Now look at other leaves
	for(int leaf2 = 0; leaf2 < graph->leaves.size(); leaf2++){
		if(leaf2 == leaf){continue;}
		pair<list<Vertex>, double> bestLeafSolution = findShortestEdge(graph, length, 'p', leaf, leaf2, sortedVertices, graphArr);
		if(bestLeafSolution.second < bestSolutionDouble.second){
			bestSolutionDouble = bestLeafSolution;
		}
	}
	if(bestSolutionDouble.second <= weight*bestSolutionSingle.second){return bestSolutionDouble;}
	return bestSolutionSingle;
}

pair<list<Vertex>, double> single_leaf_pair_search(Graph* graph, double length, vector<Vertex*> &sortedVertices, Poly_Arr &graphArr){
	pair<list<Vertex>, double> bestSolutionSingle, bestSolutionDouble;
	bestSolutionDouble.second = DBL_MAX;
	bestSolutionSingle.second = DBL_MAX;
	// First find single edge out of smallest index leaf
	int leaf = -1;
	for(int i = 0; i < graph->nodeLeaves.size(); i++){
		if(graph->nodeLeaves[i].size() > 0){
			leaf = graph->nodeLeaves[i][0];
			break;
		}
	}
	if(leaf == -1){return bestSolutionSingle;}
	// Now look at other leaves
	for(int leaf2 = 0; leaf2 < graph->leaves.size(); leaf2++){
		if(leaf2 == leaf){continue;}
		pair<list<Vertex>, double> bestLeafSolution = findShortestEdge(graph, length, 'p', leaf, leaf2, sortedVertices, graphArr);
		if(bestLeafSolution.second < bestSolutionDouble.second){
			bestSolutionDouble = bestLeafSolution;
		}
	}
	if(bestSolutionDouble.second == DBL_MAX){return single_leaf_search(graph, length, sortedVertices, graphArr);}
	return bestSolutionDouble;
}


pair<list<Vertex>, double> random_leaf_pair_search(Graph* graph, double length, vector<Vertex*> &sortedVertices, Poly_Arr &graphArr, int seed){
	pair<list<Vertex>, double> bestSolution;
	bestSolution.second = DBL_MAX;
	double sumSolutionsCost = 0;
	vector<pair<list<Vertex>, double> > foundSolutions;
	foundSolutions.reserve(graph->leaves.size() * (graph->leaves.size() - 1));
	// Look at every pair of leaves
	for(int leaf1 = 0; leaf1 < graph->leaves.size(); leaf1++){
		for(int leaf2 = leaf1 + 1; leaf2 < graph->leaves.size(); leaf2++){
			// printf("%d %d\n", graph->leaves[leaf1][0], graph->leaves[leaf2][0]);
			// if(graph->leaves[leaf1][0] != 0 || graph->leaves[leaf2][0] != 3){continue;}
			pair<list<Vertex>, double> bestLeafSolution = findShortestEdge(graph, length, 'p', leaf1, leaf2, sortedVertices, graphArr);
			// If if it is not an impossible solution (i.e. if cost is not DBL_MAX)
			if(abs(bestLeafSolution.second - DBL_MAX) < TOLERANCE){continue;}
			// Store the solution and add the cost
			sumSolutionsCost += bestLeafSolution.second;
			foundSolutions.push_back(bestLeafSolution);
		}
	}
	// If still no edge, might be that not two pairs of edges exist, so look for a single one
	if(foundSolutions.size() == 0){return random_leaf_search(graph, length, sortedVertices, graphArr, seed);}
	// Pick a random number and choose edge accordingly
	mt19937 mt_rand(seed);
	double chosenValue = sumSolutionsCost * ((double) mt_rand() / mt_rand.max());
	// Cycle through solutions and find the chosen one
	for(int i = 0; i < foundSolutions.size(); i++){
		chosenValue -= foundSolutions[i].second;
		if(chosenValue <= 0){return foundSolutions[i];}
	}
	return bestSolution;
}


tuple<double, double, double> randomAugment(Graph* graph, double length, double weight, bool printInfo){
	// Want to do the same as normal augment, but want to repeat it a few times
	makeAdjacencyMatrix(graph);
	if(!isConnected(graph)){
		printf("Graph is not connected! Ending now.\n");
		return make_tuple(NAN, 0, 0);
	}
	double bestCost = DBL_MAX;
	Graph* bestGraphAddress;
	int numChecks = 5;
	clock_t t;
	double timeStruct = 0;
	double timeAugment = 0;

	for(int i = 1; i <= numChecks; i++){
		double cost = 0;
		// Create a new copy of the graph to work with
		Graph* newGraph = new Graph();
		newGraph->vertices.reserve(graph->vertices.size());
		newGraph->edges.reserve(graph->edges.size());
		for(int i = 0; i < graph->vertices.size(); i++){
			Vertex* newVertex = new Vertex(graph->vertices[i]->id, graph->vertices[i]->x, graph->vertices[i]->y, graph->vertices[i]->usable);
			newGraph->vertices.push_back(newVertex);
		}
		for(int i = 0; i < graph->edges.size(); i++){
			Line* newEdge = new Line(graph->edges[i]->id, newGraph->vertices[graph->edges[i]->start->id], newGraph->vertices[graph->edges[i]->end->id]);
			newGraph->edges.push_back(newEdge);
		}
		// Augment as usual
		int count = 0;
		pair<list<Vertex>, double> prevEdge;
		do{
			if(printInfo){printf("\n\n================== COUNT IS %d ================== \n\n", count);}
			count++;
			t = clock();
			resetGraph(newGraph);
			makeAdjacencyMatrix(newGraph);
			l_cut_finder(newGraph, length);
			// Skip the next one 
			// findLabelsAndMinCuts(newGraph);
			findLeaves(newGraph);
			timeStruct += ((float)(clock() - t))/CLOCKS_PER_SEC;
			if(newGraph->minCuts.begin() == newGraph->minCuts.end()){break;}

			t = clock();
			// Adding stuff from find_best_edge function
			// Will need the vector of vertices sorted for reference
			vector<Vertex*> sortedVertices = newGraph->vertices;
			sort(sortedVertices.begin(), sortedVertices.end(), compareVertices);

			// Make an arrangement out of the graph for visibility purposes
			Poly_Arr graphArr;
			// Add the edges to the arrangement
			for(int i = 0; i < newGraph->edges.size(); i++){
				Rat_Segment newEdge = Rat_Segment(Rat_Point(newGraph->edges[i]->start->x, newGraph->edges[i]->start->y),
											Rat_Point(newGraph->edges[i]->end->x, newGraph->edges[i]->end->y));
				insert(graphArr, newEdge);
			}
			// Will need to add "boundary" box for outside face, will need min and max x,y coordinates
			double x_min = DBL_MAX, y_min = DBL_MIN, x_max = -DBL_MIN, y_max = -DBL_MIN;
			for(int i = 0; i < newGraph->vertices.size(); i++){
				if(newGraph->vertices[i]->x > x_max){x_max = newGraph->vertices[i]->x;}
				if(newGraph->vertices[i]->x < x_min){x_min = newGraph->vertices[i]->x;}
				if(newGraph->vertices[i]->y > y_max){y_max = newGraph->vertices[i]->y;}
				if(newGraph->vertices[i]->y < y_min){x_max = newGraph->vertices[i]->y;}
			}
			// Add four "box" edges
			insert(graphArr, Rat_Segment(Rat_Point(x_min - 2 * length, y_min - 2 * length),
														Rat_Point(x_max + 2 * length, y_min - 2 * length)));
			insert(graphArr, Rat_Segment(Rat_Point(x_max + 2 * length, y_min - 2 * length),
														Rat_Point(x_max + 2 * length, y_max + 2 * length)));
			insert(graphArr, Rat_Segment(Rat_Point(x_max + 2 * length, y_max + 2 * length),
														Rat_Point(x_min - 2 * length, y_max + 2 * length)));
			insert(graphArr, Rat_Segment(Rat_Point(x_min - 2 * length, y_max + 2 * length),
														Rat_Point(x_min - 2 * length, y_min - 2 * length)));
			// Go through each face and label it as worth analysing (i.e. might be able to add an edge here)
			// Exclude unbounded (i.e. outside of box).
			for (auto fit = graphArr.faces_begin(); fit != graphArr.faces_end(); ++fit){fit->set_data(!(fit->is_unbounded()));}
			pair<list<Vertex>, double> edge = random_combined_leaf_search(newGraph, length, sortedVertices, graphArr, weight, i);
			
			if(count > 0){
				bool repeated = true;
				auto it = edge.first.begin();
				auto it2 = prevEdge.first.begin();
				while(it != edge.first.end() && it2 != prevEdge.first.end()){
					if(it->id != it2->id){repeated = false; break;}
					++it;
					++it2;
				}
				repeated = repeated && it == edge.first.end() && it2 == prevEdge.first.end();
				if(repeated){
					edge.second = DBL_MAX;
					printf("Stuck in an infite loop!\n");
				}
				
			}
			// prevEdge = edge;

			if(edge.second == DBL_MAX){
				printf("Graph cannot be made l-resilient with current search type!\n");
				break;
			}
			// Add new edge to graph	
			addSegment(newGraph, edge.first);
			timeAugment += ((float)(clock() - t))/CLOCKS_PER_SEC;
			// printf("Best candidate edge has cost %.2f and nodes ", edge.second);
			// for(auto it = edge.first.begin(); it != edge.first.end(); ++it){
			// 	printf("%d ", it->id);
			// }
			// printf("\n");
			// I feel like storing previous edge needs to happen here instead
			prevEdge = edge;	
			cost += edge.second;
		
		}while(true);
		// printf("Cost %.2f\n", cost);
		// Check if found a better solution!
		if(abs(bestCost - DBL_MAX) < TOLERANCE){
			// printf("Improved first!\n");
			bestGraphAddress = newGraph;
			bestCost = cost;
		}else if(cost < bestCost){
			// printf("Improved!\n");
			delete bestGraphAddress;
			bestGraphAddress = newGraph;
			bestCost = cost;
		}else{
			delete newGraph;
		}
	}
	// Clear graph and assign all values, so that it can be plotted outside this function
	for(int i = 0; i < graph->vertices.size(); i++){
		delete graph->vertices[i];
		graph->vertices[i] = NULL;
	}
	for(int i = 0; i < graph->edges.size(); i++){
		delete graph->edges[i];
		graph->edges[i] = NULL;
	}
	graph->vertices.erase(graph->vertices.begin(), graph->vertices.end());
	graph->edges.erase(graph->edges.begin(), graph->edges.end());

	graph->vertices.reserve(bestGraphAddress->vertices.size());
	graph->edges.reserve(bestGraphAddress->edges.size());

	for(int i = 0; i < bestGraphAddress->vertices.size(); i++){
		Vertex* newVertex = new Vertex(bestGraphAddress->vertices[i]->id, bestGraphAddress->vertices[i]->x, bestGraphAddress->vertices[i]->y, bestGraphAddress->vertices[i]->usable);
		graph->vertices.push_back(newVertex);
	}
	for(int i = 0; i < bestGraphAddress->edges.size(); i++){
		Line* newEdge = new Line(bestGraphAddress->edges[i]->id, graph->vertices[bestGraphAddress->edges[i]->start->id], graph->vertices[bestGraphAddress->edges[i]->end->id]);
		graph->edges.push_back(newEdge);
	}
	
	delete bestGraphAddress;
	resetGraph(graph);
	// printf("One\n");
	makeAdjacencyMatrix(graph);
	// printf("Two\n");
	l_cut_finder(graph, length);
	// printf("Three\n");
	findLabelsAndMinCuts(graph);
	// printf("Four\n");
	findLeaves(graph);
	
	// return bestCost;
	return make_tuple(bestCost, timeStruct, timeAugment);
}




tuple<double, double, double> augment(Graph* graph, double length, char searchType, double weight, bool printInfo){
	if(searchType == 'u'){return randomAugment(graph, length, weight, printInfo);}
	makeAdjacencyMatrix(graph);
	if(!isConnected(graph)){
		printf("Graph is not connected! Ending now.\n");
		return make_tuple(NAN, 0, 0);
	}
	int count = 0;
	double cost = 0;
	pair<list<Vertex>, double> prevEdge;
	pair<list<Vertex>, double> edge;
	clock_t t;
	double timeStruct = 0;
	double timeAugment = 0;

	do{
		if(printInfo){printf("\n\n================== COUNT IS %d ================== \n\n", count);}
		
		// Reset values
		resetGraph(graph);
		// Make adjacency matrix
		makeAdjacencyMatrix(graph);

		// Slight variation based on the search,
		// to incorporate fully polynomial scheme
		if(searchType == 'v' || searchType == 'w'){
			t = clock();
			// Run l_cut_finder but only to find leaves
			l_cut_finder(graph, length, true);
			// Now need to save data in nodeLeaves
			storeNodeLeaves(graph);
			timeStruct += ((float)(clock() - t))/CLOCKS_PER_SEC;
			t = clock();
			// Finally if the graph has no leaves, it is l-resilient and can stop!
			if(graph->leaves.begin() == graph->leaves.end()){break;}
			// Find best edge based on search type
			if(searchType == 'v'){edge = find_best_edge(graph, length, 'o', weight);}
			else if(searchType == 'w'){edge = find_best_edge(graph, length, 'q', weight);}
			timeAugment += ((float)(clock() - t))/CLOCKS_PER_SEC;
			// printf("Have leaves (when only looking for leaves):\n");
			// for(int i = 0; i < graph->leaves.size(); i++){
			// 	for(int j = 0; j < graph->leaves[i].size(); j++){
			// 		printf("%d (%.2f, %.2f)", graph->leaves[i][j], graph->vertices[graph->leaves[i][j]]->x, graph->vertices[graph->leaves[i][j]]->y);
			// 	}
			// 	printf(" - with cut ");
			// 	for(int j = 0; j < graph->leafMinCut[i].size(); j++){
			// 		printf("%d (%d-%d) ", graph->leafMinCut[i][j], graph->edges[graph->leafMinCut[i][j]]->start->id, graph->edges[graph->leafMinCut[i][j]]->end->id);
			// 	}
			// 	printf("\n");
			// }

			// Reset values
			// resetGraph(graph);
			// Make adjacency matrix
			// makeAdjacencyMatrix(graph);
		// }
			
		}else{
			t = clock();
			// Find all cuts
			// printf("Cuts\n");
			l_cut_finder(graph, length);
			// Find min cuts
			// printf("Labels\n");
			findLabelsAndMinCuts(graph);
			// Find leaves
			// printf("Leaves\n");
			findLeaves(graph);
			timeStruct += ((float)(clock() - t))/CLOCKS_PER_SEC;
			t = clock();
			// if(count ==0){break;}
			// If have no cuts, the graph is l-resilient and can stop!
			if(graph->minCuts.begin() == graph->minCuts.end()){break;}
			// Find best edge based on search type
			// if(searchType == 'v'){edge = find_best_edge(graph, length, 'o', weight);}
			// else if(searchType == 'w'){edge = find_best_edge(graph, length, 'q', weight);}
			// else{edge = find_best_edge(graph, length, searchType, weight);}
			edge = find_best_edge(graph, length, searchType, weight);
			timeAugment += ((float)(clock() - t))/CLOCKS_PER_SEC;
			// printf("Have leaves (when only looking for everything):\n");
			// for(int i = 0; i < graph->leaves.size(); i++){
			// 	for(int j = 0; j < graph->leaves[i].size(); j++){
			// 		printf("%d (%.2f, %.2f)", graph->leaves[i][j], graph->vertices[graph->leaves[i][j]]->x, graph->vertices[graph->leaves[i][j]]->y);
			// 	}
			// 	printf(" - with cut ");
			// 	for(int j = 0; j < graph->leafMinCut[i].size(); j++){
			// 		printf("%d (%d-%d) ", graph->leafMinCut[i][j], graph->edges[graph->leafMinCut[i][j]]->start->id, graph->edges[graph->leafMinCut[i][j]]->end->id);
			// 	}
			// 	printf("\n");
			// }
		}
		
		
		if(count > 0){
			bool repeated = true;
			auto it = edge.first.begin();
			auto it2 = prevEdge.first.begin();
			while(it != edge.first.end() && it2 != prevEdge.first.end()){
				if(it->id != it2->id){repeated = false; break;}
				++it;
				++it2;
			}
			repeated = repeated && it == edge.first.end() && it2 == prevEdge.first.end();
			if(repeated){
				printf("Stuck in an infite loop! New edge ");
				it = edge.first.begin();
				while(it != edge.first.end()){
					printf("%d ", it->id);
					++it;
				}
				printf("\n");
				printf("Previous edge ");
				it = prevEdge.first.begin();
				while(it != prevEdge.first.end()){
					printf("%d ", it->id);
					++it;
				}
				printf("\n");
				return make_tuple(NAN, 0, 0);
			}	
		}
		// prevEdge = edge;
		// If could not find an edge by this stage, graph cannot be made l-resilient
		// with current technique
		if(edge.second == DBL_MAX){
			printf("Graph cannot be made l-resilient with current search type!\n");
			return make_tuple(DBL_MAX, 0, 0);
		}
		// Add new edge to graph	
		addSegment(graph, edge.first);
		if(printInfo){
			printf("Best candidate edge has cost %.2f and nodes ", edge.second);
			for(auto it = edge.first.begin(); it != edge.first.end(); ++it){
				// printf("%d (%.4f, %.4f)", it->id, graph->vertices[it->id]->x, graph->vertices[it->id]->y);
				printf("%d ", it->id);
			}
			printf("\n");
		}
		// I feel like storing previous edge needs to happen here instead
		prevEdge = edge;
		count++;
		cost += edge.second;
		
	}while(true);
	// Here want to do one final search for l-cuts to make sure none are left!
	return make_tuple(cost, timeStruct, timeAugment);
}

void connectGraph(Graph* graph){
	// First of all create the adjacency matrix
	makeAdjacencyMatrix(graph);
	// Check if it is connected, if so can stop
	vector<int> visited(graph->vertices.size(), 0);
	vector< pair <int,int> > forbidden(0);
	int group = 1;
	for(int i = 0; i < graph->vertices.size(); i++){
		if(visited[i] != 0){continue;}
		depthFirstSearch(graph, forbidden, visited, i, group);
		group++;
	}
	group--;
	// Graph is connected when there is only one group
	while(group != 1){
		// printf("Group %d\n", group);
		double bestEdgeCost = DBL_MAX;
		pair<int, int> edge;
		// Find shortest edge which does not intersect with current edges
		for(int i = 0; i < visited.size(); i++){
			for(int j = i + 1; j < visited.size(); j++){
				if(visited[i] == visited[j]){continue;}
				bool intersect = false;
				for(int k = 0; k < graph->edges.size() && !intersect; k++){
					// If share endpoints, will intersect there and nowhere else
					// Can ignore
					if(graph->edges[k]->start->id == graph->vertices[i]->id ||
						graph->edges[k]->start->id == graph->vertices[j]->id ||
						graph->edges[k]->end->id == graph->vertices[i]->id ||
						graph->edges[k]->end->id == graph->vertices[j]->id){continue;}
					// Check if intersects
					intersect = doIntersect(graph->vertices[i], graph->vertices[j],
					 						graph->edges[k]->start, graph->edges[k]->end);
				}
				// If intersected any of the edges, ignore
				if(intersect){continue;}
				// By this stage might have a better edge
				if(bestEdgeCost > edgeLength(graph->vertices[i], graph->vertices[j])){
					bestEdgeCost = edgeLength(graph->vertices[i], graph->vertices[j]);
					edge = make_pair(i,j);
				}
			}
		}
		// Add edge to graph
		Line* l = new Line(graph->edges.size(), graph->vertices[edge.first], graph->vertices[edge.second]);
		graph->edges.push_back(l);
		// Need to rewrite "visit" vector
		int before = visited[edge.first];
		int after = visited[edge.second];
		for(int i = 0; i < visited.size(); i++){
			if(visited[i] == before){visited[i] = after;}
		}
		group--;

	}
	resetGraph(graph);
}


vector<string> processInputString(string input){
	vector<string> splitString;
	// First check the format is right, i.e. starts and ends with {}
	if(input.front() != '{' || input.back() != '}'){return splitString;}
	// Now get subset
	string inside = input.substr(1, (int)input.size() - 2);
	stringstream ss(inside);
	string token;
	while(getline(ss, token, ',')){
		splitString.push_back(token);
	}
	return splitString;
}


int main(int argc, char** argv){


	// Have section which takes input. Will not let this run without any inputs, so start with an error message
	if(argc == 1){
		printf("Error: This program needs to run with at least one input! The options are...\n");
		return 0;
	}
	// At this point have at least one input, read it to know what's going on
	string mode = argv[1];

	// Just realised I might get serious problems if not all graph files already exist, so now need a section which does not augment anything, only creates graph files
	if(mode.compare("graphCreation") == 0){
		if(argc != 3){
			printf("Error: Called program with mode graphCreation (i.e. read from file with each of the required runs, and create graphs), which should look like ");
			printf("./augmentation graphCreation runName\nWanted 3 inputs but got %d! Stopping now...\n", argc);
			printf("Found inputs:\n");
			for(int i = 0; i < argc; i++){
				printf("%s\n", argv[i]);
			}
			return 0;
		}
		string runName = argv[2];
		string inputFilename = "../data/runScripts/" + runName + "-runScript.txt";
		ifstream inputFile(inputFilename, ios::in);
		if(!inputFile.is_open()){printf("Error: Could not open input file %s to read experiment specifications! Stopping now...\n", inputFilename.c_str()); return 0;}
		string line;
		Graph* graph;
		// Skip first line
		getline(inputFile, line);
		while(getline(inputFile, line)){
			// Read info
			stringstream ss(line);
			string token;
			getline(ss, token, ' ');
			int size = stoi(token);
			getline(ss, token, ' ');
			stringstream ss2(token);
			getline(ss2, token, '-');
			int seedStart = stoi(token);
			getline(ss2, token, '-');
			int seedEnd = stoi(token);
			getline(ss, token, ' ');
			char inputType = token[0];
			// Skip next, it is disaster
			getline(ss, token, ' ');
			// Skip next, it is search
			getline(ss, token, ' ');
			getline(ss, token, ' ');
			int xDim = stoi(token);
			getline(ss, token, ' ');
			int yDim = stoi(token);
			for(int seed = seedStart; seed <= seedEnd; seed++){
				graph = new Graph();
				string graphFilename;
				if(inputType == 'M'){graphFilename = "../data/graphs/graph_n" + to_string(size) + "_s" + to_string(seed) + "_MST_" + to_string(xDim) + "x" + to_string(yDim) + ".txt";}
				else{graphFilename = "../data/graphs/graph_n" + to_string(size) + "_s" + to_string(seed) + "_TopBottom_" + to_string(xDim) + "x" + to_string(yDim) + ".txt";}
				if(read_graph(graph, graphFilename) == 1){
					printf("Creating graph filename %s\n", graphFilename.c_str());
					if(inputType == 'M'){generate_graph(graph, size, 0, xDim, yDim, seed, 0);}
					else{generate_graph_top_bottom(graph, size, xDim, yDim, seed);}
					connectGraph(graph);
					write_graph(graph, graphFilename);
				}
				delete graph;
			}
		}

	}

	if(mode.compare("run") == 0){
		if(argc != 6){
			printf("Error: Called program with mode run (i.e. read from file with each of the required runs, and run specified line), which should look like ");
			printf("./augmentation run runName lineNum printInfo plotGraph\nWanted 6 inputs but got %d! Stopping now...\n", argc);
			printf("Found inputs:\n");
			for(int i = 0; i < argc; i++){
				printf("%s\n", argv[i]);
			}
			return 0;
		}
		// Open file where 
		string runName = argv[2];
		int arrayNum = stoi(argv[3]);
		bool printInfo = (bool)stoi(argv[4]);
		bool plotGraph = (bool)stoi(argv[5]);
		// Check arrayNum is not 0, as this will get the header
		if(arrayNum == 0){printf("Error: Specified first line of file as inputs to run, but this is the header line! Stopping now...\n"); return 0;}
		// Now need to get to the file in read mode, and go to the right row
		string inputFilename = "../data/runScripts/" + runName + "-runScript.txt";
		ifstream inputFile(inputFilename, ios::in);
		if(!inputFile.is_open()){printf("Error: Could not open input file %s to read experiment specifications! Stopping now...\n", inputFilename.c_str()); return 0;}
		string line;
		int lineNum = 0;
		bool found = false;
		while(getline(inputFile, line)){
			if(lineNum == arrayNum){found = true; break;}
			lineNum++;
		}
		// Check that for the right line insted of being finished
		if(!found){printf("Error: Did not get to the desired array job, specified job %d but file only has %d specifications! Stopping now...\n", arrayNum, lineNum-1); return 0;}
		// line should contain all the desired information, let's extract it
		stringstream ss(line);
		string token;
		getline(ss, token, ' ');
		int size = stoi(token);
		getline(ss, token, ' ');
		stringstream ss2(token);
		getline(ss2, token, '-');
		int seedStart = stoi(token);
		getline(ss2, token, '-');
		int seedEnd = stoi(token);
		getline(ss, token, ' ');
		char inputType = token[0];
		getline(ss, token, ' ');
		double disaster = stof(token);
		getline(ss, token, ' ');
		char search = token[0];
		double weight = 0;
		if(search == 'o' || search == 'q' || search == 'u' || search == 'v' || search == 'w'){
			string stringWeight = token.substr(1, (int)token.size() - 1);
			weight = stof(stringWeight);
		}
		getline(ss, token, ' ');
		int xDim = stoi(token);
		getline(ss, token, ' ');
		int yDim = stoi(token);
		inputFile.close();

		// Get work going 
		initialize();
		Graph* graph;
		// Initialise output file
		string outputFilename = "../data/runResults/" + runName + "_arrayJob" + to_string(arrayNum) + ".txt";
		ofstream outputFile(outputFilename);
		if(!outputFile.is_open()){printf("Error: Could not open output file %s to output experiment results! Stopping now...\n", outputFilename.c_str()); return 0;}
		outputFile << "size seed input_type disaster search xDim yDim augment_time structure_time edge_time cost\n";
		outputFile.close();

		// Cycle through the seeds
		for(int seed = seedStart; seed <= seedEnd; seed++){
			printf("Got specified size %d, seed %d, inputType %c, disasterLength %.2f, search %c, weight %.1f, xDim %d, yDim %d\n", size, seed, inputType, disaster, search, weight, xDim, yDim);
			graph = new Graph();
			string graphFilename;
			// Get graph
			if(inputType == 'M'){graphFilename = "../data/graphs/graph_n" + to_string(size) + "_s" + to_string(seed) + "_MST_" + to_string(xDim) + "x" + to_string(yDim) + ".txt";}
			else{graphFilename = "../data/graphs/graph_n" + to_string(size) + "_s" + to_string(seed) + "_TopBottom_" + to_string(xDim) + "x" + to_string(yDim) + ".txt";}
			if(read_graph(graph, graphFilename) == 1){
				printf("All good though, creating it now.\n");
				generate_graph(graph, size, 0, xDim, yDim, seed, 0);
				connectGraph(graph);
				write_graph(graph, graphFilename);
			}
			if(graph->vertices.size() == 0){
				printf("Wrong read of graph!\n");
				exit(0);
			}
		
			// Augment
			clock_t t = clock();
			tuple<double, double, double> result = augment(graph, disaster, search, weight, printInfo);
			t = clock() - t;
			// Output to file
			outputFile.open(outputFilename, ios_base::app);
			outputFile << to_string(size) << " " << to_string(seed);
			if(inputType == 'M'){outputFile << " MST ";}
			else{outputFile << " TopBottom ";}
			string weightName;
			if(abs(round(weight) - weight) < TOLERANCE){
				weightName = to_string(weight).substr(0,1);
			}else{
				weightName = to_string(weight).substr(0,3);
			}
			outputFile << to_string(disaster) << " " << search << "(" << weightName << ") " << to_string(xDim) << " " << to_string(yDim) << " ";
			outputFile << setprecision(5) << to_string(((float)t)/CLOCKS_PER_SEC) << " " << to_string(get<1>(result)) << " " << to_string(get<2>(result)) << " " << get<0>(result) << "\n";
			outputFile.close();

			string plotFilename;
			if(inputType == 'M'){plotFilename = "../data/plots/n" + to_string(size) + "_s" + to_string(seed) + "_l" + to_string(disaster) + "_MST_search_" + search + "(" + to_string((int)weight) + ").svg";}
			else{plotFilename = "../data/plots/n" + to_string(size) + "_s" + to_string(seed) + "_l" + to_string(disaster) + "_TopBottom_search_" + search + "(" + to_string((int)weight) + ").svg";}

			// plotFilename = "../data/plots/AA_custom_run_n" + to_string(n) + "_s" + to_string(s) + "_l" + to_string(length) + "_search_" + search + "(" + to_string((int)weight) + ").svg";

			if(plotGraph){write_graph_to_svg(graph, disaster, plotFilename, xDim, yDim, true, false, false, true, true, false, false);}
			delete graph;
		}
	

	}else if(mode.compare("setup") == 0){
		if(argc != 10){
			printf("Error: Called program with mode setup (i.e. create file with each of the required runs), which should look like ");
			printf("./augmentation setup runName {sizeStart,sizeEnd,stepSize} {seedsStart,seedsEnd} {M,T} {l_1, ..., l_k} {s_1,...,s_n} xDim yDim\nWanted 9 inputs but got %d! Stopping now...\n", argc);
			printf("Found inputs:\n");
			for(int i = 0; i < argc; i++){
				printf("%s\n", argv[i]);
			}
			return 0;
		}
		// Extract info
		string runName = argv[2];
		vector<string> sizes = processInputString(argv[3]);
		vector<string> seeds = processInputString(argv[4]);
		vector<string> inputTypes = processInputString(argv[5]);
		vector<string> disasterLengths = processInputString(argv[6]);
		vector<string> searchStrategies = processInputString(argv[7]);
		int xDim = stoi(argv[8]);
		int yDim = stoi(argv[9]);

		// Check all went well
		if((int)sizes.size() != 3){printf("Received sizes information, expected {sizeStart,sizeEnd,stepSize}, but got %s (processed to size %d)! Stopping now...\n", argv[2], (int)sizes.size()); return 0;}
		if((int)seeds.size() != 2){printf("Received seeds information, expected {seedsStart,seedsEnd}, but got %s! Stopping now...\n", argv[3]); return 0;}
		if((int)inputTypes.size() == 0 || (int)inputTypes.size() > 2){printf("Received inputType information, expected {M,T}, {M} or {T}, but got %s! Stopping now...\n", argv[4]); return 0;}
		if((int)disasterLengths.size() == 0){printf("Received disaster length information, expected {l_1, ..., l_k}, but got %s! Stopping now...\n", argv[5]); return 0;}
		if((int)searchStrategies.size() == 0){printf("Received search strategies information, expected {s_1,...,s_n}, but got %s! Stopping now...\n", argv[6]); return 0;}

		// Now check inputs
		if(stoi(sizes[0]) <= 0 || stoi(sizes[1]) < stoi(sizes[0]) || stoi(sizes[2]) <= 0){
			printf("Something weird is going on with sizes, want something like {10,100,5}, but got {%d,%d,%d}! Stopping now...\n", stoi(sizes[0]), stoi(sizes[1]), stoi(sizes[2]));
			return 0;
		}
		if(stoi(seeds[0]) <= 0 || stoi(seeds[1]) < stoi(seeds[0])){
			printf("Something weird is going on with seeds, want something like {1,50}, but got {%d,%d}! Stopping now...\n", stoi(seeds[0]), stoi(seeds[1]));
			return 0;
		}
		for(int i = 0; i < (int)inputTypes.size(); i++){
			if((int)inputTypes[i].size() != 1){
				printf("Something weird is going on with inputTypes, want a single character, but got %s! Stopping now...\n", inputTypes[i].c_str());
				return 0;
			}
			if(inputTypes[i][0] != 'M' && inputTypes[i][0] != 'T'){
				printf("Something weird is going on with inputTypes, want either M or T, but got %s! Stopping now...\n", inputTypes[i].c_str());
				return 0;
			}
		}
		for(int i = 0; i < (int)disasterLengths.size(); i++){
			if(stof(disasterLengths[i]) <= 0){
				printf("Something weird is going on with disaster lengths, want >= 0, but got %.2f! Stopping now...\n", stof(disasterLengths[i]));
				return 0;
			}
		}
		for(int i = 0; i < (int)searchStrategies.size(); i++){
			if((int)searchStrategies[i].size() == 0){
				printf("Something weird is going on with search strategies, got a strategy of length 0! Stopping now...\n");
				return 0;
			
			}else if((int)searchStrategies[i].size() == 1 && searchStrategies[i][0] != 'a' &&
																searchStrategies[i][0] != 'b' &&
																searchStrategies[i][0] != 'c' &&
																searchStrategies[i][0] != 'd' &&
																searchStrategies[i][0] != 'e' &&
																searchStrategies[i][0] != 'f' &&
																searchStrategies[i][0] != 'g' &&
																searchStrategies[i][0] != 'h' &&
																searchStrategies[i][0] != 'i' &&
																searchStrategies[i][0] != 'j' &&
																searchStrategies[i][0] != 'k' &&
																searchStrategies[i][0] != 'l' &&
																searchStrategies[i][0] != 'm' &&
																searchStrategies[i][0] != 'n' &&
																searchStrategies[i][0] != 'p' &&
																searchStrategies[i][0] != 'r' &&
																searchStrategies[i][0] != 's' &&
																searchStrategies[i][0] != 't'){
				printf("Something weird is going on with search strategies, got a strategy %s of length 1, but should be one of a, b, c, d, e, f, g, h, i, j, k, l, m, n, p, r, s, t! Stopping now...\n", searchStrategies[i].c_str());
				return 0;
			
			}else if((int)searchStrategies[i].size() > 1 &&
						searchStrategies[i][0] != 'o' &&
						searchStrategies[i][0] != 'q' &&
						searchStrategies[i][0] != 'u' &&
						searchStrategies[i][0] != 'v' &&
						searchStrategies[i][0] != 'w'){
				printf("Something weird is going on with search strategies, got a strategy %s of length > 1, but should be one of o, q, u followed by int weight! Stopping now...\n", searchStrategies[i].c_str());
				return 0;
			}
		}
		
		if(xDim <= 0){
			printf("Something weird is going on with xDim, want > 0, but got %d! Stopping now...\n", xDim);
			return 0;
		}
		if(yDim <= 0){
			printf("Something weird is going on with yDim, want > 0, but got %d! Stopping now...\n", yDim);
			return 0;
		}

		string setupFile = "../data/runScripts/" + runName + "-runScript.txt";
		ofstream outputFile;
		outputFile.open(setupFile);
		if(!outputFile.is_open()){printf("Error: Could not open output file %s to output experiment results! Stopping now...\n", setupFile.c_str()); return 0;}
		outputFile << "size seeds input_type disaster search xDim yDim\n";
		outputFile.close();
		int jobs = 0;
		for(int i = 0; i < (int)searchStrategies.size(); i++){
			for(int j = 0; j < (int)inputTypes.size(); j++){
				for(int k = 0; k < (int)disasterLengths.size(); k++){
					for(int size = stoi(sizes[0]); size <= stoi(sizes[1]); size += stoi(sizes[2])){
						// for(int seed = stoi(seeds[0]); seed <= stoi(seeds[1]); seed++){
						outputFile.open(setupFile, ios_base::app);
						if(!outputFile.is_open()){printf("Error: Could not open output file %s to output experiment results! Stopping now...\n", setupFile.c_str()); return 0;}
						// Adding something here to split up two searches as they take forever to complete
						if((searchStrategies[i][0] == 'a' ||
							searchStrategies[i][0] == 'u') &&
							size > 50){
							outputFile << to_string(size) << " " << seeds[0] << "-" << (int)((stoi(seeds[1]) - stoi(seeds[0]))/5 + stoi(seeds[0])) << " " << inputTypes[j] << " " << disasterLengths[k] << " " << searchStrategies[i] << " " << to_string(xDim) << " " << to_string(yDim) << "\n";
							outputFile << to_string(size) << " " << (int)((stoi(seeds[1]) - stoi(seeds[0]))/5 + stoi(seeds[0]) + 1) << "-" << (int)(2 * (stoi(seeds[1]) - stoi(seeds[0]))/5 + stoi(seeds[0])) << " " << inputTypes[j] << " " << disasterLengths[k] << " " << searchStrategies[i] << " " << to_string(xDim) << " " << to_string(yDim) << "\n";
							outputFile << to_string(size) << " " << (int)(2*(stoi(seeds[1]) - stoi(seeds[0]))/5 + stoi(seeds[0]) + 1) << "-" << (int)(3 * (stoi(seeds[1]) - stoi(seeds[0]))/5 + stoi(seeds[0])) << " " << inputTypes[j] << " " << disasterLengths[k] << " " << searchStrategies[i] << " " << to_string(xDim) << " " << to_string(yDim) << "\n";
							outputFile << to_string(size) << " " << (int)(3*(stoi(seeds[1]) - stoi(seeds[0]))/5 + stoi(seeds[0]) + 1) << "-" << (int)(4 * (stoi(seeds[1]) - stoi(seeds[0]))/5 + stoi(seeds[0])) << " " << inputTypes[j] << " " << disasterLengths[k] << " " << searchStrategies[i] << " " << to_string(xDim) << " " << to_string(yDim) << "\n";
							outputFile << to_string(size) << " " << (int)(4*(stoi(seeds[1]) - stoi(seeds[0]))/5 + stoi(seeds[0]) + 1) << "-" << (int)(5 * (stoi(seeds[1]) - stoi(seeds[0]))/5 + stoi(seeds[0])) << " " << inputTypes[j] << " " << disasterLengths[k] << " " << searchStrategies[i] << " " << to_string(xDim) << " " << to_string(yDim) << "\n";
							jobs += 5;
						}else{
							outputFile << to_string(size) << " " << seeds[0] << "-" << seeds[1] << " " << inputTypes[j] << " " << disasterLengths[k] << " " << searchStrategies[i] << " " << to_string(xDim) << " " << to_string(yDim) << "\n";
							jobs++;
						}
						outputFile.close();
						
						// }
					}
				}
			}
		}
		printf("Successfully setup runScript for job %s with %d jobs!\n", runName.c_str(), jobs);
	}

	return 0;

	// double xDim, yDim, cutSize, maxEdgeLength;
	// int size, edges, seed;
	// int sizeStart, sizeEnd, sizeInterval;
	// double disasterStart, disasterEnd, disasterInterval;
	// int seedStart, seedEnd;
	// clock_t t;







	// if(argc == 1){

	// 	// FOR TESTING
	// 	// Graph* graph = new Graph();
	// 	// read_graph(graph, "../data/graphs/graph_custom.txt");
	// 	// augment(graph, 10, 'a', 0.0);
	// 	// string svgFile = "../data/plots/custom_graph_search_a.svg";
	// 	// write_graph_to_svg(graph, 10, svgFile, 50, 50, true, false, false, true, true, true, true);

	// 	// delete graph;
	// 	// graph = new Graph();
	// 	// read_graph(graph, "../data/graphs/graph_custom.txt");
	// 	// augment(graph, 10, 'b', 0.0);
	// 	// svgFile = "../data/plots/custom_graph_search_b.svg";
	// 	// write_graph_to_svg(graph, 10, svgFile, 50, 50, true, false, false, true, true, true, true);

	// 	// delete graph;

	// 	// return 1;

	// 	xDim = 500.0;
	// 	yDim = 500.0;
	// 	size = 500;
	// 	edges = 1000;
	// 	cutSize = 1.0;

	// 	maxEdgeLength = std::min(xDim,yDim) * 0.3;
	// 	seed = 0;
	
	// }else if(argc == 6){

	// 	xDim = stoi(argv[1]);
	// 	yDim = stoi(argv[2]);
	// 	size = stoi(argv[3]);
	// 	edges = stoi(argv[4]);
	// 	cutSize = stoi(argv[5]);
	// 	maxEdgeLength = std::min(xDim,yDim) * 0.3;
	// 	seed = 0;

	// }else if(argc == 8){

	// 	xDim = stoi(argv[1]);
	// 	yDim = stoi(argv[2]);
	// 	size = stoi(argv[3]);
	// 	edges = stoi(argv[4]);
	// 	cutSize = stoi(argv[5]);
	// 	maxEdgeLength = stoi(argv[6]);
	// 	seed = stoi(argv[7]);

	// }else if(argc == 11){

	// 	xDim = stoi(argv[1]);
	// 	yDim = stoi(argv[2]);
	// 	sizeStart = stoi(argv[3]);
	// 	sizeEnd = stoi(argv[4]);
	// 	sizeInterval = stoi(argv[5]);
	// 	disasterStart = stod(argv[6]);
	// 	disasterEnd = stod(argv[7]);
	// 	disasterInterval = stod(argv[8]);
	// 	seedStart = stoi(argv[9]);
	// 	seedEnd = stoi(argv[10]);
		
	// }else{
	// 	printf("Not right usage\n");
	// 	return 0;
	// }

	// // Initialise lookup table
	// initialize();



	// // HARDCODING TESTING 
	// // vector<int> sizes{60, 60, 60, 60, 60, 60, 60, 70, 70, 70, 70, 70, 70, 70, 70, 80, 80, 80, 80, 80, 80, 80, 80, 80, 80, 80, 80, 80, 80, 80, 80, 80, 90, 90, 90, 90, 90, 90, 90, 90, 90, 90, 90, 90, 90, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 60, 60, 60, 70, 70, 70, 70, 70, 80, 80, 80, 80, 80, 80, 90, 90, 90, 90, 90, 90, 90, 90, 90, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100};
	// // vector<int> seeds{20, 20, 20, 20, 20, 20, 20, 1, 2, 2, 13, 18, 20, 20, 20, 1, 1, 1, 1, 2, 3, 3, 3, 3, 13, 13, 15, 15, 17, 18, 20, 20, 2, 3, 3, 5, 12, 13, 13, 13, 13, 13, 13, 18, 20, 3, 3, 3, 3, 3, 9, 9, 9, 12, 13, 13, 13, 13, 18, 18, 18, 20, 20, 2, 7, 14, 14, 14, 17, 17, 17, 17, 17, 18, 18, 19, 19, 17, 20, 20, 1, 2, 8, 13, 18, 1, 1, 3, 3, 3, 15, 1, 1, 2, 13, 13, 13, 13, 15, 15, 1, 3, 3, 3, 8, 8, 9, 13, 18, 18, 18, 20, 20};
	// // vector<char> inputTypes{'M', 'M', 'M', 'M', 'M', 'M', 'T', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'T', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'T', 'T', 'M', 'M', 'T', 'T', 'M', 'M', 'M', 'M', 'T', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'T', 'M', 'M', 'M', 'M', 'T', 'T', 'T', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'T', 'T', 'T', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'T', 'M', 'M', 'M', 'M', 'M', 'T', 'T', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'T', 'T', 'T', 'M', 'M'};
	// // vector<double> lengths{0.001, 0.001, 0.001, 0.001, 100, 100, 300, 300, 300, 300, 300, 300, 100, 100, 300, 100, 100, 100, 100, 300, 300, 300, 300, 300, 300, 300, 300, 300, 100, 300, 300, 300, 300, 300, 300, 100, 300, 0.001, 0.001, 0.001, 100, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 0.001, 0.001, 5, 100, 100, 100, 100, 300, 300, 300, 300, 300, 300, 300, 5, 5, 5, 100, 100, 300, 300, 300, 300, 100, 0.001, 100, 300, 300, 300, 300, 300, 100, 300, 300, 300, 300, 300, 300, 300, 300, 0.001, 0.001, 100, 300, 300, 300, 300, 300, 300, 300, 100, 100, 300, 0.001, 100, 100, 100, 300, 300};
	// // vector<char> searches{'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'q', 'o', 'o', 'o', 'o', 'o', 'o', 'q', 'o', 'o', 'o', 'o', 'q', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'q', 'o', 'o', 'o', 'o', 'q', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'q', 'o', 'o', 'q', 'o', 'q', 'q', 'q', 'q', 'q', 'q', 'q', 'q', 'q', 'q', 'q', 'q', 'q', 'q', 'q', 'q', 'q', 'q', 'q', 'q', 'q', 'q', 'q', 'q', 'q', 'q', 'q', 'q', 'q', 'q', 'q', 'q', 'q', 'q', 'q', 'q', 'q', 'q', 'q', 'q', 'q', 'q', 'q', 'q', 'q', 'q', 'q', 'q', 'q', 'q'};
	// // vector<char> weights{1, 2, 3, 4, 1, 2, 1, 1, 2, 1, 1, 3, 1, 2, 1, 2, 1, 2, 3, 1, 2, 2, 3, 4, 1, 3, 1, 2, 1, 3, 1, 2, 1, 2, 3, 4, 1, 2, 3, 4, 2, 1, 3, 1, 1, 2, 1, 2, 3, 4, 2, 2, 3, 1, 1, 2, 4, 2, 2, 1, 2, 2, 1, 1, 1, 1, 3, 4, 1, 3, 4, 1, 3, 3, 4, 3, 4, 3, 1, 1, 1, 1, 1, 1, 4, 1, 1, 1, 3, 4, 1, 1, 1, 1, 3, 4, 1, 1, 3, 4, 1, 1, 3, 4, 3, 4, 1, 1, 1, 3, 4, 3, 4};
	// // vector<int> sizes{100, 200, 300, 400, 500};
	// // vector<int> seeds{1, 1, 1, 1, 1};
	// // vector<char> inputTypes{'M', 'M', 'M', 'M', 'M'};
	// // vector<double> lengths{100, 100, 100, 100, 100};
	// // vector<char> searches{'q', 'q', 'q', 'q', 'q'};
	// // vector<char> weights{2, 2, 2, 2, 2};
	// vector<int> sizes{200, 200};
	// vector<int> seeds{1, 1};
	// vector<char> inputTypes{'M', 'M'};
	// vector<double> lengths{100, 100};
	// vector<char> searches{'q', 'o'};
	// vector<char> weights{2, 3};
	
	// string exec_times_filename = "../data/exec_times/custom_run.txt";
	// ofstream file;
	// file.open(exec_times_filename);
	// if(!file.is_open()){
	// 	printf("Couldn't open exec_times_filename write to file\n");
	// 	return 0;
	// }
	// file << "size seed input_type disaster search search_weight augment_time cost\n";
	// file.close();	
	// Graph* graph = new Graph();
	// // Cycle through all the inputs
	// for(int i = 0; i < sizes.size(); i++){
	// 	int n = sizes[i];
	// 	int s = seeds[i];
	// 	double length = lengths[i];
	// 	double cost = 0.0;
	// 	char search = searches[i];
	// 	double weight = weights[i];
	// 	char inputType = inputTypes[i];
	// 	string filename;

	// 	printf("\nn = %d s = %d l = %.3f search = %c (%.2f) type = %c\n\n", n, s, length, search, weight, inputType);

	// 	file.open(exec_times_filename, ios_base::app);
	// 	file << n << " " << s;
	// 	if(inputType == 'M'){file << " MST ";}
	// 	else{file << " TopBottom ";}
	// 	file << length << " " << search << " " << weight << " ";

	// 	delete graph;
	// 	graph = new Graph();

	// 	if(inputType == 'M'){filename = "../data/graphs/graph_n" + to_string(n) + "_s" + to_string(s) + "_MST.txt";}
	// 	else{filename = "../data/graphs/graph_n" + to_string(n) + "_s" + to_string(s) + "_TopBottom.txt";}

	// 	// Attempt to read graph, if could not create and write
	// 	if(read_graph(graph, filename) == 1){
	// 		printf("All good though, creating it now.\n");
	// 		generate_graph(graph, n, 0, xDim, yDim, s, 0);
	// 		connectGraph(graph);
	// 		write_graph(graph, filename);
	// 	}

	// 	// read_graph(graph, filename);
	// 	if(graph->vertices.size() == 0){
	// 		printf("Wrong read of graph!\n");
	// 		exit(0);
	// 	}


	// 	t = clock();
	// 	cost = augment(graph, length, search, weight);
	// 	t = clock() - t;
	// 	file << setprecision(5) << ((float)t)/CLOCKS_PER_SEC << " " << cost << "\n";

	// 	if(inputType == 'M'){filename = "../data/plots/n" + to_string(n) + "_s" + to_string(s) + "_l" + to_string(length) + "_MST_search_" + search + "(" + to_string((int)weight) + ").svg";}
	// 	else{filename = "../data/plots/n" + to_string(n) + "_s" + to_string(s) + "_l" + to_string(length) + "_TopBottom_search_" + search + "(" + to_string((int)weight) + ").svg";}

	// 	filename = "../data/plots/AA_custom_run_n" + to_string(n) + "_s" + to_string(s) + "_l" + to_string(length) + "_search_" + search + "(" + to_string((int)weight) + ").svg";

	// 	write_graph_to_svg(graph, length, filename, xDim, yDim, true, false, false, true, true, false, false);
	// 	file.close();
	// }

	// return 0;

	// UNCOMMENT FROM HERE
	// // vector<char> searches{'n', 'h', 'l', 's', 'p', 'c', 'x'};
	// // vector<double> weights{-1, 1000, -1, -1, -1, 2, 2};
	// // vector<char> searches{'l', 's', 'p', 'c', 'x'};
	// // vector<double> weights{-1, -1, -1, 2, 2};
	// // vector<char> searches{'n', 'h'};
	// // vector<double> weights{-1, 1000};
	// // vector<char> searches{'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o', 'p', 'q', 'r', 's', 't'};
	// // vector<double> weights{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
	// // vector<char> searches{'k', 'n'};
	// // vector<double> weights{0.0, 0.0};


	// // vector<char> searches{'a', 'b', 'c', 'd', 'e', 'i', 'j', 'k', 'p', 'o', 's'};
	// // vector<double> weights{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

	// // vector<char> searches{'q', 'q', 'q', 'q', 'q'};
	// // vector<double> weights{1.0, 2.0, 3.0, 4.0, 5.0};

	// // vector<char> searches{'u', 'q', 'o', 'o', 'o', 'o'};
	// // vector<double> weights{2.0, 2.0, 1.0, 2.0, 3.0, 4.0};

	// //vector<double> lengths{0.001, 5.0, 100.0, 300.0};
	// vector<double> lengths{300.0};

	// vector<char> searches{'q', 'q', 'q', 'u', 'u', 'u'};
	// vector<double> weights{1.0, 3.0, 4.0, 1.0, 3.0, 4.0};

	

	// // Create file in which info gets store
	// string exec_times_filename = "../data/exec_times/run_n" + to_string(sizeStart) + "<-" + to_string(sizeInterval) + "->" + to_string(sizeEnd) + 
	// 									"_s" + to_string(seedStart) + "<->" + to_string(seedEnd) + "_l"  + to_string(disasterStart) + "<-" + 
	// 									to_string(disasterInterval) + "->" + to_string(disasterEnd) + "_xDim" + to_string(xDim) + "_yDim" +
	// 									to_string(yDim) + ".txt";
	// ofstream file;
	// file.open(exec_times_filename);
	// if(!file.is_open()){
	// 	printf("Couldn't open exec_times_filename write to file\n");
	// 	return 0;
	// }
	// file << "size seed input_type disaster search search_weight augment_time cost\n";
	// file.close();	
	// Graph* graph = new Graph();
	// for(int n = sizeStart; n <= sizeEnd; n += sizeInterval){
	// 	for(int s = seedStart; s <= seedEnd; s++){
	// 		// Create graphs and store if they do not exist yet
	// 		// Start with top bottom
	// 		string filename = "../data/graphs/graph_n" + to_string(n) + "_s" + to_string(s) + "_TopBottom.txt";
	// 		delete graph;
	// 		graph = new Graph();
	// 		// Attempt to read graph, if could not create and write
	// 		if(read_graph(graph, filename) == 1){
	// 			printf("All good though, creating it now.\n");
	// 			generate_graph_top_bottom(graph, n, xDim, yDim, s);
	// 			write_graph(graph, filename);
	// 		}
	// 		string svgStartFilename = "../data/plots/n" + to_string(n) + "_s" + to_string(s) + "_TopBottom_START.svg";
	// 		write_graph_to_svg(graph, 0, svgStartFilename, xDim, yDim, true, false, false, false, false, false, false);
	// 		// Next is Minimum Spanning Tree
	// 		filename = "../data/graphs/graph_n" + to_string(n) + "_s" + to_string(s) + "_MST.txt";
	// 		delete graph;
	// 		graph = new Graph();
	// 		if(read_graph(graph, filename) == 1){
	// 			printf("All good though, creating it now.\n");
	// 			generate_graph(graph, n, 0, xDim, yDim, s, 0);
	// 			connectGraph(graph);
	// 			write_graph(graph, filename);
	// 		}
	// 		svgStartFilename = "../data/plots/n" + to_string(n) + "_s" + to_string(s) + "_MST_START.svg";
	// 		write_graph_to_svg(graph, 0, svgStartFilename, xDim, yDim, true, false, false, false, false, false, false);
	// 		// Sweet, now can augment
	// 		// for(double length = disasterStart; length < disasterEnd || abs(length - disasterEnd) < TOLERANCE; length += disasterInterval){
	// 		for(int l = 0; l < lengths.size(); l++){
	// 			double length = lengths[l];
	// 			for(int i = 0; i < searches.size(); i++){
	// 				double cost = 0.0;
	// 				char search = searches[i];
	// 				double weight = weights[i];

	// 				printf("\nn = %d s = %d l = %.3f search = %c (%.2f) type = MST\n\n", n, s, length, search, weight);

	// 				file.open(exec_times_filename, ios_base::app);
	// 				file << n << " " << s << " MST " << length << " " << search << " " << weight << " ";

	// 				delete graph;
	// 				graph = new Graph();
	// 				filename = "../data/graphs/graph_n" + to_string(n) + "_s" + to_string(s) + "_MST.txt";
	// 				read_graph(graph, filename);
	// 				if(graph->vertices.size() == 0){
	// 					printf("Wrong read of graph!\n");
	// 					exit(0);
	// 				}

	// 				t = clock();
	// 				cost = augment(graph, length, search, weight);
	// 				t = clock() - t;
	// 				file << setprecision(5) << ((float)t)/CLOCKS_PER_SEC << " " << cost << "\n";
	// 				filename = "../data/plots/n" + to_string(n) + "_s" + to_string(s) + "_l" + to_string(length) + "_MST_search_" + search + ".svg";
	// 				write_graph_to_svg(graph, length, filename, xDim, yDim, true, false, false, true, true, true, true);
	// 				file.close();

	// 				printf("\nn = %d s = %d l = %.3f search = %c (%.2f) type = TopBottom\n\n", n, s, length, search, weight);

	// 				file.open(exec_times_filename, ios_base::app);
	// 				file << n << " " << s << " TopBottom " << length << " " << search << " " << weight << " ";
	// 				delete graph;
	// 				graph = new Graph();
	// 				filename = "../data/graphs/graph_n" + to_string(n) + "_s" + to_string(s) + "_TopBottom.txt";
	// 				read_graph(graph, filename);
	// 				if(graph->vertices.size() == 0){
	// 					printf("Wrong read of graph!\n");
	// 					exit(0);
	// 				}

	// 				t = clock();
	// 				cost = augment(graph, length, search, weight);
	// 				t = clock() - t;
	// 				file << setprecision(5) << ((float)t)/CLOCKS_PER_SEC << " " << cost << "\n";
	// 				filename = "../data/plots/n" + to_string(n) + "_s" + to_string(s) + "_l" + to_string(length) + "_TopBottom_search_" + search + ".svg";
	// 				write_graph_to_svg(graph, length, filename, xDim, yDim, true, false, true, true, true, true, true);
	// 				file.close();

	// 			}

	// 		}
	// 	}
	// }

	// delete graph;

	// return 0;

	// TO HERE


	// Uncomment this code for a single run with info
	// graph = new Graph();
	// generate_graph(graph, size, edges, xDim, yDim, seed, maxEdgeLength);
	// connectGraph(graph);
	// string startFilename = "../testRuns/augmentationRun_n" + to_string(size) + "_e" + to_string(edges) + "_l" + to_string((int)cutSize) + 
	// 									"_maxEdgeLenght" + to_string((int)maxEdgeLength) + "_X(" + to_string((int)xDim) + ")_Y(" + to_string((int)yDim) + ")_START.svg";
	// string endFilename = "../testRuns/augmentationRun_n" + to_string(size) + "_e" + to_string(edges) + "_l" + to_string((int)cutSize) + 
	// 									"_maxEdgeLenght" + to_string((int)maxEdgeLength) + "_X(" + to_string((int)xDim) + ")_Y(" + to_string((int)yDim) + ")_END.svg";
	// write_graph_to_svg(graph, cutSize, startFilename, xDim, yDim, true, false, false, true, true, false, true);
	// augment(graph, cutSize, 'p', 10);
	// write_graph_to_svg(graph, cutSize, endFilename, xDim, yDim, true, false, true, true, true, false, true);
	
	// printf("Block ids: ");
	// for(int i = 0; i < graph->blockIds.size(); i++){
	// 	printf("%d ", graph->blockIds[i]);
	// }
	// printf("\n");

	// write_graph_to_svg(graph, cutSize, "test.svg", xDim, yDim, true, false, true, true, true);

	// add_edge(graph, cutSize);
	// Print cuts
	// printf("Have %d cuts\n", (int)graph->cuts.size());
	// for(auto it = graph->cuts.begin(); it != graph->cuts.end(); ++it){
	// 	printf("Cut ");
	// 	for(int j = 0; j < it->size(); j++){
	// 		printf("(%d, %d) ", graph->edges[(*it)[j]]->start->id,
	// 						graph->edges[(*it)[j]]->end->id);
	// 	}
	// 	printf("\n");
	// }
	// // Print min cuts
	// printf("Have %d min cuts\n", (int)graph->minCuts.size());
	// for(auto it = graph->minCuts.begin(); it != graph->minCuts.end(); ++it){
	// 	printf("Min cut ");
	// 	for(int j = 0; j < it->size(); j++){
	// 		printf("(%d, %d) ", graph->edges[(*it)[j]]->start->id,
	// 						graph->edges[(*it)[j]]->end->id);
	// 	}
	// 	printf("\n");
	// }


	// // Print leaves
	// printf("Have %d leaves\n", (int)graph->leaves.size());
	// for(int i = 0; i < graph->leaves.size(); i++){
	// 	printf("Leaf ");
	// 	for(int j = 0; j < graph->leaves[i].size(); j++){
	// 		printf("%d ", graph->leaves[i][j]);
	// 	}
	// 	printf("\n");
	// }

	// Uncomment this code for multiple runs with seeds
	// Have used this to find bugs
	// for(int i = 0; i < 100; i++){
	// 	graph = new Graph();
	// 	printf("%d\n",i);
	// 	generate_graph(graph, size, edges, xDim, yDim, i, maxEdgeLength);
	// 	write_graph_to_svg(graph, "test.svg", xDim, yDim);
	// 	makeAdjacencyMatrix(graph);
	// 	l_cut_finder(graph, cutSize);
	// }


}

#endif