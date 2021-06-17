#ifndef INPUT_OUTPUT_CPP
#define INPUT_OUTPUT_CPP

#include "input_output.h"


bool isConnected(vector<vector<bool> > &edges){
	vector<int> visited(edges.size(), 0);
	int group = 1;
	for(int i = 0; i < edges.size(); i++){
		if(visited[i] != 0){continue;}
		depthFirstSearch(edges, visited, i, group);
		group++;
	}
	return group == 2;
}

void depthFirstSearch(vector<vector<bool> > &edges, vector<int> &visited, int vertexId, int group){
	visited[vertexId] = group;
	for(int i = 0; i < edges[vertexId].size(); i++){
		if(edges[vertexId][i] && visited[i] == 0){
			depthFirstSearch(edges, visited, i, group);
		}
	}
}

void generate_graph_top_bottom(Graph* graph, int n, double X, double Y, int seed){
	mt19937 mt_rand(seed);
	Vertex* v;
	Line* l;
	// Start with the vertices
	bool found = true;
	for(int i = 0; i < n && found; i++){
		found = false;
		int tryNum = 0;
		while(tryNum < 100000 && !found){
			found = true;
			tryNum++;
			v = new Vertex();
			v->id = i;
			v->x = X * ((double) mt_rand() / mt_rand.max());
			v->y = Y * ((double) mt_rand() / mt_rand.max());
			// Check that a vertex with the same coordinates has not already been created
			for(int j = 0; j < graph->vertices.size() && found; j++){
				if(std::abs(graph->vertices[j]->x - v->x) < TOLERANCE &&
					std::abs(graph->vertices[j]->y - v->y) < TOLERANCE){
					found = false;
				}
				// Adding another check to make sure the two points are not too close in the y coordinate,
				// as this creates floating point errors
				if(std::abs(graph->vertices[j]->y - v->y) < 0.1){
					found = false;
				}
			}
			if(found){graph->vertices.push_back(v);}
		}
		if(!found){
			printf("Note: With seed %d, wanted %d vertices but could only allocate %d\n",seed,n,(int)graph->vertices.size());
		}
	}
	// Excellent, now assume start with complete graph, and remove edges until there are no intersecting edges
	vector< vector<bool> > edges;
	edges.reserve(n);
	for(int i = 0; i < n; i++){
		vector<bool> edge(n, true);

		edges.push_back(edge);
	}
	// Remove edge with same endpoints
	for(int i = 0; i < n; i++){
		edges[i][i] = false;
	}
	// Remove any edges which have a point in them
	for(int i = 0; i < n; i++){
		for(int j = i + 1; j < n; j++){
			for(int k = 0; k < n; k++){
				if(k == i || k == j){continue;}
				if(orientation(graph->vertices[i], graph->vertices[j], graph->vertices[k]) == 0 &&
					onSegment(graph->vertices[i], graph->vertices[j], graph->vertices[k])){
					edges[i][j] = false;
					edges[j][i] = false;
				}
			}
		}
	}
	// Sweet now remove edges which intersect with another one
	// Choose where to start randomly! Choose edge v1 - v2, 
	// compare with edge v3 - v4
	bool removed = true;
	while(removed){
		removed = false;
		int v1Start = mt_rand() % n;
		int v2Start = mt_rand() % n;
		int v3Start = mt_rand() % n;
		int v4Start = mt_rand() % n;
		for(int i = 0; i < n && !removed; i++){
			int v1 = (v1Start + i) % n;
			
			for(int j = 0; j < n && !removed; j++){
				int v2 = (v2Start + j) % n;
				if(v1 == v2 || edges[v1][v2] == false){continue;}
				
				for(int k = 0; k < n && !removed; k++){
					int v3 = (v3Start + k) % n;
					if(v3 == v1 || v3 == v2){continue;}

					for(int l = 0; l < n && !removed; l++){
						int v4 = (v4Start + l) % n;
						if(v4 == v1 || v4 == v2 || v4 == v3 || edges[v3][v4] == false){continue;}
						// If intersect, remove and stop
						if(doIntersect(graph->vertices[v1], graph->vertices[v2], 
										graph->vertices[v3], graph->vertices[v4])){
							removed = true;
							edges[v1][v2] = false;
							edges[v2][v1] = false;
						}
					}
				}
			} 
		}
	}
	// Need to check if graph is disconnected, if this is the case need to add some edges
	// to reconnect
	vector<int> visited(edges.size(), 0);
	int group = 1;
	for(int i = 0; i < edges.size(); i++){
		if(visited[i] != 0){continue;}
		depthFirstSearch(edges, visited, i, group);
		group++;
	}
	group--;
	while(group != 1){
		bool added = false;
		int v1Start = mt_rand() % n;
		int v2Start = mt_rand() % n;
		int v1, v2;
		for(int i = 0; i < n && !added; i++){
			v1 = (v1Start + i) % n;
			for(int j = 0; j < n && !added; j++){
				v2 = (v2Start + j) % n;
				if(visited[v1] == visited[v2]){continue;}
				// Have two nodes in different connected components,
				// Check if adding this edge intersects with existing edges
				bool intersects = false;
				for(int k = 0; k < n && !intersects; k++){
					if(v1 == k || v2 == k){continue;}
					for(int l = 0; l < n && !intersects; l++){
						if(v1 == l || v2 == l || k == l || edges[k][l] == false){continue;}

						// If intersect, edge is not viable
						if(doIntersect(graph->vertices[v1], graph->vertices[v2], 
										graph->vertices[k], graph->vertices[l])){
							intersects = true;
						}
					}
				}
				if(!intersects){
					edges[v1][v2] = true;
					edges[v2][v1] = true;
					added = true;
				}
			} 
		}
		// Relable groups
		int oldGroup = visited[v2];
		for(int i = 0; i < visited.size(); i++){
			if(visited[i] == oldGroup){visited[i] = visited[v1];}
		}
		group--;
	}

	// Sould have a connected graph
	if(!isConnected(edges)){printf("Graph disconnected!!\n"); exit(0);}
	// Should not have intersecting edges
	for(int i = 0; i < n; i++){
		for(int j = i + 1; j < n; j++){
			for(int k = 0; k < n; k++){
				for(int l = k + 1; l < n; l++){
					if(i == k || i == l || j == k || j == l || edges[i][j] == false || edges[k][l] == false){continue;}
					if(doIntersect(graph->vertices[i], graph->vertices[j], 
										graph->vertices[k], graph->vertices[l])){
						printf("Edges intersect (%d,%d) - (%d, %d), problem!!\n", i, j, k, l);
						exit(0);
					}
				}
			}
		}
	}

	// Count how many edges have, randomly choose how many want to end up with
	int edgeCount = 0;
	for(int i = 0; i < n; i++){
		for(int j = 0; j < n; j++){
			if(edges[i][j]){edgeCount++;}
		}
	}
	if(edgeCount % 2 != 0){printf("Something wrong with edges!!\n"); exit(0);}
	edgeCount = edgeCount / 2;
	// The min number of edges is n - 1, a min spanning tree
	// Max edges is the current number of edges
	int maxEdges = edgeCount;
	int minEdges = n - 1;
	int chosenEdgeNumber = minEdges + (mt_rand() % (maxEdges - minEdges));
	// Sweet now need to remove edges until get the desired number. 
	while(edgeCount != chosenEdgeNumber){
		// So now must choose an edge to remove, note cannot disconnect the graph
		removed = false;
		int v1Start = mt_rand() % n;
		int v2Start = mt_rand() % n;
		// printf("START %d %d\n", v1Start, v2Start);
		for(int i = 0; i < n && !removed; i++){
			for(int j = 0; j < n && !removed; j++){
				int v1 = (v1Start + i) % n;
				int v2 = (v2Start + j) % n;
				if(edges[v1][v2] == false){continue;}

				// printf("Check %d %d\n", v1, v2);
				// Here found a potential edge to remove, remove it and check 
				// graph is still connected
				removed = true;
				edges[v1][v2] = false;
				edges[v2][v1] = false;
				// Check if graph is still connected
				if(isConnected(edges)){continue;}
				// printf("NO\n");
				// If here removing the edge disconnected the graph, so add it
				// again and continue the search
				removed = false;
				edges[v1][v2] = true;
				edges[v2][v1] = true;
			}
		}
		if(!removed){
			printf("Couldn't find edge to remove!! %d %d\n", edgeCount, chosenEdgeNumber); 
			for(int i = 0; i < n; i++){
				for(int j = i + 1; j < n; j++){
					if(edges[i][j]){printf("%d %d\n", i, j);}
				}
			}
			exit(0);}
		edgeCount--;
	}
	// By this stage should be good! Save edges and be done!
	int edgeId = 0;
	for(int i = 0; i < n; i++){
		for(int j = i + 1; j < n; j++){
			if(edges[i][j] == false){continue;}
			l = new Line(edgeId, graph->vertices[i], graph->vertices[j]);
			graph->edges.push_back(l);
			edgeId++;
		}
	}
	if(edgeId != chosenEdgeNumber){printf("Not as many edges as wanted!\n"); exit(0);}
}


void generate_graph(Graph* graph, int n, int e, double X, double Y, int seed, double maxEdgeLength){
	mt19937 mt_rand(seed);
	Vertex* v;
	Line* l;
	bool found = true;
	for(int i = 0; i < n && found; i++){
		found = false;
		int tryNum = 0;
		while(tryNum < 100000 && !found){
			found = true;
			tryNum++;
			v = new Vertex();
			v->id = i;
			v->x = X * ((double) mt_rand() / mt_rand.max());
			v->y = Y * ((double) mt_rand() / mt_rand.max());
			// v->x =  (double) (mt_rand() % (int)X);
			// v->y =  (double) (mt_rand() % (int)Y);
			// Check that a vertex with the same coordinates has not already been created
			for(int j = 0; j < graph->vertices.size() && found; j++){
				if(std::abs(graph->vertices[j]->x - v->x) < TOLERANCE &&
					std::abs(graph->vertices[j]->y - v->y) < TOLERANCE){
					found = false;
				}
				// Adding another check to make sure the two points are not too close in the y coordinate,
				// as this creates floating point errors
				if(std::abs(graph->vertices[j]->y - v->y) < 0.1){
					found = false;
				}
			}
			if(found){graph->vertices.push_back(v);}
		}
		if(!found){
			printf("Note: With seed %d, wanted %d vertices but could only allocate %d\n",seed,n,(int)graph->vertices.size());
		}
	}
	n = (int)graph->vertices.size();
	found = true;
	for(int i = 0; i < e && found; i++){
		int tryNum = 0;
		found = false;
		while(tryNum < 100000 && !found){
			found = true;
			tryNum++;
			// Pick a random pair of nodes, check if edge is valid
			int p1 = mt_rand() % n;
			int p2 = mt_rand() % n;
			// Check it is no longer than specified
			if(edgeLength(graph->vertices[p1],graph->vertices[p2]) > maxEdgeLength){
				found = false;
				continue;
			}
			// Make sure p1 is smaller than p2
			if(p1 == p2){
				found = false;
				continue;
			}
			if(p1 > p2){
				int temp = p1;
				p1 = p2;
				p2 = temp;
			}
			// Make sure other vertex is not on this edge
			for(int j = 0; j < graph->vertices.size() && found; j++){
				if(p1 == graph->vertices[j]->id || p2 == graph->vertices[j]->id){continue;}
				if(orientation(graph->vertices[p1], graph->vertices[p2], graph->vertices[j]) == 0 &&
					onSegment(graph->vertices[p1], graph->vertices[j], graph->vertices[p2])){
					found = false;
				}
			}

			// Check do not have this edge, or does not intersect other edges
			for(int j = 0; j < graph->edges.size() && found; j++){
				if(graph->edges[j]->start->id == p1 && graph->edges[j]->end->id == p2){
					// Got this one already
					found = false;
				}else if(graph->edges[j]->start->id == p1 || graph->edges[j]->end->id == p2){
					// Will intersect at the endpoints, it is ok
					continue;
				}else if(doIntersect(graph->vertices[p1], graph->vertices[p2], 
								graph->edges[j]->start, 
								graph->edges[j]->end)){
					found = false;
				}
			}
			// If "done", save edge, continue
			if(found){
				l = new Line(i, graph->vertices[p1], graph->vertices[p2]);
				graph->edges.push_back(l);
			}
		}
		if(!found){
			printf("Note: With seed %d, wanted %d edges but could only allocate %d\n",seed,e,(int)graph->edges.size());
		}
	}
}


int read_graph(Graph* graph, string filename){
	FILE * fp;
	char buffer[1001];
	char cstr[filename.length() + 1];
	strcpy(cstr, filename.c_str());
	fp = fopen(cstr, "r");
	if(!fp){
		printf("Couldn't open file!\n");
        return 1;
    }
    if(fgets(buffer, 1000, fp) == NULL){printf("Could not read from buffer when skipping line!\n");}
	double x, y;
    int id, vertices, edges, p1, p2;
    sscanf(buffer, "V %d E %d\n\n",&vertices, &edges);

    if(fgets(buffer, 1000, fp) == NULL){printf("Could not read from buffer when skipping line!\n");}
    // if(fgets(buffer, 1000, fp) == NULL){printf("Could not read from buffer when skipping line!\n");}
    Vertex* v;
    for(int i = 0; i < vertices; i++){
    	if(fgets(buffer, 1000, fp) == NULL){printf("Could not read from buffer when reading vertices!\n");}
    	sscanf(buffer, "%d %lf %lf\n",&id, &x, &y);
    	v = new Vertex();
		v->id = id;
		v->x = x;
		v->y = y;
		graph->vertices.push_back(v);
    }
    Line* l;
    if(fgets(buffer, 1000, fp) == NULL){printf("Could not read from buffer when skipping line!\n");}   
    for(int i = 0; i < edges; i++){
    	if(fgets(buffer, 1000, fp) == NULL){printf("Could not read from buffer when reading edges!\n");}
    	sscanf(buffer, "%d %d\n",&p1, &p2);
    	l = new Line(i, graph->vertices[p1], graph->vertices[p2]);
    	graph->edges.push_back(l);
    }

    fclose(fp);

	return 0;
}


void write_graph(Graph* graph, string filename){
	ofstream file;
	// Start with graph file
	file.open(filename);
	if(!file.is_open()){
		printf("Couldn't open graph write to file\n");
	}
	file << "V " << (int)graph->vertices.size() << " E " << (int)graph->edges.size() << "\n\n";
	for(int i = 0; i < graph->vertices.size(); i++){
		file << graph->vertices[i]->id << " " << graph->vertices[i]->x << " " << graph->vertices[i]->y << "\n";
	}
	file << "\n";
	for(int i = 0; i < graph->edges.size(); i++){
		file << graph->edges[i]->start->id << " " << graph->edges[i]->end->id << "\n";
	}
	file.close();
}


void write_graph_to_svg(Graph* graph, double length, string filename, double dimX, double dimY, 
								bool showLeaves, bool showCuts, bool showMinCuts,
								bool showVulnerableRegions,
								bool showDestructionRegions,
								bool showAuxEdges, bool showDisaster){
	
	// Modiy to fit large edges due to disaster
	dimX += 2*length; 
	// Set the defaults so graph looks nice regardless of input
	double shiftX = length + dimX * 0.05;
	double shiftY = dimY * 0.05;
	double limit = std::max(dimX, dimY);
	double pointSize = limit/1000;
	double edgeSize = pointSize / 2;
	int font = (int)limit/100;
	int fontShift = (int)limit/500;
	double leafSize = pointSize * 4;
	double cutEdgeSize = edgeSize*2;
	double vulnerableSize = edgeSize*2;
	double destructionSize = edgeSize*2;
	double destructionLineSize = edgeSize*3;
	
	// Initialise document
	svg::Dimensions dimensions(dimX * 1.1, dimY * 1.1);
    svg::Document doc(filename, svg::Layout(dimensions, svg::Layout::BottomLeft));

    // Write in the vertices
    for(int i = 0; i < graph->vertices.size(); i++){
    	if(!showAuxEdges && !graph->vertices[i]->usable){continue;}
    	doc << svg::Circle(svg::Point(graph->vertices[i]->x + shiftX, graph->vertices[i]->y + shiftY), 4*pointSize,
    	 svg::Color::Black);
    	// doc << svg::Text(svg::Point(graph->vertices[i]->x + shiftX + fontShift, graph->vertices[i]->y + shiftY + fontShift),
    	// 	 to_string(graph->vertices[i]->id), svg::Color::Black, svg::Font(font, "Verdana"));
    }

    // Write edges
    for(int i = 0; i < graph->edges.size(); i++){
    	// int multiplier = (1+ 35*(int)(graph->edges[i]->id > 198));
    	// doc << (svg::Polyline(svg::Stroke(multiplier*edgeSize/6, svg::Color::Black)) << 
    	doc << (svg::Polyline(svg::Stroke(edgeSize/3, svg::Color::Black)) << 
    				svg::Point(graph->edges[i]->start->x + shiftX, graph->edges[i]->start->y + shiftY) << 
    				svg::Point(graph->edges[i]->end->x + shiftX, graph->edges[i]->end->y + shiftY));
    }

    // Add leaves if want to
    for(int i = 0; i < graph->leaves.size() && showLeaves; i++){
    	for(int j = 0; j < graph->leaves[i].size(); j++){

    		doc << svg::Circle(svg::Point(graph->vertices[graph->leaves[i][j]]->x + shiftX,
    									graph->vertices[graph->leaves[i][j]]->y + shiftY), 
    									leafSize,
    	 								svg::Color::Red);
    	}
    }
    // Add cuts if want to
    for(auto it = graph->cuts.begin(); it != graph->cuts.end() && showCuts; ++it){
    	for(int i = 0; i < it->size(); i++){
    		doc << (svg::Polyline(svg::Stroke(cutEdgeSize, svg::Color::Brown)) << 
    				svg::Point(graph->edges[(*it)[i]]->start->x + shiftX, graph->edges[(*it)[i]]->start->y + shiftY) << 
    				svg::Point(graph->edges[(*it)[i]]->end->x + shiftX, graph->edges[(*it)[i]]->end->y + shiftY));
    	}
    	// Add vulnerable regions if want to
    	Rat_Polygon vulnerable = cut_to_vulnerable(graph, (*it), length);
    	for (auto ei = vulnerable.edges_begin(); ei != vulnerable.edges_end() && showVulnerableRegions; ++ei){
    		doc << (svg::Polyline(svg::Stroke(vulnerableSize, svg::Color::Red)) << 
    				svg::Point(CGAL::to_double((*ei).source().x()) + shiftX, CGAL::to_double((*ei).source().y()) + shiftY) << 
    				svg::Point(CGAL::to_double((*ei).target().x()) + shiftX, CGAL::to_double((*ei).target().y()) + shiftY));
    	}
    	// Add destruction regions if want to
    	Rat_Polygon destruction = cut_to_destruction(graph, (*it), length);
    	for (auto ei = destruction.edges_begin(); ei != destruction.edges_end() && showDestructionRegions; ++ei){
    		doc << (svg::Polyline(svg::Stroke(destructionSize, svg::Color::Orange)) << 
    				svg::Point(CGAL::to_double((*ei).source().x()) + shiftX, CGAL::to_double((*ei).source().y()) + shiftY) << 
    				svg::Point(CGAL::to_double((*ei).target().x()) + shiftX, CGAL::to_double((*ei).target().y()) + shiftY));
    	}
    }
    // Add min cuts if want to
     for(auto it = graph->minCuts.begin(); it != graph->minCuts.end() && showMinCuts; ++it){
     	for(int i = 0; i < it->size(); i++){
    		doc << (svg::Polyline(svg::Stroke(cutEdgeSize, svg::Color::Brown)) << 
    				svg::Point(graph->edges[(*it)[i]]->start->x + shiftX, graph->edges[(*it)[i]]->start->y + shiftY) << 
    				svg::Point(graph->edges[(*it)[i]]->end->x + shiftX, graph->edges[(*it)[i]]->end->y + shiftY));
    	}
    	// Add vulnerable regions if want to
    	Rat_Polygon vulnerable = cut_to_vulnerable(graph, (*it), length);
    	for (auto ei = vulnerable.edges_begin(); ei != vulnerable.edges_end() && showVulnerableRegions; ++ei){
    		doc << (svg::Polyline(svg::Stroke(vulnerableSize, svg::Color::Red)) << 
    				svg::Point(CGAL::to_double((*ei).source().x()) + shiftX, CGAL::to_double((*ei).source().y()) + shiftY) << 
    				svg::Point(CGAL::to_double((*ei).target().x()) + shiftX, CGAL::to_double((*ei).target().y()) + shiftY));
    	}
    	// Add destruction regions if want to
    	Rat_Polygon destruction = cut_to_destruction(graph, (*it), length);
    	for (auto ei = destruction.edges_begin(); ei != destruction.edges_end() && showDestructionRegions; ++ei){
    		doc << (svg::Polyline(svg::Stroke(destructionSize, svg::Color::Orange)) << 
    				svg::Point(CGAL::to_double((*ei).source().x()) + shiftX, CGAL::to_double((*ei).source().y()) + shiftY) << 
    				svg::Point(CGAL::to_double((*ei).target().x()) + shiftX, CGAL::to_double((*ei).target().y()) + shiftY));
    	}
    }
    // Add disaster if want to
    if(showDisaster){
    	doc << (svg::Polyline(svg::Stroke(destructionLineSize, svg::Color::Red)) << 
    				svg::Point(dimX, shiftY/2.0) << svg::Point(dimX - length, shiftY/2.0));
    }

    doc.save();

}

void visibilityToSVG(string filename, Face_Handle &face, Poly_Arr visibility){
	// Start by finding the min and max values
	double xMin = DBL_MAX;
	double xMax = DBL_MIN;
	double yMin = DBL_MAX;
	double yMax = DBL_MIN;
	Arr_Edge_Circulator startCirc, circ, circ2;
	startCirc = face->outer_ccb();
	circ = startCirc;
	do{
		if(CGAL::to_double(circ->source()->point().x()) > xMax){xMax = CGAL::to_double(circ->source()->point().x());}
		if(CGAL::to_double(circ->source()->point().x()) < xMin){xMin = CGAL::to_double(circ->source()->point().x());}
		if(CGAL::to_double(circ->source()->point().y()) > yMax){yMax = CGAL::to_double(circ->source()->point().y());}
		if(CGAL::to_double(circ->source()->point().y()) < yMin){yMin = CGAL::to_double(circ->source()->point().y());}
		circ++;
	}while(startCirc != circ);

	// Draw the face dealing with
	svg::Dimensions dimensions(xMax - xMin, yMax - yMin);
    svg::Document doc(filename, svg::Layout(dimensions, svg::Layout::BottomLeft));
    // Add in the vertices and edges
    do{
    	doc << svg::Circle(svg::Point(CGAL::to_double(circ->source()->point().x()) - xMin,
    									CGAL::to_double(circ->source()->point().y()) - yMin), 
    						3, svg::Color::Black);

    	doc << (svg::Polyline(svg::Stroke(1, svg::Color::Black)) << 
    				svg::Point(CGAL::to_double(circ->source()->point().x()) - xMin,
    									CGAL::to_double(circ->source()->point().y()) - yMin) << 
    				svg::Point(CGAL::to_double(circ->target()->point().x()) - xMin,
    									CGAL::to_double(circ->target()->point().y()) - yMin));

    	circ++;
    }while(startCirc != circ);

    for(auto it = face->holes_begin(); it != face->holes_end(); ++it){
    	startCirc = *it;
    	circ = startCirc;
    	do{
    		doc << svg::Circle(svg::Point(CGAL::to_double(circ->source()->point().x()) - xMin,
    										CGAL::to_double(circ->source()->point().y()) - yMin), 
    							3, svg::Color::Black);

    		doc << (svg::Polyline(svg::Stroke(1, svg::Color::Black)) << 
    					svg::Point(CGAL::to_double(circ->source()->point().x()) - xMin,
    										CGAL::to_double(circ->source()->point().y()) - yMin) << 
    					svg::Point(CGAL::to_double(circ->target()->point().x()) - xMin,
    										CGAL::to_double(circ->target()->point().y()) - yMin));

    		circ++;
    	}while(startCirc != circ);

    }

    // Ok now do the same but with he visibility arrangement
    for(auto it = visibility.edges_begin(); it != visibility.edges_end(); ++it){

    	doc << svg::Circle(svg::Point(CGAL::to_double(it->source()->point().x()) - xMin,
    									CGAL::to_double(it->source()->point().y()) - yMin), 
    						3, svg::Color::Red);

    	doc << (svg::Polyline(svg::Stroke(1, svg::Color::Red)) << 
    				svg::Point(CGAL::to_double(it->source()->point().x()) - xMin,
    									CGAL::to_double(it->source()->point().y()) - yMin) << 
    				svg::Point(CGAL::to_double(it->target()->point().x()) - xMin,
    									CGAL::to_double(it->target()->point().y()) - yMin));
    }

    doc.save();

	
}







#endif

