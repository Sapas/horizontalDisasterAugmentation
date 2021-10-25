#ifndef CGAL_INTERFACE_CPP
#define CGAL_INTERFACE_CPP

#include "cgal_interface.h"

Rat_Polygon edge_minkowski_sum(Line* edge, double length){
	Rat_Polygon parallelogram;
	// Two options, either flat or not. Deal with flat first
	if(std::abs(edge->start->y - edge->end->y) < TOLERANCE){
		parallelogram.push_back(Rat_Point(std::min(edge->start->x, edge->start->x) - length/2.0, edge->start->y));
		parallelogram.push_back(Rat_Point(std::min(edge->start->x, edge->start->x) + length/2.0, edge->start->y));
	}else{
		// Not flat, add four points
		parallelogram.push_back(Rat_Point(edge->start->x - length/2.0, edge->start->y));
		parallelogram.push_back(Rat_Point(edge->start->x + length/2.0, edge->start->y));
		parallelogram.push_back(Rat_Point(edge->end->x + length/2.0, edge->end->y));
		parallelogram.push_back(Rat_Point(edge->end->x - length/2.0, edge->end->y));
	}
	return parallelogram;
}

Rat_Polygon cut_to_vulnerable(Graph* graph, vector<int> &cut, double length){
	// Set up initial polygon
	Rat_Polygon currentVulnerableRegion = edge_minkowski_sum(graph->edges[cut[0]], length);
	if(currentVulnerableRegion.orientation() == CGAL::COLLINEAR){return cut_to_vulnerable_horizontal(graph, cut, length, 0);}
	if(currentVulnerableRegion.orientation() != CGAL::CLOCKWISE){currentVulnerableRegion.reverse_orientation();}
	// For each extra edge, do Minkowski Sum and find intersection
	std::list<Rat_Polygon_Holes> polyIterator;
	for(int i = 1; i < cut.size(); i++){
		Rat_Polygon parallelogram = edge_minkowski_sum(graph->edges[cut[i]], length);
		// Intersections must be clockwise for intersection
		if(parallelogram.orientation() == CGAL::COLLINEAR){return cut_to_vulnerable_horizontal(graph, cut, length, i);}
		if(parallelogram.orientation() != CGAL::CLOCKWISE){parallelogram.reverse_orientation();}
		// Find intersection
		CGAL::intersection(currentVulnerableRegion, parallelogram, std::back_inserter(polyIterator));
		// Should only get one polygon back as should be dealing with convex polygons
		if(polyIterator.size() != 1){
			printf("Got multiple polygons from intersections of parallelograms %d (convex). Weird. Stopping.\n", (int)polyIterator.size());
			return currentVulnerableRegion;
		}
		currentVulnerableRegion = polyIterator.begin()->outer_boundary();
		// Empty list
		polyIterator.clear();
		if(currentVulnerableRegion.orientation() != CGAL::CLOCKWISE){currentVulnerableRegion.reverse_orientation();}	
	}
	return currentVulnerableRegion;
}

Rat_Polygon cut_to_vulnerable_horizontal(Graph* graph, vector<int> &cut, double length, int flatIndex){
	Rat_Polygon vulnearbleZone;
	// First of all need to know the height of the horizontal edge
	double height = graph->edges[cut[flatIndex]]->start->y;
	double xMin = std::min(graph->edges[cut[flatIndex]]->start->x, graph->edges[cut[flatIndex]]->end->x) - length/2.0;
	double xMax = std::max(graph->edges[cut[flatIndex]]->start->x, graph->edges[cut[flatIndex]]->end->x) + length/2.0;
	// Since the intersection will end up being a single line, need to find the line
	// relevant to each of the edges
	// Cycle through the edges and find this information
	for(int i = 0; i < cut.size(); i++){
		double xCurrentMin;
		double xCurrentMax;
		if(i == flatIndex){continue;}
		// First check the edge is not horizontal
		if(std::abs(graph->edges[cut[i]]->start->y - graph->edges[cut[i]]->end->y) < TOLERANCE){
			// Edge is horizontal, first check the height is the same as the other flat edge
			if(std::abs(graph->edges[cut[i]]->start->y - graph->edges[cut[flatIndex]]->start->y) > TOLERANCE){
				printf("Warning: Wrong behaviour, input to cut_to_culnerable_horizontal is not a cut!");
				return vulnearbleZone;
			}
			xCurrentMin = std::min(graph->edges[cut[i]]->start->x, graph->edges[cut[i]]->end->x) - length/2.0;
			xCurrentMax = std::max(graph->edges[cut[i]]->start->x, graph->edges[cut[i]]->end->x) + length/2.0;
		}else{
			// Not horizontal, work out the x position at this height
			double changeHeight = graph->edges[cut[i]]->start->y - graph->edges[cut[i]]->end->y;
    		double changeWidth = graph->edges[cut[i]]->start->x - graph->edges[cut[i]]->end->x;
    		double currentX = graph->edges[cut[i]]->end->x + (height - graph->edges[cut[i]]->end->y)*changeWidth/changeHeight;
    		xCurrentMin = currentX - length/2.0;
			xCurrentMax = currentX + length/2.0;
		}
		if(xMin < xCurrentMin){xMin = xCurrentMin;}
		if(xMax > xCurrentMax){xMax = xCurrentMax;}
	}
	if(xMin > xMax){
		printf("Warning: Wrong behaviour, found line of negative length (%.2f - %.2f)\n", xMin, xMax);
		for(int i = 0; i < cut.size(); i++){
			printf("Edge (%d,%d)\n", graph->edges[cut[i]]->start->id, graph->edges[cut[i]]->end->id);
		}
		return vulnearbleZone;
	}
	vulnearbleZone.push_back(Rat_Point(xMin, height));
	vulnearbleZone.push_back(Rat_Point(xMax, height));
	return vulnearbleZone;
}

Rat_Polygon vulnerable_to_destruction(Rat_Polygon vulnerable_zone, double length){
	Rat_Polygon destruction_zone;
	// First of all check vulnerable zone is not just a line; possible if cut 
	// contains a horizontal edge
	if(vulnerable_zone.orientation() == CGAL::COLLINEAR){
		// Flat, so only need extend horizontally
		auto iterator = vulnerable_zone.vertices_begin();
		double height = CGAL::to_double((*iterator).y());
		double xMin = CGAL::to_double((*iterator).x());
		double xMax = CGAL::to_double((*iterator).x());
		++iterator;
		int count = 0;
		while(iterator != vulnerable_zone.vertices_end()){
			if(std::abs(height - CGAL::to_double((*iterator).y())) > TOLERANCE){
				printf("Warning, assumed vulnearbleZone is horizontal but this is incorrect!\n");
				return destruction_zone;
			}
			if(CGAL::to_double((*iterator).x()) > xMax){xMax = CGAL::to_double((*iterator).x());}
			if(CGAL::to_double((*iterator).x()) < xMin){xMin = CGAL::to_double((*iterator).x());}
			count++;
			++iterator;
		}
		if(count != 1){
			printf("Warning, assumed vulnearbleZone is a line but this is incorrect\n");
			return destruction_zone;
		}
		destruction_zone.push_back(Rat_Point(xMin - length/2.0, height));
		destruction_zone.push_back(Rat_Point(xMax + length/2.0, height));
		return destruction_zone;

		
	}
	// Minkowski Sum requires both polygons to be counter clockwise
	if(vulnerable_zone.orientation() != CGAL::COUNTERCLOCKWISE){vulnerable_zone.reverse_orientation();}
	// Create disaster polygon
	Rat_Polygon disaster;
    disaster.push_back(Rat_Point(-length/2.0,0));
    disaster.push_back(Rat_Point(length/2.0,0));
    // Find Minkowski Sum
    destruction_zone = CGAL::minkowski_sum_by_full_convolution_2(disaster, vulnerable_zone).outer_boundary();
    return destruction_zone;
}

Rat_Polygon cut_to_destruction(Graph* graph, vector<int> &cut, double length){
	// First of all need to find vulnerable zone
	Rat_Polygon vulnerable_zone = cut_to_vulnerable(graph, cut, length);
	return vulnerable_to_destruction(vulnerable_zone, length);
}

Poly_Arr destructionArrangement(Graph* graph, double length, char searchType, int label1, int label2, vector<pair<int, int> > horizontal_forbidden){
	Poly_Arr resultArr;
	Overlay_traits overlay_traits;
	// Depending on the type, need to do different things.
	if(searchType == 'l'){
		// Leaf to something else
		// Label "forbidden" edges in the cut if horizontal
		for(int i = 0; i < graph->leafMinCut[label1].size(); i++){
			if(abs(graph->edges[graph->leafMinCut[label1][i]]->start->y - graph->edges[graph->leafMinCut[label1][i]]->end->y) < TOLERANCE){
				horizontal_forbidden.push_back(make_pair(graph->edges[graph->leafMinCut[label1][i]]->start->id, graph->edges[graph->leafMinCut[label1][i]]->end->id));
			}
		}
		// if(horizontal_forbidden.size() > 0){
		// 	for(int i = 0; i < horizontal_forbidden.size(); i++){
		// 		printf("Forbidden single %d - %d\n", horizontal_forbidden[i].first, horizontal_forbidden[i].second);
		// 	}
		// }
		// Destruction zone of a single min cut. Get it and turn it into an arrangement
		Rat_Polygon destruction = cut_to_destruction(graph, graph->leafMinCut[label1], length);
		for (auto ei = destruction.edges_begin(); ei != destruction.edges_end(); ++ei){
			// insert_non_intersecting_curve(resultArr, *ei);
			insert(resultArr, *ei);
 		}
 		// Label the outside face as the only valid one
 		CGAL_assertion (resultArr.number_of_faces() == 2);
  		for (auto fit = resultArr.faces_begin(); fit != resultArr.faces_end(); ++fit){
  			fit->set_data(fit == resultArr.unbounded_face());
  		}
  		return resultArr;

	}else if(searchType == 'p'){
		// Pair of leaves
		// Need two min cuts
		Poly_Arr arr1, arr2;
		// Label "forbidden" edges in the cut if horizontal
		for(int i = 0; i < graph->leafMinCut[label1].size(); i++){
			if(abs(graph->edges[graph->leafMinCut[label1][i]]->start->y - graph->edges[graph->leafMinCut[label1][i]]->end->y) < TOLERANCE){
				horizontal_forbidden.push_back(make_pair(graph->edges[graph->leafMinCut[label1][i]]->start->id, graph->edges[graph->leafMinCut[label1][i]]->end->id));
			}
		}
		for(int i = 0; i < graph->leafMinCut[label2].size(); i++){
			if(abs(graph->edges[graph->leafMinCut[label2][i]]->start->y - graph->edges[graph->leafMinCut[label2][i]]->end->y) < TOLERANCE){
				horizontal_forbidden.push_back(make_pair(graph->edges[graph->leafMinCut[label2][i]]->start->id, graph->edges[graph->leafMinCut[label2][i]]->end->id));
			}
		}
		// if(horizontal_forbidden.size() > 0){
		// 	for(int i = 0; i < horizontal_forbidden.size(); i++){
		// 		printf("Forbidden double %d - %d\n", horizontal_forbidden[i].first, horizontal_forbidden[i].second);
		// 	}
		// }
		
		Rat_Polygon destruction = cut_to_destruction(graph, graph->leafMinCut[label1], length);
		for (auto ei = destruction.edges_begin(); ei != destruction.edges_end(); ++ei){
			insert(arr1, *ei);
 		}
 		Rat_Polygon destruction2 = cut_to_destruction(graph, graph->leafMinCut[label2], length);
		for (auto ei = destruction2.edges_begin(); ei != destruction2.edges_end(); ++ei){
			insert(arr2, *ei);
 		}
 		// Assign outside face as the only valid one
 		CGAL_assertion (arr1.number_of_faces() == 2);
  		for (auto fit = arr1.faces_begin(); fit != arr1.faces_end(); ++fit){
  			fit->set_data(fit == arr1.unbounded_face());
  		}
  		CGAL_assertion (arr2.number_of_faces() == 2);
  		for (auto fit = arr2.faces_begin(); fit != arr2.faces_end(); ++fit){
  			fit->set_data(fit == arr2.unbounded_face());
  		}
 		// Overlay the two
 		overlay(arr1, arr2, resultArr, overlay_traits);
 		return resultArr;

	}else if(searchType == 'n'){
		// Need to look for every cut that separates the two blocks
		// Need a node from each of the blocks to "ask questions"
		int node1 = -1, node2 = -1;
		for(int i = 0; i < graph->blockIds.size() && (node1 == -1 || node2 == -1); i++){
			if(graph->blockIds[i] == label1){node1 = i;}
			if(graph->blockIds[i] == label2){node2 = i;}
		}

		Poly_Arr currentDestruction;
		bool first = true;
		auto itCut = graph->minCuts.begin();
		auto itCutLabels = graph->minCutLabels.begin();
		while(itCut != graph->minCuts.end() && itCutLabels != graph->minCutLabels.end()){
			Poly_Arr newDestruction, tempDestruction;
			// If this is a cut which separates the two blocks, new edge cannot
			// intersect destruction region
			if((*itCutLabels)[node1] != (*itCutLabels)[node2]){
				// printf("Looking at cut ");
				// for(int i = 0; i < itCut->size(); i++){
				// 	printf("(%d, %d),  ", graph->edges[(*itCut)[i]]->start->id, graph->edges[(*itCut)[i]]->end->id);
				// }
				// printf("\n");
				// Label "forbidden" edges in the cut if horizontal
				for(int i = 0; i < (*itCut).size(); i++){
					if(abs(graph->edges[(*itCut)[i]]->start->y - graph->edges[(*itCut)[i]]->end->y) < TOLERANCE){
						horizontal_forbidden.push_back(make_pair(graph->edges[(*itCut)[i]]->start->id, graph->edges[(*itCut)[i]]->end->id));
					}
				}
				Rat_Polygon destruction = cut_to_destruction(graph, (*itCut), length);
				for (auto ei = destruction.edges_begin(); ei != destruction.edges_end(); ++ei){
					insert(newDestruction, *ei);
 				}
 				CGAL_assertion (newDestruction.number_of_faces() == 2);
  				for (auto fit = newDestruction.faces_begin(); fit != newDestruction.faces_end(); ++fit){
  					fit->set_data(fit == newDestruction.unbounded_face());
  				}
  				if(first){
  					first = false;
  					currentDestruction = newDestruction;
  				}else{
  					overlay(currentDestruction, newDestruction, tempDestruction, overlay_traits);
  					currentDestruction = tempDestruction;
  				}
			}
			++itCut;
			++itCutLabels;
		}
		return currentDestruction;
	}else{
		printf("Usage of destruction zone function wrong!\n");
		return resultArr;
	}	
}

void makeVisibilityGraph(Graph* graph, Face_Handle &face, TEV &tev, vector<Poly_Arr::Vertex_handle> &faceVertices, vector<list< pair<int,double> > > &adjacentMatrix, vector<pair<int, int> > horizontal_forbidden){
	Arr_Edge_Circulator startCirc, circ, circ2;
	// Looks like getting repeated vertices, must avoid this 
	// Need to have condition here for outside face
	if(face->holes_begin() == face->holes_end()){
		// No holes, so use outer boundary
		startCirc = face->outer_ccb();
	}else{
		// Hole, so this should be real outside face, and hole is outer boundary
		startCirc = (*face->holes_begin());
	}
	circ = startCirc;
	// Go once and reset labels from faces
	do{
		circ->target()->set_data(make_pair(circ->target()->data().first, -1));
		circ++;
	}while(startCirc != circ);
	
	// Want to store points slightly inside the face to differentiate
	// between points that are on either side of an edge
	// Store the Nodes, the points and the id of the point (as some can be repeated)
	int idCount = 0;
	// Find all vertices on the face
	faceVertices.reserve(graph->vertices.size());
	vector<Vertex> viewPoints;
	vector<int> viewPointsIds;
	viewPoints.reserve(graph->vertices.size());
	viewPointsIds.reserve(graph->vertices.size());
	int currentId;
	do{
		// Ignore zero edge length edge
		// Can't just check ids, need more
		Vertex start, end;
		start.x = CGAL::to_double(circ->source()->point().x());
		start.y = CGAL::to_double(circ->source()->point().y());
		end.x = CGAL::to_double(circ->target()->point().x());
		end.y = CGAL::to_double(circ->target()->point().y());
		if(edgeLength(start, end) < TOLERANCE){circ++; continue;}



		// if(circ->target()->data().first == circ->source()->data().first){circ++; continue;}
		// Need the next non zero length edge
		circ2 = circ;
		do{
			circ2++;
			start.x = CGAL::to_double(circ2->source()->point().x());
			start.y = CGAL::to_double(circ2->source()->point().y());
			end.x = CGAL::to_double(circ2->target()->point().x());
			end.y = CGAL::to_double(circ2->target()->point().y());
		}while(edgeLength(start, end) < TOLERANCE);

		// Sanity check
		// if(circ2->source()->data().first != circ->target()->data().first){printf("What you thought is wrong!\n");}
		// Sweet now want the viewpoint, slightly inside face

		// Make vertices out of the relevant points
		Vertex newViewPoint, source, middle, target;
		source.x = CGAL::to_double(circ->source()->point().x());
		source.y = CGAL::to_double(circ->source()->point().y());
		middle.x = CGAL::to_double(circ->target()->point().x());
		middle.y = CGAL::to_double(circ->target()->point().y());
		target.x = CGAL::to_double(circ2->target()->point().x());
		target.y = CGAL::to_double(circ2->target()->point().y());
		if(circ2->target()->data().first == circ->source()->data().first){newViewPoint = extendEdge(source, middle);}
		else{newViewPoint = moveInsideFace(source, middle, target);}

		// if(circ->target()->data().first == 29 || circ->target()->data().first == 30){
		// 	printf("%d (%.2f, %.2f), %d (%.2f, %.2f), %d(%.2f, %.2f)\n", circ->source()->data().first,
		// 							CGAL::to_double(circ->source()->point().x()),
		// 							CGAL::to_double(circ->source()->point().y()),
		// 							circ->target()->data().first,
		// 							CGAL::to_double(circ->target()->point().x()),
		// 							CGAL::to_double(circ->target()->point().y()),
		// 							circ2->target()->data().first,
		// 							CGAL::to_double(circ2->target()->point().x()),
		// 							CGAL::to_double(circ2->target()->point().y()));
		// }

		// Work out index
		if(circ->target()->data().second == -1){
			// Have not seen this yet, or it is not a usable node.
			// Either way need to treat it as a new node
			faceVertices.push_back(circ->target());

			// Must store this id if it is a usable node
			// Changed my mind, usable nodes cannot be used to to "go through" edges
			// if(circ->target()->data().first < graph->vertices.size() &&
			// 	graph->vertices[circ->target()->data().first]->usable){

			// 	circ->target()->set_data(make_pair(circ->target()->data().first, idCount));
			// }
			currentId = idCount;
			idCount++;
		}else{
			currentId = circ->target()->data().second;
		}

		// Have view point as well as relevant id, store
		viewPoints.push_back(newViewPoint);
		viewPointsIds.push_back(currentId);
		circ++;
	}while(startCirc != circ);

	// At this stage should just look at every pair once, but do both and make sure 
	// they can see each other (sanity check)
	// Save space for visibility
	vector< vector<bool> > canSee;
	canSee.reserve(faceVertices.size());
	adjacentMatrix.reserve(faceVertices.size());
	for(int i = 0; i < faceVertices.size(); i++){
		list< pair<int, double> > newList(0);
		vector<bool> vision(faceVertices.size(), false); 
		adjacentMatrix.push_back(newList);
		canSee.push_back(vision);
	}

	// Look at every pair
	for(int p1 = 0; p1 < viewPoints.size(); p1++){
		// Create visibility arrangement and searcher
		Poly_Arr visibility;
		Rat_Point queryPoint = Rat_Point(viewPoints[p1].x, viewPoints[p1].y);
		tev.compute_visibility(queryPoint, face, visibility);
		Point_Searcher point_searcher(visibility);

		// if(faceVertices[viewPointsIds[p1]]->data().first == 30 ||
		// 	faceVertices[viewPointsIds[p1]]->data().first == 29){printf("%d (%d) Real vs moved\n(%.10f, %.10f)\n(%.10f, %.10f)\n", p1, faceVertices[viewPointsIds[p1]]->data().first,
		// 	 CGAL::to_double(faceVertices[viewPointsIds[p1]]->point().x()),
		// 	 CGAL::to_double(faceVertices[viewPointsIds[p1]]->point().y()),
		// 	 viewPoints[p1].x, viewPoints[p1].y);

		// string filename = "../testRuns/visibility" + to_string(p1) + "-" + to_string(faceVertices[viewPointsIds[p1]]->data().first) + ".svg";
		// visibilityToSVG(filename, face, visibility);
		// }

		for(int p2 = 0; p2 < viewPoints.size(); p2++){
			if(faceVertices[viewPointsIds[p1]]->data().first == faceVertices[viewPointsIds[p2]]->data().first){continue;}
			// Here need to add that if both points are vertices of the graph, horizontal,
			// an edge already exists and it is part of a cut, these two cannot "see" each other
			if(faceVertices[viewPointsIds[p1]]->data().first < graph->vertices.size() &&
				faceVertices[viewPointsIds[p2]]->data().first < graph->vertices.size() &&
				abs(graph->vertices[faceVertices[viewPointsIds[p2]]->data().first]->y - graph->vertices[faceVertices[viewPointsIds[p1]]->data().first]->y) < TOLERANCE){

				// Cycle through "forbidden" edges, check if this one already exists
				bool forbidden = false;
				for(int i = 0; i < horizontal_forbidden.size(); i++){
					if((horizontal_forbidden[i].first == faceVertices[viewPointsIds[p1]]->data().first &&
						horizontal_forbidden[i].second == faceVertices[viewPointsIds[p2]]->data().first) || 
						(horizontal_forbidden[i].first == faceVertices[viewPointsIds[p2]]->data().first &&
						horizontal_forbidden[i].second == faceVertices[viewPointsIds[p1]]->data().first)){
						
						forbidden = true;
						break;

					}
				}
				if(forbidden){continue;}
			}
			
			Rat_Point lookingPoint = Rat_Point(viewPoints[p2].x, viewPoints[p2].y);
			Poly_Arr::Face_const_handle * foundFace;
  			CGAL::Arr_point_location_result<Poly_Arr>::Type obj = point_searcher.locate(lookingPoint);
			bool canSeeResult = !((foundFace = boost::get<Poly_Arr::Face_const_handle> (&obj)) && (*foundFace)->is_unbounded());
			
			// Another sanity check, same point should not see itself (from different view points)
			if(faceVertices[viewPointsIds[p1]]->data().first == faceVertices[viewPointsIds[p2]]->data().first &&
				canSeeResult){
				printf("Should not be happening!! %d\n(%.12f,%.12f)\n(%.12f,%.12f)\n", faceVertices[viewPointsIds[p1]]->data().first,
																					viewPoints[p1].x, viewPoints[p1].y, viewPoints[p2].x, viewPoints[p2].y);
				visibilityToSVG("../data/plots/visibility.svg", face, visibility);
				exit(0);
			}
			

			if(canSeeResult){
  				canSee[viewPointsIds[p1]][viewPointsIds[p2]] = true;
  			}
		}
	}

	// So now should go through each pair of the nodes, and if they can see each other save the point
	for(int i = 0; i < faceVertices.size(); i++){
		for(int j = i + 1; j < faceVertices.size(); j++){
			if((canSee[i][j] && !canSee[j][i]) || (!canSee[i][j] && canSee[j][i])){printf("Weird, seeing is not symmetric!! %d %d\n", faceVertices[i]->data().first, faceVertices[j]->data().first);}

			if(!canSee[i][j] || !canSee[j][i]){continue;}

			Vertex p1 = Vertex(-1, CGAL::to_double(faceVertices[i]->point().x()), CGAL::to_double(faceVertices[i]->point().y()));
  			Vertex p2 = Vertex(-1, CGAL::to_double(faceVertices[j]->point().x()), CGAL::to_double(faceVertices[j]->point().y()));
  			adjacentMatrix[i].push_back(make_pair(j, edgeLength(p1,p2)));
  			adjacentMatrix[j].push_back(make_pair(i, edgeLength(p1,p2)));
		}
	}
	// Now need to add the "free" edges which are edges between "available" nodes
	// REMOVING THIS AS I DO NOT THINK IT MAKES SENSE lOGIC WISE
	// do{
	// 	if(circ->target()->data().first != circ->source()->data().first &&
	// 		circ->target()->data().second != -1 &&
	// 		circ->source()->data().second != -1){
	// 		// This is a "free" edge, so add it!
	// 		adjacentMatrix[circ->target()->data().second].push_back(make_pair(circ->source()->data().second, 0.0));
	// 		adjacentMatrix[circ->source()->data().second].push_back(make_pair(circ->target()->data().second, 0.0));
	// 	}
	// 	circ++;
	// }while(startCirc != circ);

}

double dijkstraMultipleSourceGoal(vector<int> &startGoalMid, vector<list< pair<int,double> > > &adjacentMatrix, vector<int> &path){
	double best = DBL_MAX;
	int size = startGoalMid.size();
	// Use every available source
	for(int source = 0; source < size; source++){
		if(startGoalMid[source] != 1){continue;}
		// Define the data structures that will be needed 
		vector<double> cost(size, DBL_MAX);
		vector<int> predecessor(size, -1);
		predecessor[source] = source;
		cost[source] = 0;
		priority_queue< pair<int, double>, vector <pair< int, double> >, dijkstraCustomComparator > pq; 
		pq.push(make_pair(source, 0.0));
		int current;
		while(!pq.empty()){
			current = pq.top().first;
			// printf("%d cost %.4f pred %d - ", current, cost[current], predecessor[current]);
			// Check if this is one of the goals
			if(startGoalMid[current] == 2){break;}
			// Otherwise check neighbours and update
			for(auto it = adjacentMatrix[current].begin(); it != adjacentMatrix[current].end(); ++it){
				if(cost[it->first] > pq.top().second + it->second && abs(cost[it->first] - pq.top().second + it->second) > TOLERANCE){
					cost[it->first] = pq.top().second + it->second;
					predecessor[it->first] = pq.top().first;
					pq.push(make_pair(it->first, pq.top().second + it->second));
				}
			}
			pq.pop();
		}
		// printf("\n");
		if(startGoalMid[current] != 2){
			// Due to floating point error, can't find a path. Assume it does not exist
			continue;
		}
		// If what is found is worse than the best found so far, stop!
		if(cost[current] > best){continue;}
		best = cost[current];
		// Here current will be the last node, need to get the path back!
		path.clear();
		path.reserve(size);
		path.push_back(current);
		while(current != predecessor[current]){
			current = predecessor[current];
			path.push_back(current);
		}
		// Reverse
		reverse(path.begin(), path.end());
	}
	return best;
}

pair<list<Vertex>, double> findShortestEdge(Graph* graph, double length, char searchType, int label1, int label2, vector<Vertex*> &sortedVertices, Poly_Arr &graphArr){
	list<Vertex> placeholder;
	pair<list<Vertex>, double> bestSolution = make_pair(placeholder, DBL_MAX);
	vector<pair<int, int> > horizontal_forbidden;
	// Find the arrangement of the destruction zones and overlay it with the graph arrangement
  	Poly_Arr destructionArr = destructionArrangement(graph, length, searchType, label1, label2, horizontal_forbidden);
  	Poly_Arr resultArr;
  	Overlay_traits overlay_traits;
  	overlay(graphArr, destructionArr, resultArr, overlay_traits);
  	// Sweet! Now don't know how to keep track of ids of vertices, so re ad ids to all nodes
  	// in arrangement, and add id larger than existing vertices if not an original vertex of the graph
  	// Also keep track of poisitions of new vertices
  	vector<Vertex*> newVertices;
  	newVertices.reserve(graph->vertices.size());
  	int newIds = graph->vertices.size();
  	for(auto vit = resultArr.vertices_begin(); vit != resultArr.vertices_end(); ++vit){
  		Vertex* vertex = new Vertex(-1, CGAL::to_double(vit->point().x()), CGAL::to_double(vit->point().y()));
  		int result = findVertexId(sortedVertices, vertex);
  		if(result != -1){
  			delete vertex;
  			vit->set_data(make_pair(result, -1));
  		}
  		else{
  			vit->set_data(make_pair(newIds,-1));
  			vertex->id = newIds;
  			newVertices.push_back(vertex);
  			newIds++;
  		}	
  	}
  	// Want to associate point searcher to arrangement now for preprocessing
	TEV tev;
	tev.attach(resultArr);
	int specCount = 0;
  	// Ok so should be pretty sweet now, should be able to pass onto two functions to take care of the rest!
  	// Cycle through the faces and check whether an edge can be added here
  	bool foundOutsideFace = false;
  	for (Face_Handle fit = resultArr.faces_begin(); fit != resultArr.faces_end(); ++fit){
  		// Check this one is acceptable (i.e. does not intersect destruction zone)
  		if(!fit->data()){continue;}
		// Need to have at least two candidate endpoints for the edge
		bool vertexIn = false;
		bool vertexOut = false;
		Poly_Arr::Ccb_halfedge_circulator circStart, circ;
		// If have holes, this is actual outside face, so should use the whole
		// Should have a check to make sure this only happens once
		if(fit->holes_begin() == fit->holes_end()){
			// No holes, so use outer boundary
			circStart = fit->outer_ccb();
		}else{
			// Hole, so this should be real outside face, and hole is outer boundary
			// Check this has only happened once
			if(foundOutsideFace){printf("Looks like duplicate outside face, assumptions incorrect!!\n"); return bestSolution;}
			foundOutsideFace = true;
			circStart = (*fit->holes_begin());
		}
  		circ = circStart;
		do{
			if(circ->source()->data().first < graph->vertices.size()){
				if(searchType == 'l'){
					// Adding a leaf to something else, so need something in the leaf and something outside
					vertexIn = vertexIn || intInVect(graph->nodeLeaves[circ->source()->data().first], label1);
					vertexOut = vertexOut || (!intInVect(graph->nodeLeaves[circ->source()->data().first], label1) && graph->vertices[circ->source()->data().first]->usable);
				}else if(searchType == 'p'){
					// Adding an edge between two leaves, so need a node from both leaves
					vertexIn = vertexIn || intInVect(graph->nodeLeaves[circ->source()->data().first], label1);
					vertexOut = vertexOut || intInVect(graph->nodeLeaves[circ->source()->data().first], label2);
				}else if(searchType == 'n'){
					// Adding an edge between two blocks, so need node from the block which is usable
					vertexIn = vertexIn || (graph->blockIds[circ->source()->data().first] == label1 && graph->vertices[circ->source()->data().first]->usable);
					vertexOut = vertexOut || (graph->blockIds[circ->source()->data().first] == label2 && graph->vertices[circ->source()->data().first]->usable);
				}else{
					printf("This has not been implemented yet!\n");
				}
				
			}
			++circ;
		}while(circ != circStart && (!vertexIn || !vertexOut));
		if(!vertexIn || !vertexOut){continue;}
		// If got this far, know what face is interesting. Print to show it all works
		// printf("Face ");
		// circ = circStart;
		// do{
		// 	if(circ->source()->data().first >= graph->vertices.size()){printf("*");}
		// 	printf("%d ", circ->source()->data().first);
		// 	++circ;
		// }while(circ != circStart);
		// printf(" should be looked at\n");
		// Display face so I can see it
		// if(label1 == 0 && label2 == -1){
		// 	printf("Print!\n");
		// 	string filename = "../data/plots/visibility" + to_string(specCount) + ".svg";
		// 	Poly_Arr visibility;
		// 	visibilityToSVG(filename, fit, visibility);
		// }
		// Ok so now comes the point at which the visibility graph is needed
		vector<Poly_Arr::Vertex_handle> faceVertices;
		vector<list< pair<int,double> > > adjacentMatrix;
		makeVisibilityGraph(graph, fit, tev, faceVertices, adjacentMatrix, horizontal_forbidden);
		// printf("Vertices: ");
		// for(int i = 0; i < faceVertices.size(); i++){
		// 	printf("%d ", faceVertices[i]->data().first);
		// }
		// printf("\n");
		// printf("Edges:\n");
		// for(int i = 0; i < adjacentMatrix.size(); i++){
		// 	printf("%d: ", faceVertices[i]->data().first);
		// 	for(auto it = adjacentMatrix[i].begin(); it != adjacentMatrix[i].end(); ++it){
		// 		printf("%d (%d)(%.2f)", faceVertices[it->first]->data().first, faceVertices[it->first]->data().second, it->second);
		// 	}
		// 	printf("\n");
		// }
		// Ok! Have visibility graph, now want to find shortest connection between leaf and not leaf!
		// Create a vector which denotes a vertex 0 if not a start or goal, 1 if start, and 2 if goal
		vector<int> startGoalMid(faceVertices.size(), 0);
		for(int i = 0; i < faceVertices.size(); i++){
			if(searchType == 'l'){
				// Adding an edge from a leaf to outside, 
				// so a node in the leaf is a start, and any other node (which is usable) is an end
				if(faceVertices[i]->data().first >= graph->vertices.size() || !(graph->vertices[faceVertices[i]->data().first]->usable)){continue;}
				if(intInVect(graph->nodeLeaves[faceVertices[i]->data().first], label1)){startGoalMid[i] = 1;}
				else{startGoalMid[i] = 2;}
			}else if(searchType == 'p'){
				// Adding an edge between leaves, so a node in the first leaf
				// is a start, and a node in the second leaf is a goal
				if(faceVertices[i]->data().first >= graph->vertices.size() || !(graph->vertices[faceVertices[i]->data().first]->usable)){continue;}
				if(intInVect(graph->nodeLeaves[faceVertices[i]->data().first], label1) && 
					!intInVect(graph->nodeLeaves[faceVertices[i]->data().first], label2)){startGoalMid[i] = 1;}
				if(!intInVect(graph->nodeLeaves[faceVertices[i]->data().first], label1) && 
					intInVect(graph->nodeLeaves[faceVertices[i]->data().first], label2)){startGoalMid[i] = 2;}
			}else if(searchType == 'n'){
				// Adding an edge between two nodes, 
				if(faceVertices[i]->data().first >= graph->vertices.size() || !(graph->vertices[faceVertices[i]->data().first]->usable)){continue;}
				if(graph->blockIds[faceVertices[i]->data().first] == label1){startGoalMid[i] = 1;}
				if(graph->blockIds[faceVertices[i]->data().first] == label2){startGoalMid[i] = 2;}
			}else{
				printf("This has not been implemented yet!\n");
			}
			
		}
		// for(int i = 0; i < startGoalMid.size(); i++){
		// 	printf("%d ", startGoalMid[i]);
		// }
		// printf("\n");
		vector<int> path;
		double pathCost = dijkstraMultipleSourceGoal(startGoalMid, adjacentMatrix, path);
		
		// printf("Locations:\n");
		// for(int i = 0; i < faceVertices.size(); i++){
		// 	printf("%d (%.4f, %.4f)\n", i, CGAL::to_double(faceVertices[i]->point().x()), CGAL::to_double(faceVertices[i]->point().y()));
		// }
		// printf("Visibility:\n");
		// for(int i = 0; i < adjacentMatrix.size(); i++){
		// 	printf("%d: ", i);
		// 	for(auto it = adjacentMatrix[i].begin(); it != adjacentMatrix[i].end(); ++it){
		// 		printf("%d ", it->first);
		// 	}
		// 	printf("\n");
		// }
		// NEED TO SAVE IT SOMEHOW
		if(pathCost < bestSolution.second){
			// Make path of vertices
			list<Vertex> newPath;
			// Show it
			// for(int i = 0; i < path.size(); i++){
			// 	printf("%d ", faceVertices[path[i]]->data().first);
			// }
			// printf("\n");

			// printf("Found a better face, the best edge has cost %.2f and is: ", pathCost);
			for(int i = 0; i < path.size(); i++){
				// printf("%d ", faceVertices[path[i]]->data().first);
				// Check no repeated vertices
				if(i > 0){
					Vertex start, end;
					start.x = CGAL::to_double(faceVertices[path[i]]->point().x());
					start.y = CGAL::to_double(faceVertices[path[i]]->point().y());
					end.x = CGAL::to_double(faceVertices[path[i-1]]->point().x());
					end.y = CGAL::to_double(faceVertices[path[i-1]]->point().y());
					if(edgeLength(start, end) < TOLERANCE){circ++; continue;}
				}
				if(faceVertices[path[i]]->data().first >= graph->vertices.size()){
					newPath.push_back(*(newVertices[faceVertices[path[i]]->data().first - graph->vertices.size()]));
				}else{
					newPath.push_back(*(graph->vertices[faceVertices[path[i]]->data().first]));
				}
			}
			// printf("\n");
			bestSolution = make_pair(newPath, pathCost);
		}
		specCount++;
	}

	// Need to do something about all the stored data. Start by deleting all "newVertices"
	for(int i = 0; i < newVertices.size(); i++){delete newVertices[i];}

	return bestSolution;
}


#endif