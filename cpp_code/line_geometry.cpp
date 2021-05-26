#ifndef LINE_GEOMETRY_CPP
#define LINE_GEOMETRY_CPP

#include "line_geometry.h"

void addSegment(Graph* graph, list<Vertex> &edge){
    // It is possible two horizontal edges are added when one would do the same
    // that is with a new node in between. This creates a leaf with a 
    // non usable node which can be destroyed
    // Will simply remove a mid point in collinear edges if the mid point is new
    // as this looks to be a cgal error
    // for(auto it = edge.begin(); it != edge.end(); ++it){
    //     printf("%d ", it->id);
    // }
    // printf("\n");


    auto edgeIt = edge.begin();
    auto edgeIt2 = edgeIt;
    ++edgeIt2;
    auto edgeIt3 = edgeIt2;
    ++edgeIt3;
    while(edgeIt3 != edge.end()){
        // Check if collinear and node is new
        if(edgeIt2->id >= graph->vertices.size() && 
            orientation((*edgeIt), (*edgeIt2), (*edgeIt3)) == 0 &&
                onSegment((*edgeIt), (*edgeIt2), (*edgeIt3))){

            // Gotta remove useless node, move third iterator forward
            edgeIt2 = edge.erase(edgeIt2);
            ++edgeIt3;
        }else{
            ++edgeIt;
            ++edgeIt2;
            ++edgeIt3;
        } 
    }



    int currGraphSize = graph->vertices.size();
    // Add vertices and edges! Start with vertices, add if not in graph yet
    for(auto it = edge.begin(); it != edge.end(); ++it){
        if(it->id < currGraphSize){continue;}
        // New vertex, need to add it! Mark false as not usable
        // printf("Add %d (%d)\n", (int)graph->vertices.size(), it->id);
        it->id = graph->vertices.size();
        Vertex* newVertex = new Vertex(it->id, it->x, it->y, false);
        graph->vertices.push_back(newVertex);

        // Need to check new node being placed down is not on an existing edge.
        // If this is the case, need to break up existing edge
        for(int i = 0; i < graph->edges.size(); i++){
            if(graph->edges[i]->start->id == it->id || graph->edges[i]->end->id == it->id){continue;}
            if(orientation(*(graph->edges[i]->start), *(graph->edges[i]->end), (*it)) == 0 &&
                onSegment(*(graph->edges[i]->start), (*it), *(graph->edges[i]->end))){
                // Must break up old line
                // Add a new one
                // printf("Break %d %d by %d\n", graph->edges[i]->start->id, graph->edges[i]->end->id, it->id);
                Line* newLine = new Line(graph->edges.size(), graph->edges[i]->end, graph->vertices[it->id]);
                graph->edges.push_back(newLine);
                // Rewrite old one
                graph->edges[i]->end = graph->vertices[it->id];

            }
        }


    }
    // Now add edges, must check they do not already exist
    auto it = edge.begin();
    auto it2 = it;
    ++it2;
    while(it2 != edge.end()){
        // printf("Edge %d %d\n", (*it).id, (*it2).id);
        // Need to check this edge does not have an existing node as a midpoint.
        // If so, "break up" the edge into two
        bool breakEdge = false;
        for(int i = 0; i < graph->vertices.size() && !breakEdge; i++){
            if(graph->vertices[i]->id == it->id || graph->vertices[i]->id == it2->id){continue;}
            // Check existing node is not on the line, this is so if colinear and on the line
            if(orientation((*it), (*it2), (*graph->vertices[i])) == 0 &&
                onSegment((*it), (*graph->vertices[i]), (*it2))){
                // Must break up new line
                edge.insert(it2, (*graph->vertices[i]));
                breakEdge = true;
            }
        }
        if(breakEdge){
            --it2;
            continue;
        }


        // Need to check this edge does not already exist
        bool exists = false;
        // Go through all existing edges and check
        for(int i = 0; i < graph->edges.size() && !exists; i++){
            if((graph->edges[i]->start->id == it->id && graph->edges[i]->end->id == it2->id) ||
                (graph->edges[i]->start->id == it2->id && graph->edges[i]->end->id == it->id)){
                    exists = true;
                }
        }
        if(exists){
            // printf("Overlap of edges! Not horrible but worth mentioning!\n");
            ++it;
            ++it2; 
            continue;
        }

        Line* newLine;
        if(it->id < it2->id){newLine = new Line(graph->edges.size(), graph->vertices[it->id], graph->vertices[it2->id]);}
        else{ newLine = new Line(graph->edges.size(), graph->vertices[it2->id], graph->vertices[it->id]);}
        graph->edges.push_back(newLine);
        // printf("Add line (%d-%d)\n", graph->vertices[it2->id]->id, graph->vertices[it->id]->id);
        ++it;
        ++it2;
    }
}

int findVertexId(vector<Vertex*> &sortedVertices, Vertex* point){
    vector<Vertex*>::iterator iterator;
    // if(!std::binary_search(sortedVertices.begin(), sortedVertices.end(), point, compareVertices)){return -1;}
    // iterator = std::lower_bound (sortedVertices.begin(), sortedVertices.end(), point, compareVertices);
    if(!binary_search(sortedVertices.begin(), sortedVertices.end(), point, compareVertices)){return -1;}
    iterator = lower_bound (sortedVertices.begin(), sortedVertices.end(), point, compareVertices);
    if(iterator == sortedVertices.end()){return -1;}
    return (*iterator)->id;
}

bool compareVertices(Vertex* v1, Vertex* v2){
    return *v1 < *v2;
}

void makeAdjacencyMatrix(Graph* graph){

    vector< vector<int> > matrix(graph->vertices.size());

    graph->adjacencyMatrix = matrix;
    for(int i = 0; i < graph->vertices.size(); i++){
        graph->adjacencyMatrix[i].reserve(graph->vertices.size());
    }
    for(int i = 0; i < graph->edges.size(); i++){
        graph->adjacencyMatrix[graph->edges[i]->start->id].push_back(graph->edges[i]->end->id);
        graph->adjacencyMatrix[graph->edges[i]->end->id].push_back(graph->edges[i]->start->id);
    }
        
}

void resetGraph(Graph* graph){
    for(auto it = graph->cuts.begin(); it != graph->cuts.end(); ++it){
        it->erase(it->begin(), it->end());
    }
    graph->cuts.clear();
        
    for(auto it = graph->minCuts.begin(); it != graph->minCuts.end(); ++it){
        it->erase(it->begin(), it->end());
    }
    graph->minCuts.clear();

    for(auto it = graph->cutLabels.begin(); it != graph->cutLabels.end(); ++it){
        it->erase(it->begin(), it->end());
    }
    graph->cutLabels.clear();

    for(auto it = graph->minCutLabels.begin(); it != graph->minCutLabels.end(); ++it){
        it->erase(it->begin(), it->end());
    }
    graph->minCutLabels.clear();

    for(int i = 0; i < graph->leaves.size(); i++){
        graph->leaves[i].erase(graph->leaves[i].begin(), graph->leaves[i].end());
    }
    graph->leaves.erase(graph->leaves.begin(), graph->leaves.end());

    for(int i = 0; i < graph->leafMinCut.size(); i++){
        graph->leafMinCut[i].erase(graph->leafMinCut[i].begin(), graph->leafMinCut[i].end());
    }
    graph->leafMinCut.erase(graph->leafMinCut.begin(), graph->leafMinCut.end());

    vector<int> resetBlockIds(graph->vertices.size(), 0);
    graph->blockIds = resetBlockIds;

    for(int i = 0; i < graph->adjacencyMatrix.size(); i++){
        graph->adjacencyMatrix[i].erase(graph->adjacencyMatrix[i].begin(), graph->adjacencyMatrix[i].end());
    }
    graph->adjacencyMatrix.erase(graph->adjacencyMatrix.begin(), graph->adjacencyMatrix.end());
}


double distanceAtHeight(Line* line1, Line* line2, double height){
    double changeHeight1 = line1->start->y - line1->end->y;
    double changeHeight2 = line2->start->y - line2->end->y;
    double changeWidth1 = line1->start->x - line1->end->x;
    double changeWidth2 = line2->start->x - line2->end->x;
    double x1, x2, left1, left2, right1, right2;
    // If both lines are horizontal, only need to find closest points
    // Points will be right and left of different lines
    if(std::abs(changeHeight1) < TOLERANCE && std::abs(changeHeight2) < TOLERANCE){
        // Check left-right for both, return min
        left1 = std::min(line1->start->x,line1->end->x);
        left2 = std::min(line2->start->x,line2->end->x);
        right1 = std::max(line1->start->x,line1->end->x);
        right2 = std::max(line2->start->x,line2->end->x);
        return std::min(std::abs(left1 - right2), std::abs(left2 - right1));

    }
    // Again if one line is flat, work out the position of the other line
    // and then find closest point from horizontal line
    else if(std::abs(changeHeight1) < TOLERANCE){
        left1 = std::min(line1->start->x,line1->end->x);
        right1 = std::max(line1->start->x,line1->end->x);
        x2 = line2->end->x + (height - line2->end->y)*changeWidth2/changeHeight2;
        return std::min(std::abs(x2 - left1), std::abs(x2 - right1));
    
    }
    // Do the same if the other one is horizontal
    else if(std::abs(changeHeight2) < TOLERANCE){
        left2 = std::min(line2->start->x,line2->end->x);
        right2 = std::max(line2->start->x,line2->end->x);
        x1 = line1->end->x + (height - line1->end->y)*changeWidth1/changeHeight1;
        return std::min(std::abs(x1 - left2), std::abs(x1 - right2));
    }
    // Neither line is horizontal, easy thing.
    // Find position of both lines at this height and return.
    else{
        x1 = line1->end->x + (height - line1->end->y)*changeWidth1/changeHeight1;
        x2 = line2->end->x + (height - line2->end->y)*changeWidth2/changeHeight2;
        return std::abs(x1-x2);
    }
}   

bool onSegment(Vertex* p, Vertex* q, Vertex* r){
    return onSegment(*p, *q, *r);
}

bool onSegment(Vertex p, Vertex q, Vertex r){ 
    // Note this assumes the lines are colinear
    if (q.x <= max(p.x, r.x) && q.x >= min(p.x, r.x) && 
        q.y <= max(p.y, r.y) && q.y >= min(p.y, r.y)) 
       return true; 
    return false; 
} 

int orientation(Vertex* p, Vertex* q, Vertex* r){
    return orientation(*p, *q, *r);
}

int orientation(Vertex p, Vertex q, Vertex r){ 
    // See https://www.geeksforgeeks.org/orientation-3-ordered-points/ 
    // for details of below formula. 
    double val = (q.y - p.y) * (r.x - q.x) - 
              (q.x - p.x) * (r.y - q.y); 
    if (abs(val) < TOLERANCE) return 0;  // colinear 
  
    return (val > 0)? 1: 2; // clock or counterclock wise 
}

bool doIntersect(Vertex* p1, Vertex* q1, Vertex* p2, Vertex* q2){
    return doIntersect(*p1, *q1, *p2, *q2);
}

bool doIntersect(Vertex p1, Vertex q1, Vertex p2, Vertex q2){
     
    // Find the four orientations needed for general and 
    // special cases 
    int o1 = orientation(p1, q1, p2); 
    int o2 = orientation(p1, q1, q2); 
    int o3 = orientation(p2, q2, p1); 
    int o4 = orientation(p2, q2, q1); 
  
    // General case 
    if (o1 != o2 && o3 != o4) 
        return true; 
  
    // Special Cases 
    // p1, q1 and p2 are colinear and p2 lies on segment p1q1 
    if (o1 == 0 && onSegment(p1, p2, q1)) return true; 
  
    // p1, q1 and q2 are colinear and q2 lies on segment p1q1 
    if (o2 == 0 && onSegment(p1, q2, q1)) return true; 
  
    // p2, q2 and p1 are colinear and p1 lies on segment p2q2 
    if (o3 == 0 && onSegment(p2, p1, q2)) return true; 
  
     // p2, q2 and q1 are colinear and q1 lies on segment p2q2 
    if (o4 == 0 && onSegment(p2, q1, q2)) return true; 
  
    return false; // Doesn't fall in any of the above cases 
}  

double when_d_appart(Line* line1, Line* line2, double d){
    // Note it is assumed both lines are not horizontal
    double changeHeight1 = line1->start->y - line1->end->y;
    double changeHeight2 = line2->start->y - line2->end->y;
    double changeWidth1 = line1->start->x - line1->end->x;
    double changeWidth2 = line2->start->x - line2->end->x;
    // Assuming lines are not parallel, need to deal with that as well
    // If parallel, return NAN as error
    if(std::abs(changeWidth1/changeHeight1 - changeWidth2/changeHeight2) < TOLERANCE){return NAN;}

    // Note there are two possible heights for this exact value
    double height_1 = (d - line1->end->x + line1->end->y * changeWidth1/changeHeight1 + line2->end->x - line2->end->y * changeWidth2/changeHeight2)/
                        (changeWidth1/changeHeight1 - changeWidth2/changeHeight2);
    double height_2 = (-d - line1->end->x + line1->end->y * changeWidth1/changeHeight1 + line2->end->x - line2->end->y * changeWidth2/changeHeight2)/
                        (changeWidth1/changeHeight1 - changeWidth2/changeHeight2);
    
    // Both cannot be in both lines, as otherwise this means the lines cross
    bool height1_good, height2_good;
    // Check which one (if either) is possible
    // Start with the first one; if at level with start or end point cannot use it
    // since endpoint is protected
    if(std::abs(height_1 - line1->start->y) < TOLERANCE || 
        std::abs(height_1 - line1->end->y) < TOLERANCE || 
        std::abs(height_1 - line2->start->y) < TOLERANCE || 
        std::abs(height_1 - line2->end->y) < TOLERANCE){
        height1_good = false;

    }else if(height_1 > std::min(line1->start->y,line1->end->y) &&
        height_1 > std::min(line2->start->y,line2->end->y) &&
        height_1 < std::max(line1->start->y,line1->end->y) &&
        height_1 < std::max(line2->start->y,line2->end->y)){
        height1_good = true;
    }else{
        height1_good = false;
    }
    // Repeat with second one
    if(std::abs(height_2 - line1->start->y) < TOLERANCE || 
        std::abs(height_2 - line1->end->y) < TOLERANCE || 
        std::abs(height_2 - line2->start->y) < TOLERANCE || 
        std::abs(height_2 - line2->end->y) < TOLERANCE){
        height2_good = false;

    }else if(height_2 > std::min(line1->start->y,line1->end->y) &&
        height_2 > std::min(line2->start->y,line2->end->y) &&
        height_2 < std::max(line1->start->y,line1->end->y) &&
        height_2 < std::max(line2->start->y,line2->end->y)){
        height2_good = true;
    }else{
        height2_good = false;
    }

    if(height1_good && height2_good){
        printf("Cannot have both good, perhaps edges intersect?\n");
    }
    if(height1_good){
        return height_1;
    }else if(height2_good){
        return height_2;
    }else{
        return NAN;
    }
}

double edgeLength(Vertex* v1, Vertex* v2){
    return edgeLength(*v1, *v2);    
}

double edgeLength(Vertex v1, Vertex v2){
    double changeX = v1.x - v2.x;
    double changeY = v1.y - v2.y;
    return sqrt(changeX*changeX + changeY*changeY);    
}

Vertex extendEdge(Vertex source, Vertex target){
    Vertex extendedPoint;
    // Deal with vertical or horizontal edge
    if(edgeLength(source, target) < TOLERANCE){
        printf("Asking to extend a zero length edge!!\n");
        exit(0);
    }
    // Deal with vertical
    if(abs(source.x - target.x) < TOLERANCE){
        if(source.y > target.y){
            // Go slightly down
            extendedPoint.x = target.x;
            extendedPoint.y = target.y - TOLERANCE;
        }else{
            // Go slightly up
            extendedPoint.x = target.x;
            extendedPoint.y = target.y + TOLERANCE;
        }
    }
    // Deal with horizontal
    else if(abs(source.y - target.y) < TOLERANCE){
        if(source.x > target.x){
            // Go slightly left
            extendedPoint.x = target.x - TOLERANCE;
            extendedPoint.y = target.y;
        }else{
            // Go slightly up
            extendedPoint.x = target.x + TOLERANCE;
            extendedPoint.y = target.y;
        }
    }else{
        // Neither horizontal nor vertical. Find gradient
        double gradient = (target.y - source.y)/(target.x - source.x);
        // Work out in which direction to go
        double product = 1.0;
        if(source.x > target.x){
            product = -1.0;
        }
        extendedPoint.x = target.x + product * TOLERANCE;
        extendedPoint.y = target.y + product * TOLERANCE * gradient;
    }
    return extendedPoint;

}

Vertex moveInsideFace(Vertex source, Vertex middle, Vertex target){

    Vertex result;
    // Get the vectors
    double x2 = source.x - middle.x;
    double y2 = source.y - middle.y;
    double x1 = target.x - middle.x;
    double y1 = target.y - middle.y;
    
    // Get counterclockwise angle of the angles
    double dot = x1*x2 + y1*y2;      
    double det = x1*y2 - y1*x2;
    double angleVectors = atan2(det, dot);
    if(angleVectors < 0){angleVectors += 2*M_PI;}

    dot = 1*x1;
    det = y1;
    double angleBase = atan2(det, dot);
    if(angleBase < 0){angleBase += 2*M_PI;}

    double newAngle = angleBase + angleVectors/2;

    // Finally work out change 
    double changeX = TOLERANCE * 0.001 * cos(newAngle);
    double changeY = TOLERANCE * 0.001 * sin(newAngle);
    result.x = middle.x + changeX;
    result.y = middle.y + changeY;

    return result;

}

bool intInVect(vector<int> &vector, int id){
    bool in = false;
    for(int i = 0; i < vector.size() && !in; i++){
        in = (vector[i] == id);
    }
    return in;
}


 



#endif