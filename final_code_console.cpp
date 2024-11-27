#include <iostream>
#include <cmath>
#include <vector>
#include <bits/stdc++.h>
#include <random>
using namespace std;

//the modifiable attributes for handling the program
#define V 256 //total size of the universe
#define v2 16 //square root of the universe size
const int MAX_SIZE = 256; // Maximum size for the adjacency matrix
#define INF 99999 //infinity value for the graph
#define k 4 //range for the nearby drivers
#define graph_file "graph.txt" //path of the file containing the graph
#define driver_file "drivers.txt" //path of the file containing the driver locations

//class to manage the vEB tree operations
class vEBTree {
public:
    int u,u2; // Universe size
    int min, max; // Minimum and Maximum
    vEBTree* summary; // Summary for clusters
    vEBTree** cluster; // Pointer array for clusters

    // Constructor
    vEBTree(int size) {
        u = size; // Universe size
        u2 = std::sqrt(size); // Square root of universe size
        min = -1; //setting min and max of the tree to -1
        max = -1; 

        // If universe size is greater than 2, create summary and cluster
        if (size > 2) {
            int upperSqrt = upperSquareRoot(size); 
            summary = new vEBTree(upperSqrt);
            cluster = new vEBTree*[upperSqrt];
            for (int i = 0; i < upperSqrt; ++i) {
                cluster[i] = new vEBTree(upperSquareRoot(size));
            }
        } else {
            summary = nullptr;
            cluster = nullptr;
        }
    }

    // Destructor to free up memory
    ~vEBTree() {
        if (cluster) {
            int upperSqrt = upperSquareRoot(u);
            for (int i = 0; i < upperSqrt; ++i) {
                delete cluster[i];
            }
            delete[] cluster;
        }
        delete summary;
    }

    // Helper Functions
    int high(int x) { return x / lowerSquareRoot(u); }
    int low(int x) { return x % lowerSquareRoot(u); }
    int index(int high, int low) { return high * lowerSquareRoot(u) + low; }
    int lowerSquareRoot(int x) { return (int)pow(2, floor(log2(x) / 2)); }
    int upperSquareRoot(int x) { return (int)pow(2, ceil(log2(x) / 2)); }

    // Insert Function for the vEB Tree. We use Lazy Propagation for insertions
    void insert(int x) {
        // If tree is empty, set min and max to x and we are done
        // this is done for lazy propagation
        if (min == -1) {
            min = max = x;
        } 
        // If the tree is not empty, we need to insert the element
        else {
            // if x is less than min, then we set min as x and insert the previously min value in the tree
            if (x < min) swap(x, min);
            if (u > 2) {
                // If cluster is empty, insert x in summary and cluster
                if (cluster[high(x)]->min == -1) {
                    summary->insert(high(x));
                    cluster[high(x)]->insert(low(x));
                } 
                // If cluster is not empty, just insert x in cluster                
                else {
                    cluster[high(x)]->insert(low(x));
                }
            }
            //set the max value of the tree
            if (x > max) max = x;
        }
    }

    // Insert function which takes a pair of integers which are the coordinates of the point
    void insertpair(int x, int y) {
        int key = x * u2 + y; // Adjusted key generation formula for universe size 256 (16x16 grid)
        insert(key);
    }

    // Delete Function
    void remove(int x) {
        // If the structure contains only one element, reset both min and max to -1
        if (min == max) {
            min = max = -1;
        }
        // If u (universe size) is 2, only two elements are possible (0 and 1)
        else if (u == 2) {
            // Update min and max based on the removed element
            if (x == 0) min = 1;
            else min = 0;
            max = min;
        } 
        else {
            // If the element to be removed is the current min, find the next minimum
            if (x == min) {
                int firstCluster = summary->min;  // Find the first non-empty cluster
                x = index(firstCluster, cluster[firstCluster]->min);  // Calculate the new min element
                min = x;
            }
            
            // Remove the element from its cluster
            cluster[high(x)]->remove(low(x));
            
            // If the cluster becomes empty after removal
            if (cluster[high(x)]->min == -1) {
                summary->remove(high(x));  // Remove the cluster from the summary
                
                // If the removed element was also the max, update the max
                if (x == max) {
                    int summaryMax = summary->max;
                    if (summaryMax == -1) {  // If all clusters are empty, set max to min
                        max = min;
                    } else {
                        max = index(summaryMax, cluster[summaryMax]->max);  // Update max to the largest element in non-empty clusters
                    }
                }
            } 
            // If the removed element was max but the cluster isn't empty, update max within the cluster
            else if (x == max) {
                max = index(high(x), cluster[high(x)]->max);
            }
        }
    }

    // Successor Function
    int successor(int x) {
        // For a universe size of 2, check if x has a successor
        if (u == 2) {
            if (x == 0 && max == 1) return 1;  // If x is 0 and max is 1, successor is 1
            else return -1;  // No successor available
        }
        // If x is less than the current minimum, return the minimum as the successor
        else if (min != -1 && x < min) {
            return min;
        }
        else {
            // Get the maximum element in the cluster where x is located
            int maxLow = cluster[high(x)]->max;
            
            // If there is a successor within the same cluster
            if (maxLow != -1 && low(x) < maxLow) {
                int offset = cluster[high(x)]->successor(low(x));  // Find successor in the current cluster
                return index(high(x), offset);  // Return the full index of the successor
            }
            else {
                // Find the next cluster with elements
                int succCluster = summary->successor(high(x));
                
                // If there's no cluster with a successor, return -1
                if (succCluster == -1) return -1;
                else {
                    // Find the smallest element in the successor cluster
                    int offset = cluster[succCluster]->min;
                    return index(succCluster, offset);  // Return the full index of the successor
                }
            }
        }
    }

    // function to update the location of a driver
    void update_location(int x, int y, int new_x, int new_y) {
        int key = x * u2 + y;
        remove(key);
        insert(new_x * u2 + new_y);
    }

};


// Function to Run Dijkstra's Algorithm on the Graph and return the distance and parent vectors
std::pair<std::vector<int>, std::vector<int>> dijkstra(int graph[V][V], int src) {
    // Priority queue to select the minimum distance vertex.
    priority_queue<pair<int, int>, vector<pair<int, int>>, greater<pair<int, int>>> pq;

    vector<int> dist(V, INF); // Initialize all distances to INF
    vector<int> parents(V, -1); // Track the parent of each node

    // Start with the source node, with a distance of 0
    pq.push(make_pair(0, src));
    dist[src] = 0;

    // Process the queue until all reachable nodes have been handled
    while (!pq.empty()) {
        // Extract the node with the smallest distance
        int u = pq.top().second;
        pq.pop();

        // Explore all adjacent vertices v of node u
        for (int v = 0; v < V; v++) {
            // Check if there's an edge and if we found a shorter path to v
            if (graph[u][v] != INF && dist[u] != INF && dist[u] + graph[u][v] < dist[v]) {
                dist[v] = dist[u] + graph[u][v];
                parents[v] = u; // Set the parent of v to u
                pq.push(make_pair(dist[v], v));
            }
        }
    }

    // Return the distance and parent vectors
    return make_pair(dist, parents);
}

// Function to get the shortest path from the parent vector
std::vector<int> Get_Shortest_Path(std::vector<int> parents, int destVertex) {
    std::vector<int> path;
    if (parents[destVertex] == -1) return path;

    int current = destVertex;
    // Traverse the parent vector to get the path
    while (current != -1) {
        path.push_back(current);
        current = parents[current];
    }
    // Reverse the path to get the correct order
    reverse(path.begin(), path.end());
    return path;
}

// Function to calculate the fare based on distance and number of drivers
int fare(int dist, int num_drivers){
    // Base fare of 50, fare per km is 3, and fare reduction factor for available drivers is 0.035
    return 50 + 3*dist - 0.035*dist*num_drivers;
}

// Function to find non-empty points within a given range k
std::vector<std::pair<int, int>> findNearbyDrivers(int x, int y, vEBTree* tree1) {
    std::vector<std::pair<int, int>> nonEmptyPoints;
    int u = tree1->u2;
    // Define the boundaries of the square.
    int rowLow = std::max(0, x - k);
    int rowHigh = std::min(u-1, x + k); // For a 16x16 grid, row index is [0, 15]

    // Iterate through each row in the range [rowLow, rowHigh].
    for (int row = rowLow; row <= rowHigh; ++row) {
        // Define the column range for the square.
        int colLow = std::max(0, y - k);
        int colHigh = std::min(u-1, y + k); // Column index is also [0, 15] in a 16x16 grid

        // Convert the lower bound of the range to a key and search for the first successor.
        int startKey = row * u + colLow;
        int col = tree1->successor(startKey - 1);  // Start looking after the previous key

        // Iterate through all successors within the column range [colLow, colHigh]
        while (col != -1 && col / u == row && col % u <= colHigh) {
            if (col/u!=x || col%u!=y) nonEmptyPoints.push_back({col / u, col % u});
            col = tree1->successor(col);  // Move to the next element in the tree
        }
    }
    return nonEmptyPoints;
}   

// Main function to run the Taxi Management System
int main(){
    int universe_size = V; // Universe size for a 16x16 grid
    vEBTree tree(universe_size); // Insertion and queries row-wise
    cout<<"\n------------------------ TAXI MANAGEMENT SYSTEM ------------------------\n\n";
    
    // Input the Pick-Up and Destination Coordinates
    cout<<"Enter the Pick-Up Coordinates (x 'space' y): ";
    int x,y;
    cin>>x>>y;
    if (x<0 || x>=v2 || y<0 || y>=v2){
        cout<<"Invalid Coordinates\n";
        return 0;
    }

    cout<<"Enter the Destination Coordinates (x 'space' y): ";
    int dest_x,dest_y;
    cin>>dest_x>>dest_y;
    if (dest_x<0 || dest_x>=v2 || dest_y<0 || dest_y>=v2){
        cout<<"Invalid Coordinates\n";
        return 0;
    }

    // Insert the driver locations into the vEB Tree by reading data from the drivers file
    fstream file1(driver_file);
    if (!file1.is_open()) {
        cout << "Error opening drivers file" << endl;
        return 1;
    }
    int d1,d2;
    while (file1 >> d1 >> d2) {
        tree.insertpair(d1, d2);
    }

    // Find the nearby drivers within the range k from the source location
    vector<pair<int, int>> drivers = findNearbyDrivers(x,y, &tree);

    // If no drivers are available, return
    if (drivers.size()==0){
        cout<<"No Drivers available in your area. Please try again later.\n";
        return 0;
    }

    // Generate the graph from the graph file
    int graph[V][V];
    for (int i = 0; i < V; i++) {
        for (int j = 0; j < V; j++) {
            graph[i][j] = INF;
        }
    }
    fstream file2(graph_file);
    if (!file2.is_open()) {
        cout << "Error opening graph file" << endl;
        return 1;
    }
    
    // Read the graph data from the file
    int u, v, w;
    while (file2 >> u >> v >> w) {
        graph[u][v] = w;
        graph[v][u] = w;
    }
    
    // Run Dijkstra's Algorithm to find the shortest path and calculate the fare
    pair<vector<int>, vector<int>> dijkstra_result = dijkstra(graph, x*v2+y);
    vector<int> dist = dijkstra_result.first;
    vector<int> parents = dijkstra_result.second;

    // Sort the drivers based on their distance from the source
    vector <pair<int,int>> distances;
    for (const auto& point : drivers){
        distances.push_back({dist[point.first*v2+point.second],point.first*v2+point.second});
    }
    sort(distances.begin(),distances.end());

    // Display the list of drivers nearby
    cout<<"\nList of Drivers available nearby to pick a ride:\n";
    for (const auto& point : distances){
        cout<<"Location: ("<<point.second/v2<<","<<point.second%v2<<")    Distance: "<<point.first<<endl;
    }
    
    // Calculate the fare and ETA for the ride
    vector<int> shortestPath = Get_Shortest_Path(parents, dest_x*v2+dest_y);
    cout<<"\nFare Calculated: Rs. "<<fare(dist[dest_x*v2+dest_y],drivers.size())<<endl;
    cout<<"ETA: "<<dist[dest_x*v2+dest_y]<<" minutes"<<endl;

    // Display the shortest path from source to destination
    if (shortestPath.size()>0){
        cout<<"\nShortest Path from Source to Destination:\n";
        cout<<"("<<shortestPath[0]/v2<<","<<shortestPath[0]%v2<<")";
        for (int i=1;i<shortestPath.size();i++){
            cout<<" -> "<<"("<<shortestPath[i]/v2<<","<<shortestPath[i]%v2<<")";
        }
        cout<<endl;   
    }
    cout<<endl;
    return 0;
}
