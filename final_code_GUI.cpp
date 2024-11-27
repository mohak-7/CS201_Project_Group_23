//wxWidgets libraries
#include <wx/wx.h>
#include <wx/graphics.h>
#include <wx/valnum.h>

//standard C++ libraries
#include <bits/stdc++.h>
#include <random>
#include <fstream>
#include <cmath>
using namespace std;

//modifiable paramenters to handle the code
#define V 256 //total size of the universe
#define v2 16 //square root of the universe size
#define INF 99999 //infinity value for the graph
#define k 4 //range for the nearby drivers
#define graph_file "graph.txt" //path of the file containing the graph
#define driver_file "drivers.txt" //path of the file containing the driver locations

//class to manage the vEB tree operations
class vEBTree {
public:
    int u,u2; // Universe size and its square root
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
            summary = new vEBTree(upperSqrt); //the summary of the tree will also be a vEB tree
            cluster = new vEBTree*[upperSqrt]; //the cluster will be an array of vEB trees
            for (int i = 0; i < upperSqrt; ++i) {
                cluster[i] = new vEBTree(upperSquareRoot(size)); //each cluster will have a size of upper square root of the universe size
            }
        } 
        // If universe size is 2, no need for summary and cluster
        else {
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

    // Helper Functions for vEB Tree
    int high(int x) { return x / lowerSquareRoot(u); } 
    int low(int x) { return x % lowerSquareRoot(u); }
    int index(int high, int low) { return high * lowerSquareRoot(u) + low; }
    int lowerSquareRoot(int x) { return (int)pow(2, floor(log2(x) / 2)); }
    int upperSquareRoot(int x) { return (int)pow(2, ceil(log2(x) / 2)); }

    // INSERT FUNCTION
    void insert(int x){
        // If tree is empty, set min and max to x and we are done
        // this is done for lazy propagation
        if (min == -1) {
            min = max = x;
        }
        // If tree is not empty
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
            //set the max value
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


// Function to generate the graph from the txt file
void GenerateGraph(int graph[V][V]) {
    // Initialize with INF
    for (int i = 0; i < V; i++) {
        for (int j = 0; j < V; j++) {
            graph[i][j] = INF;
        }
    }
    
    // Set diagonal to 0
    for (int i = 0; i < V; i++) {
        graph[i][i] = 0;
    }

    // taking graph input from the file
    fstream file2(graph_file);
    if (!file2.is_open()) {
        cout << "Error opening graph file" << endl;
        return;
    }

    // Read the graph from the file and set the edge weights accordingly
    int u, v, w;
    while (file2 >> u >> v >> w) {
        graph[u][v] = w;
        graph[v][u] = w;
    }
}

// Function to run Dijkstra's algorithm on the graph
pair<vector<int>, vector<int>> RunDijkstra(int src, int graph[V][V]) {
    // Priority queue to select the minimum distance vertex. 
    // Stores pairs of (distance, vertex) for ordering by distance.
    priority_queue<pair<int, int>, vector<pair<int, int>>, greater<pair<int, int>>> pq;
    
    // Distance vector initialized to INF for all vertices
    vector<int> dist(V, INF);
    
    // Parent vector to keep track of the path; initialized to -1 (no parent)
    vector<int> g_parent(V, -1);
    
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
            // Check for a valid edge from u to v and update dist[v] if a shorter path is found
            if (graph[u][v] != INF && dist[u] != INF && dist[u] + graph[u][v] < dist[v]) {
                dist[v] = dist[u] + graph[u][v];  // Update the distance to vertex v
                g_parent[v] = u;  // Update the parent of v to u for path tracking
                pq.push(make_pair(dist[v], v));  // Push the updated distance and vertex to the queue
            }
        }
    }

    // Prepare the answer as a pair of vectors: distances and parent path
    pair<vector<int>, vector<int>> ans;
    ans.first = dist;
    ans.second = g_parent;
    
    return ans;
}

// Class to create a dialogue box for entering pickup and drop coordinates
class CoordinateDialog : public wxDialog {
private:
    // Text controls for coordinate input
    wxTextCtrl *sourceX, *sourceY, *destX, *destY;
    
    // Variables to store coordinates internally
    int srcX, srcY, dstX, dstY;

public:
    // Constructor to initialize the dialogue box UI elements
    CoordinateDialog(wxWindow* parent)
        : wxDialog(parent, wxID_ANY, "Taxi Management System - CS201 Project", 
                  wxDefaultPosition, wxSize(460, 440)) {

        // Set background and text colors for the dialogue
        this->SetBackgroundColour(wxColor(37,38,40));
        this->SetForegroundColour(*wxWHITE);
        
        wxPanel* panel = new wxPanel(this);
        wxBoxSizer* vbox = new wxBoxSizer(wxVERTICAL);

        // Title text with custom font settings
        wxStaticText* title = new wxStaticText(panel, wxID_ANY, "Taxi Management System");
        wxFont titleFont = title->GetFont();
        titleFont.SetPointSize(titleFont.GetPointSize() + 5);  // Increase font size
        titleFont.SetWeight(wxFONTWEIGHT_BOLD);  // Set to bold
        title->SetFont(titleFont);
        title->SetForegroundColour(*wxWHITE);
        vbox->Add(title, 0, wxALL | wxALIGN_CENTER, 10);

        // Description text with custom font settings
        wxTextCtrl* desc = new wxTextCtrl(panel, wxID_ANY, 
            "Enter the source and destination coordinates to get your Nearest Drivers, Shortest Route, and Fare for the taxi:",
            wxDefaultPosition, wxSize(460, 60),  wxTE_MULTILINE | wxTE_READONLY | wxBORDER_NONE);
        wxFont descFont = desc->GetFont();
        descFont.SetPointSize(descFont.GetPointSize() + 2);  // Adjust font size
        desc->SetFont(descFont);
        desc->SetBackgroundColour(wxColor(28,29,31));  // Set description background color
        desc->SetForegroundColour(*wxWHITE);
        vbox->Add(desc, 0, wxALIGN_CENTER | wxALL, 5);

        // Validators to ensure only integer inputs and set input ranges (0-15)
        wxIntegerValidator<int> validatorSrcX(&srcX);
        wxIntegerValidator<int> validatorSrcY(&srcY);
        wxIntegerValidator<int> validatorDstX(&dstX);
        wxIntegerValidator<int> validatorDstY(&dstY);
        
        validatorSrcX.SetRange(0, 15);
        validatorSrcY.SetRange(0, 15);
        validatorDstX.SetRange(0, 15);
        validatorDstY.SetRange(0, 15);
        
        // Input fields for source coordinates
        wxStaticText* srcLabel = new wxStaticText(panel, wxID_ANY, "Start Coordinates (0-15):");
        srcLabel->SetForegroundColour(*wxWHITE);
        wxBoxSizer* srcBox = new wxBoxSizer(wxHORIZONTAL);
        sourceX = new wxTextCtrl(panel, wxID_ANY, "", wxDefaultPosition, wxDefaultSize, 0, validatorSrcX);
        sourceY = new wxTextCtrl(panel, wxID_ANY, "", wxDefaultPosition, wxDefaultSize, 0, validatorSrcY);
        sourceX->SetForegroundColour(wxColor(120,120,120));
        sourceY->SetForegroundColour(wxColor(120,120,120));
        
        // Source coordinate labels
        wxStaticText* srcXLabel = new wxStaticText(panel, wxID_ANY, "X: ");
        srcXLabel->SetForegroundColour(*wxWHITE);
        srcBox->Add(srcXLabel, 0, wxALL | wxALIGN_CENTER_VERTICAL, 5);
        srcBox->Add(sourceX, 1, wxALL, 5);
        wxStaticText* srcYLabel = new wxStaticText(panel, wxID_ANY, "Y: ");
        srcYLabel->SetForegroundColour(*wxWHITE);
        srcBox->Add(srcYLabel, 0, wxALL | wxALIGN_CENTER_VERTICAL, 5);
        srcBox->Add(sourceY, 1, wxALL, 5);

        // Input fields for destination coordinates
        wxStaticText* dstLabel = new wxStaticText(panel, wxID_ANY, "Destination Coordinates (0-15):");
        dstLabel->SetForegroundColour(*wxWHITE);
        wxBoxSizer* dstBox = new wxBoxSizer(wxHORIZONTAL);
        destX = new wxTextCtrl(panel, wxID_ANY, "", wxDefaultPosition, wxDefaultSize, 0, validatorDstX);
        destY = new wxTextCtrl(panel, wxID_ANY, "", wxDefaultPosition, wxDefaultSize, 0, validatorDstY);
        destX->SetForegroundColour(wxColor(120,120,120));
        destY->SetForegroundColour(wxColor(120,120,120));
        
        // Destination coordinate labels
        wxStaticText* dstXLabel = new wxStaticText(panel, wxID_ANY, "X: ");
        dstXLabel->SetForegroundColour(*wxWHITE);
        dstBox->Add(dstXLabel, 0, wxALL | wxALIGN_CENTER_VERTICAL, 5);
        dstBox->Add(destX, 1, wxALL, 5);
        wxStaticText* dstYLabel = new wxStaticText(panel, wxID_ANY, "Y: ");
        dstYLabel->SetForegroundColour(*wxWHITE);
        dstBox->Add(dstYLabel, 0, wxALL | wxALIGN_CENTER_VERTICAL, 5);
        dstBox->Add(destY, 1, wxALL, 5);

        // Adding source and destination coordinate fields to the main sizer
        vbox->Add(srcLabel, 0, wxALL, 10);
        vbox->Add(srcBox, 0, wxEXPAND | wxLEFT | wxRIGHT, 10);
        vbox->Add(dstLabel, 0, wxALL, 10);
        vbox->Add(dstBox, 0, wxEXPAND | wxLEFT | wxRIGHT, 10);
        
        // OK and Cancel buttons
        wxButton* okButton = new wxButton(panel, wxID_OK, "&Find Drivers");
        okButton->SetForegroundColour(wxColor(120,120,120));
        wxButton* cancelButton = new wxButton(panel, wxID_CANCEL, "&Exit");
        cancelButton->SetForegroundColour(wxColor(120,120,120));
        wxBoxSizer* buttonBox = new wxBoxSizer(wxHORIZONTAL);
        buttonBox->Add(okButton, 0, wxALL, 5);
        buttonBox->Add(cancelButton, 0, wxALL, 5);
        vbox->Add(buttonBox, 0, wxALIGN_CENTER | wxALL, 10);

        // Credits section for project contributors
        wxTextCtrl* credits = new wxTextCtrl(panel, wxID_ANY, 
            "_________________________________________________________________________\nDSA (CS201) Project By: \n1. Gitansh Bansal (2023MCB1294)\n2. Mohakjot Dhiman (2023MCB1302)\n3. Kanav Puri (2023MCB1298)",
            wxDefaultPosition, wxSize(440, 100),  wxTE_MULTILINE | wxTE_READONLY | wxBORDER_NONE);
        wxFont credFont = credits->GetFont();
        credFont.SetPointSize(credFont.GetPointSize() - 3);  // Adjust font size for credits
        credits->SetFont(credFont);
        credits->SetBackgroundColour(wxColor(28,29,31));
        credits->SetForegroundColour(*wxWHITE);
        vbox->Add(credits, 0, wxALIGN_CENTER | wxALL, 10);
        
        // Set sizer for the panel and center the dialog on screen
        panel->SetSizer(vbox);
        Centre();
    }
    
    // Retrieve and calculate source vertex based on X and Y coordinates
    int GetSourceVertex() {
        int x, y;
        sourceX->GetValue().ToInt(&x);
        sourceY->GetValue().ToInt(&y);
        return x * v2 + y;  // Combine x and y to calculate source vertex
    }
    
    // Retrieve and calculate destination vertex based on X and Y coordinates
    int GetDestVertex() {
        int x, y;
        destX->GetValue().ToInt(&x);
        destY->GetValue().ToInt(&y);
        return x * v2 + y;  // Combine x and y to calculate destination vertex
    }
};

// Class to create a panel for displaying driver and ride information
class InfoPanel : public wxPanel {
private:
    vector<pair<int, int>> nearbyDrivers;    // Stores nearby driver locations
    wxTextCtrl* infoText;                    // Text control for displaying fare and ETA
    wxScrolledWindow* scrollWindow;          // Scrolled window to hold driver info panels
    wxBoxSizer* driversSizer;                // Sizer for arranging driver panels within the scrolled window
    int fare;                                // Variable for fare calculation
    int eta;                                 // Variable for estimated time of arrival (ETA)

public:
    // Constructor for initializing the panel layout
    InfoPanel(wxWindow* parent) : wxPanel(parent, wxID_ANY) {
        wxBoxSizer* mainSizer = new wxBoxSizer(wxVERTICAL);

        // Set panel background and text colors
        this->SetBackgroundColour(wxColor(37,38,40));
        this->SetForegroundColour(*wxWHITE);

        // Title for the panel
        wxStaticText* title = new wxStaticText(this, wxID_ANY, "Nearby Drivers");
        wxFont titleFont = title->GetFont();
        titleFont.SetPointSize(titleFont.GetPointSize() + 2);
        titleFont.SetWeight(wxFONTWEIGHT_BOLD);
        title->SetFont(titleFont);
        title->SetForegroundColour(*wxWHITE);
        mainSizer->Add(title, 0, wxALL | wxALIGN_CENTER, 10);

        // Scrolled window to display the list of nearby drivers
        scrollWindow = new wxScrolledWindow(this, wxID_ANY);
        scrollWindow->SetScrollRate(0, 5); // Sets vertical scroll speed
        driversSizer = new wxBoxSizer(wxVERTICAL); // Layout for driver panels
        scrollWindow->SetSizer(driversSizer);
        mainSizer->Add(scrollWindow, 1, wxEXPAND | wxALL, 5);

        // Text control for displaying additional information (fare, ETA)
        infoText = new wxTextCtrl(this, wxID_ANY, "", wxDefaultPosition, wxDefaultSize, wxTE_MULTILINE | wxTE_READONLY | wxTE_NO_VSCROLL);
        infoText->SetForegroundColour(*wxWHITE);
        infoText->SetBackgroundColour(wxColor(100,100,100));
        mainSizer->Add(infoText, 0, wxEXPAND | wxALL, 5);

        // Set main sizer for the panel
        SetSizer(mainSizer);
    }

    // Function to set the list of nearby drivers and update the display
    void SetNearbyDrivers(const vector<pair<int, int>>& drivers) {
        nearbyDrivers = drivers;
        UpdateDisplay(); // Refresh display to show the updated drivers
    }

    // Function to calculate and set fare based on ETA and number of drivers, then update display
    void SetFare(int x, int num) {
        // Base fare + (3 * ETA) - (driver-based discount)
        if (num==0){
            fare = 0;
            eta = 0;
        }
        else{
           fare = 50 + 3 * x - num*x*0.035;
            eta = x;
        }
        UpdateDisplay(); // Refresh display to show fare and ETA
    }

private:
    // Helper function to update the display with nearby drivers, fare, and ETA
    void UpdateDisplay() {
        // Clear existing content in the drivers sizer
        driversSizer->Clear(true);
        
        // If there are no drivers, show a message
        if (nearbyDrivers.empty()) {
            wxStaticText* noDrivers = new wxStaticText(scrollWindow, wxID_ANY, 
                "No drivers currently available in your area. Please try again later.");
            driversSizer->Add(noDrivers, 0, wxALL | wxALIGN_CENTER, 10);
        } else {
            // Add each driver as a panel within the scrolled window
            for (size_t i = 0; i < nearbyDrivers.size(); i++) {
                wxPanel* driverPanel = CreateDriverPanel(i + 1, nearbyDrivers[i]);
                driversSizer->Add(driverPanel, 0, wxEXPAND | wxALL, 5);
            }
        }

        // Set the summary text for fare and ETA
        if (fare==0){
            infoText->SetValue("");
        }
        else{
            wxString summary = wxString::Format("Fare Calculated: Rs. %d\nETA: %d minutes", fare, eta); 
            infoText->SetValue(summary);
        }
        // Set font style for info text (fare and ETA)
        wxFont largefont = infoText->GetFont();
        largefont.SetPointSize(largefont.GetPointSize() + 2); // Adjust size as needed
        largefont.MakeBold();
        infoText->SetFont(largefont);
        infoText->SetForegroundColour(*wxWHITE);

        // Refresh layout to show updates
        scrollWindow->FitInside(); // Adjust scroll area size
        scrollWindow->Layout();
        Layout();
        Refresh();
    }

    // Helper function to create a panel for each driver with index and location
    wxPanel* CreateDriverPanel(int index, const pair<int, int>& location) {
        wxPanel* panel = new wxPanel(scrollWindow, wxID_ANY);
        panel->SetBackgroundColour(wxColor(100, 100, 100)); // Light grey background
        
        wxBoxSizer* panelSizer = new wxBoxSizer(wxVERTICAL);
        
        // Driver header showing position number
        wxString headerText = wxString::Format("Driver #%d", index);
        wxStaticText* header = new wxStaticText(panel, wxID_ANY, headerText);
        wxFont headerFont = header->GetFont();
        headerFont.SetWeight(wxFONTWEIGHT_BOLD);
        header->SetFont(headerFont);
        
        // Location information of the driver
        wxString locationText = wxString::Format("Location: (%d, %d)", location.first, location.second);
        wxStaticText* locationInfo = new wxStaticText(panel, wxID_ANY, locationText);
        
        // Add header and location info to panel
        panelSizer->Add(header, 0, wxALL, 3);
        panelSizer->Add(locationInfo, 0, wxALL, 3);
        
        // Set sizer for the driver panel and return it
        panel->SetSizer(panelSizer);
        
        return panel;
    }
};


// Class to create a graphical canvas for visualizing a graph, such as a map of nodes and edges
class GraphCanvas : public wxPanel {
private:
    int graph[V][V];                        // Adjacency matrix to store graph weights
    vector<int> distances;                   // Stores distances from source to each node
    vector<int> parents;                     // Stores parent nodes for each node in shortest path tree
    int nodeRadius;                          // Radius for drawing each node
    int gridSize;                            // Size of the grid for layout
    int margin;                              // Margin from edges of the panel
    int sourceVertex;                        // Starting vertex for shortest path
    int destVertex;                          // Destination vertex for shortest path
    vector<int> nearbyDrivers;               // List of nearby driver nodes

    // Helper function to get the shortest path from source to destination
    vector<int> GetShortestPath() {
        vector<int> path;
        if (distances[destVertex] == INF) return path;  // No path if destination is unreachable
        
        int current = destVertex;
        while (current != -1) {
            path.push_back(current);
            current = parents[current];
        }
        reverse(path.begin(), path.end());
        return path;
    }
    
public:
    // Constructor to initialize GraphCanvas with parent window, source, destination, and parent nodes
    GraphCanvas(wxFrame* parent, int src, int dst, vector<int> g_parents) 
        : wxPanel(parent), nodeRadius(11), gridSize(v2), margin(40), 
          sourceVertex(src), destVertex(dst), parents(g_parents) {
        // Initialize graph matrix with INF (unreachable)
        for (int i = 0; i < V; i++) {
            for (int j = 0; j < V; j++) {
                graph[i][j] = INF;
            }
        }
        distances.resize(V, INF);           // Set initial distances to INF
        Bind(wxEVT_PAINT, &GraphCanvas::OnPaint, this); // Bind paint event
        SetBackgroundColour(*wxWHITE);      // Set background color
    }

    // Function to set the graph adjacency matrix and refresh display
    void SetGraph(int g[V][V]) {
        memcpy(graph, g, sizeof(int) * V * V);
        Refresh();
    }
    
    // Function to set distances from source to each node and refresh display
    void SetDistances(const vector<int>& dist) {
        distances = dist;
        Refresh();
    }

    // Function to set nearby drivers and map their positions onto the grid
    void SetNearbyDrivers(const vector<pair<int, int>>& drivers) {
        for (const auto& driver : drivers) {
            nearbyDrivers.push_back(driver.first * v2 + driver.second);
        }
        Refresh();
    }
    
private:
    // Paint handler for drawing nodes, edges, paths, and legends
    void OnPaint(wxPaintEvent& evt) {
        wxPaintDC dc(this);
        wxGraphicsContext* gc = wxGraphicsContext::Create(dc);
        if (!gc) return;

        wxSize size = GetSize();
        double cellWidth = (size.GetWidth() - 2 * margin) / gridSize;
        double cellHeight = (size.GetHeight() - 2 * margin) / gridSize;
        
        // Draw regular edges between nodes
        gc->SetPen(wxPen(wxColor(200, 200, 200), 10));  // Light gray for regular edges
        for (int i = 0; i < V; i++) {
            for (int j = i + 1; j < V; j++) {
                if (graph[i][j] != INF) {
                    // Calculate positions of the nodes
                    int x1 = margin + (i % gridSize) * cellWidth;
                    int y1 = margin + (i / gridSize) * cellHeight;
                    int x2 = margin + (j % gridSize) * cellWidth;
                    int y2 = margin + (j / gridSize) * cellHeight;
                    
                    // Draw line and weight for each edge
                    gc->StrokeLine(x1, y1, x2, y2);
                    wxString weight = wxString::Format("%d", graph[i][j]);
                    double tx = (x1 + x2) / 2;
                    double ty = (y1 + y2) / 2;
                    gc->SetFont(wxFont(8, wxFONTFAMILY_DEFAULT, wxFONTSTYLE_NORMAL, wxFONTWEIGHT_NORMAL), *wxBLACK);
                    wxDouble w, h;
                    gc->GetTextExtent(weight, &w, &h);
                    gc->DrawText(weight, tx - w/2, ty - h/2);
                }
            }
        }
        
        // Draw edges along the shortest path in light blue
        vector<int> shortestPath = GetShortestPath();
        if (shortestPath.size() > 1) {
            gc->SetPen(wxPen(wxColor(120, 170, 255), 12)); // Light blue for shortest path
            for (size_t i = 0; i < shortestPath.size() - 1; i++) {
                int v1 = shortestPath[i];
                int v22 = shortestPath[i + 1];
                int x1 = margin + (v1 % gridSize) * cellWidth;
                int y1 = margin + (v1 / gridSize) * cellHeight;
                int x2 = margin + (v22 % gridSize) * cellWidth;
                int y2 = margin + (v22 / gridSize) * cellHeight;
                
                // Draw line and weight for each edge on shortest path
                gc->StrokeLine(x1, y1, x2, y2);
                wxString weight = wxString::Format("%d", graph[v1][v22]);
                double tx = (x1 + x2) / 2;
                double ty = (y1 + y2) / 2;
                gc->SetFont(wxFont(8, wxFONTFAMILY_DEFAULT, wxFONTSTYLE_NORMAL, wxFONTWEIGHT_NORMAL), *wxBLACK);
                wxDouble w, h;
                gc->GetTextExtent(weight, &w, &h);
                gc->DrawText(weight, tx - w/2, ty - h/2);
            }
        }
        
        // Draw each node, with color depending on role (source, destination, nearby drivers, or default)
        for (int i = 0; i < V; i++) {
            int x = margin + (i % gridSize) * cellWidth;
            int y = margin + (i / gridSize) * cellHeight;
            
            wxColor color;
            if (i == sourceVertex) {
                color = wxColor(70, 255, 70);  // Green for source
            } else if (i == destVertex) {
                color = wxColor(255, 170, 60);  // Orange for destination
            } else if (find(nearbyDrivers.begin(), nearbyDrivers.end(), i) != nearbyDrivers.end()) {
                color = wxColor(255, 255, 0);  // Yellow for nearby drivers
            } else {
                color = wxColor(220, 220, 220);  // Gray for other nodes
            }
            
            // Draw node as a circle with label
            gc->SetBrush(wxBrush(color));
            gc->SetPen(wxPen(*wxBLACK, 1));
            gc->DrawEllipse(x - nodeRadius, y - nodeRadius, 2 * nodeRadius, 2 * nodeRadius);
            
            wxString number = wxString::Format("(%d,%d)", i / v2, i % v2);
            gc->SetFont(wxFont(8, wxFONTFAMILY_DEFAULT, wxFONTSTYLE_NORMAL, wxFONTWEIGHT_NORMAL), *wxBLACK);
            wxDouble w, h;
            gc->GetTextExtent(number, &w, &h);
            gc->DrawText(number, x - w/2, y - h/2);
        }

        // Draw legend at the bottom of the window
        int legendX = margin;
        int legendY = size.GetHeight() - 40; // Position near the bottom

        gc->SetFont(wxFont(10, wxFONTFAMILY_DEFAULT, wxFONTSTYLE_NORMAL, wxFONTWEIGHT_NORMAL), *wxBLACK);
        
        // Source label
        gc->SetBrush(wxBrush(wxColor(70, 255, 70))); // Green for source
        gc->DrawRectangle(legendX, legendY, 20, 12);
        gc->DrawText("Pick Up", legendX + 25, legendY);
        
        // Destination label
        legendX += 100;
        gc->SetBrush(wxBrush(wxColor(255, 170, 60))); // Orange for destination
        gc->DrawRectangle(legendX, legendY, 20, 12);
        gc->DrawText("Destination", legendX + 25, legendY);
        
        // Nearby drivers label
        legendX += 130;
        gc->SetBrush(wxBrush(wxColor(255, 255, 0))); // Yellow for nearby drivers
        gc->DrawRectangle(legendX, legendY, 20, 12);
        gc->DrawText("Taxis", legendX + 25, legendY);

        delete gc;  // Clean up graphics context
    }
};

// Main frame class for the application
class MainFrame : public wxFrame {
private:
    GraphCanvas* canvas;                   // Canvas to draw the graph and path
    int graph[V][V];                       // Graph adjacency matrix
    vEBTree* tree;                         // vEB tree for driver location storage
    InfoPanel* infoPanel;                  // Panel to display driver and fare information
    vector<int> currentDistances;          // Store distances for sorting
    vector<pair<int, int>> nearbyDrivers;  // Stores the nearby drivers' coordinates

public:
    // Constructor initializing source and destination vertices
    MainFrame(int sourceVertex, int destVertex) 
        : wxFrame(nullptr, wxID_ANY, "Taxi Route Visualizer", wxPoint(30, 30), wxSize(800, 800)) {
        
        // Generate graph and compute shortest paths using Dijkstra's algorithm
        GenerateGraph(graph);
        pair<vector<int>, vector<int>> dijkstra_output = RunDijkstra(sourceVertex, graph);
        currentDistances = dijkstra_output.first;  // Store distances
        vector<int> g_parents = dijkstra_output.second;

        // Initialize and populate the vEBTree with driver locations
        tree = new vEBTree(V);
        fstream file1(driver_file);
        if (!file1.is_open()) {
            cout << "Error opening drivers file" << endl;
            return;
        }
        int x, y;
        while (file1 >> x >> y) {
            tree->insertpair(x, y);
        }

        // Setup UI layout with main sizer
        wxBoxSizer* mainSizer = new wxBoxSizer(wxHORIZONTAL);
        SetSizer(mainSizer);

        // Canvas and info panel sizers for layout
        wxBoxSizer* canvasSizer = new wxBoxSizer(wxVERTICAL);
        wxBoxSizer* rightPane = new wxBoxSizer(wxVERTICAL);

        // Initialize canvas and info panel
        canvas = new GraphCanvas(this, sourceVertex, destVertex, g_parents);
        infoPanel = new InfoPanel(this);

        // Add canvas and info panel to respective sizers
        canvasSizer->Add(canvas, 1, wxEXPAND | wxALL, 5);
        rightPane->Add(infoPanel, 1, wxEXPAND | wxALL, 5);

        // Add sizers to the main layout
        mainSizer->Add(canvasSizer, 4, wxEXPAND | wxALL, 5);
        mainSizer->Add(rightPane, 1, wxEXPAND | wxALL, 5);

        SetSizer(mainSizer);
        this->Maximize(true);  // Maximize window on launch

        // Highlight nearby drivers based on source location
        HighlightNearbyDrivers(sourceVertex / v2, sourceVertex % v2);

        // Set graph data and distances on canvas and fare on info panel
        canvas->SetGraph(graph);
        canvas->SetDistances(currentDistances);
        //std::cout << "Number of nearby drivers: " << nearbyDrivers.size() << std::endl;
        infoPanel->SetFare(currentDistances[destVertex], nearbyDrivers.size());

    }    

    // Destructor to clean up resources
    ~MainFrame() {
        delete tree;
    }

    // Function to find non-empty points within a given radius k
    vector<pair<int, int>> findNonEmptyPoints(int x, int y, vEBTree* tree1) {
        vector<tuple<int, int, int>> points_with_distances;  // Stores points with distances
        int u = tree1->u2;

        int rowLow = max(0, x - k);
        int rowHigh = min(u - 1, x + k);

        for (int row = rowLow; row <= rowHigh; ++row) {
            int colLow = max(0, y - k);
            int colHigh = min(u - 1, y + k);

            int startKey = row * u + colLow;
            int col = tree1->successor(startKey - 1);

            // Check each column in the given range
            while (col != -1 && col / u == row && col % u <= colHigh) {
                if (col / u != x || col % u != y) {
                    int vertex = col;
                    int distance = currentDistances[vertex];
                    if (distance != INF) {  // Only include reachable drivers
                        points_with_distances.push_back({distance, col / u, col % u});
                    }
                }
                col = tree1->successor(col);
            }
        }

        // Sort points by distance
        sort(points_with_distances.begin(), points_with_distances.end());

        vector<pair<int, int>> sorted_points;
        for (const auto& point : points_with_distances) {
            sorted_points.push_back({get<1>(point), get<2>(point)});
        }

        return sorted_points;
    }

    // Function to highlight nearby drivers around the source location
    void HighlightNearbyDrivers(int srcX, int srcY) {
        nearbyDrivers = findNonEmptyPoints(srcX, srcY, tree);  // Find drivers within radius 4
        canvas->SetNearbyDrivers(nearbyDrivers);
        infoPanel->SetNearbyDrivers(nearbyDrivers);
    }
};

// Main application class to initialize the dialog and main frame
class App : public wxApp {
public:
    bool OnInit() {
        // Show the coordinate dialog to get source and destination vertices
        CoordinateDialog* dialog = new CoordinateDialog(nullptr);
        
        // If the user clicks OK, create the main frame with the selected source and destination vertices
        if (dialog->ShowModal() == wxID_OK) {
            int sourceVertex = dialog->GetSourceVertex();
            int destVertex = dialog->GetDestVertex();
            dialog->Destroy();
            
            // Create the main frame with the selected source and destination vertices
            MainFrame* frame = new MainFrame(sourceVertex, destVertex);
            // Show the main frame
            frame->Show();
            return true;
        }
        
        // If the user clicks Cancel, destroy the dialog and exit the application
        dialog->Destroy();
        return false;
    }
};

// Entry point for the application
wxIMPLEMENT_APP(App);
