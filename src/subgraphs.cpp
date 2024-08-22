//subgraphs.cpp

#include<subgraphs.h>

#include "Escape/GraphIO.h"
#include "Escape/EdgeHash.h"
#include "Escape/Digraph.h"
#include "Escape/Triadic.h"
#include "Escape/FourVertex.h"
#include "Escape/Conversion.h"
#include "Escape/GetAllCounts.h"

//------------------------
//-----Get Edges List-----
//------------------------

//Reads edges from file and converts to a nice format C++ can work with

//Assume that you can store the edge list in C++ - will be much slower if you can't
//Relabels the vertices to be 0..n-1 - this will need to be changed if we ever want to do analysis on individual nodes (or at least we will need to keep the mapping)
//Has edges in both directions (not directed)
//output[indices[i]] is the first occurrence of vertex i in the first position - helps speed up some other functions
std::vector< std::pair<int,int> > getEdgeList(std::ifstream &graphFile, int &n, std::vector<int> &indices) {
    graphFile.clear();
    graphFile.seekg(0);
    std::vector< std::pair<int, int> > output;
    std::set< std::pair<int, int> > outputSet;

    bool csv;
    int temp1;
    char temp2;
    graphFile >> temp1;
    graphFile >> temp2;

    if(temp2 == ' ') {
        csv = false;
    }
    else {
        csv = true;
    }

    graphFile.clear();
    graphFile.seekg(0);

    std::unordered_map<int, int> labelMap;
    std::unordered_set<int> oldVertexLabels;
    n = 0;

    if(csv) {
        int a,b;
        char c;

        while (graphFile >> a >> c >> b) {
            if(oldVertexLabels.count(a) == 0) {
                oldVertexLabels.insert(a);
                labelMap[a] = n;
                ++n;
            }
            if(oldVertexLabels.count(b) == 0) {
                oldVertexLabels.insert(b);
                labelMap[b] = n;
                ++n;
            }

            outputSet.insert(std::make_pair(labelMap[a],labelMap[b]));
            outputSet.insert(std::make_pair(labelMap[b],labelMap[a]));
        }
    }

    else {
        int a,b;
        while (graphFile >> a >> b) {
            if(oldVertexLabels.count(a) == 0) {
                oldVertexLabels.insert(a);
                labelMap[a] = n;
                ++n;
            }
            if(oldVertexLabels.count(b) == 0) {
                oldVertexLabels.insert(b);
                labelMap[b] = n;
                ++n;
            }

            outputSet.insert(std::make_pair(labelMap[a],labelMap[b]));
            outputSet.insert(std::make_pair(labelMap[b],labelMap[a]));
        }
    }

    int counter = 0;
    for(auto myPair : outputSet) {
        output.push_back(myPair);
        if(myPair.first == counter) {
            indices.push_back(output.size()-1);
            ++counter;
        }
    }
    return output;
}


//-----------------------------
//-----Get Subgraph Counts-----
//-----------------------------

//Gets subgraph counts necessary for localized point cloud

//edges is from getEdgeList
//n is number of vertices in the graph
//indices is from getEdgeList
//If the degree of a vertex is less than cutOff, we don't include it in the localized point cloud
//Vector format: K2 K3c P3c P3 K3 K4c ... K4 (see ESCAPE for exact ordering)
//output[i][j] = c, i is vtx #, j is subgraph number (0 \le j < 15 for up to 4 vtxs), c is subgraph density
std::vector< std::vector<double> > subgraphCounts(const std::vector<std::pair<int, int> > &edges, const int n, const std::vector<int> &indices, int cutOff) {
    std::vector< std::vector<double> > output;
    int outputIndex = -1;
    std::cout << n;

    for(int i = 0; i < n; ++i) {
        std::cout << "When getting subgraph counts, iteration " << i+1 << " out of " << n << "\n";
        int index = indices[i];
        int counter = 0;
        std::unordered_map<int, int> myMap;
        std::unordered_set<int> nbrhd;
        while(edges[index].first == i) {
            nbrhd.insert(edges[index].second);
            myMap[edges[index].second] = counter;
            ++counter;
            ++index;
        }
        
        if(nbrhd.size() >= cutOff) {
            //Write to file so ESCAPE can use
            std::ofstream tempFile;
            tempFile.open("temp.edges");
            std::vector< std::pair<int, int> > tempVec; 

            for(auto v : nbrhd) {
                index = indices[v];
                while(edges[index].first == v) {
                    if((edges[index].second > v) && (nbrhd.count(edges[index].second) == 1)) {
                        tempVec.push_back(std::make_pair(myMap[v],myMap[edges[index].second]));
                    }
                    ++index;
                }
            }
            int deg = nbrhd.size();

            tempFile << deg << " " << tempVec.size() << "\n"; //Need this first which is why we use tempVec
            for(auto myPair : tempVec) {
                tempFile << myPair.first << " " << myPair.second << "\n";
            }
            tempFile.close();


            //Actual use ESCAPE (directly taken ESCAPE code, with only minor changes)
            if(tempVec.size() != 0) {
                
                ++outputIndex;
                output.push_back({});
                output[outputIndex].push_back((tempVec.size()*2.)/(deg*(deg-1.)));

                Graph g;
                if (loadGraph("temp.edges", g, 1, IOFormat::escape)) {
                    exit(1);
                }
                
                CGraph cg = makeCSR(g);
                cg.sortById();

                CGraph cg_relabel = cg.renameByDegreeOrder();
                cg_relabel.sortById();
                CDAG dag = degreeOrdered(&cg_relabel);
                
                (dag.outlist).sortById();
                (dag.inlist).sortById();

                double nonInd_three[4], nonInd_four[11];
                double ind_three[4], ind_four[11];
                  
                getAllThree(&cg_relabel, &dag, nonInd_three);
                
                for (int j=0; j<4; j++) {
                    ind_three[j] = 0;
                    for (int k=0; k<4; k++) { 
                        ind_three[j] += A3[j][k]*nonInd_three[k];
                    }
                }    

                getAllFour(&cg_relabel, &dag, nonInd_four);
                for (int j=0; j<11; j++) {
                    ind_four[j] = 0;
                    for (int k=0; k<11; k++) {
                        ind_four[j] += A4[j][k]*nonInd_four[k];
                    }
                }
                
                for(int j = 0; j < 4; ++j) {
                    output[outputIndex].push_back(ind_three[j]*6./(deg*(deg-1.)*(deg-2.)));
                }
                for(int j = 0; j < 11; ++j) {
                    output[outputIndex].push_back(ind_four[j]*24./(deg*(deg-1.)*(deg-2.)*(deg-3.)));
                }
            }
        }
    }    
    return output;
}


//-------------------------------------
//-----Largest Connected Component-----
//-------------------------------------

//Returns largest component, used in pruning algorithm

//Could easily be modified to return all connected components
//Rather than taking in an adjacency matrix, this takes in a vector of pairs and a minimum distance, d, two vertices are connected if their distance is less than d
//Currently this is the slowest part of the algorithm
//I tried implement a boost rtree to make this run faster, but it made the code about 5x slower. Possibly more fine tuning of the tree could fix the issue, however.
std::vector<bool> largestConnectedComponent(const std::vector< std::pair<double, double> > &coordinates, double d) {
    int n = coordinates.size(); //Reduce # of calls to .size()
    double d2 = pow(d,2); //d squared, precompute as much as possible 

    std::vector<bool> S(n, false); //Vertices we haven't dealt with yet

    int largestComponentSize = 0;
    std::vector<bool> output;

    //DFS and then remove the set and repeat
    while(true) {
        int startVertex = -1;
        std::vector<bool> component(n, false); //Vertices in the current component
        std::stack<int> myStack;

        //Find a vertex we have dealt with yet to start
        for(int i = 0; i < n; ++i) {
            if(!(S[i])) {
                startVertex = i;
                S[i] = true;
                component[i] = true;
                myStack.push(startVertex);
                i = n;
            }
        }

        if(startVertex == -1) {
            break;
        }

        //DFS
        while(!myStack.empty()) {
            int v = myStack.top();
            myStack.pop();

            for(int u = 0; u < n; ++u) {
                if(!S[u]) {
                    double dis = pow(coordinates[u].first - coordinates[v].first,2) + pow(coordinates[u].second - coordinates[v].second,2);
                    if(dis < d2) {
                        myStack.push(u);
                        component[u] = true;
                        S[u] = true;
                    }
                }
            }
        }

        int componentSize = count(component.begin(), component.end(), true);

        if(componentSize > largestComponentSize) {
            largestComponentSize = componentSize;
            output = component;
        }
    }

    return output;
}


//---------------
//-----Prune-----
//---------------

//Prunes coordinates to remove outliers

//Assumes coordinates are 2D. Shouldn't be hard to extend though.
//coordinates is from subgraph counts
//d is the selected radius - this is used in the getArea function
//c is the percentage of points that aren't outliers
//epsilon is the tolerance on the output radius size
//Need to have a pairs rather than a matrix because the edge density values may not all be the same anymore. 
//Not the most efficient memory storage - improve here if we need to (could return vector of indices e.g.)
//Output[i][j] = (a,b) i is subgraph #, j is vertex index, a is edge density, b is subgraph density (different from output of subgraph counts b/c edge density might differ btwn pairs)
std::vector< std::vector< std::pair<double, double> > >  prune(const std::vector< std::vector<double> > &coordinates, std::vector<double> &d, const double c, const double eps) {
    std::vector< std::vector< std::pair<double,double > > > output;
    const int numSubgraphs = coordinates[0].size();
    const int n = coordinates.size();
    
    for(int subgraph = 1; subgraph < numSubgraphs; ++subgraph) {
        double lb = 0.; //I really don't see a situation where you wouldn't want this to be 0
        double ub = 0.1; //Careful with this, it may be too small if running on a small number of vtxs

        std::cout << "While pruning, iteration " << subgraph << " out of " << numSubgraphs-1 << "\n";
        output.push_back({});
        std::vector<bool> indices;

        std::vector< std::pair<double, double> > coordinatePairs;
        for(int i = 0; i < n; ++i) {
            coordinatePairs.push_back(std::make_pair(coordinates[i][0], coordinates[i][subgraph]));
        }

        while(ub - lb > eps) {
            std::cout << "lb = " << lb << " ub = " << ub << "\n";
            indices = largestConnectedComponent(coordinatePairs, (ub+lb)/2.);
            int indicesSize = count(indices.begin(), indices.end(), true);

            if(indicesSize < c*n) {
                lb = (ub+lb)/2.;
            }
            else {
                ub = (ub+lb)/2.;
            }
        }

        indices = largestConnectedComponent(coordinatePairs, ub);

        d.push_back(ub);
        for(int i = 0; i < n; ++i) {
            if(indices[i]) {
                output[subgraph-1].push_back(std::make_pair(coordinates[i][0], coordinates[i][subgraph]));
            }
        }
        std::cout << "Radius = " << ub << "\n"; 
    }
    return output;
}


//------------------
//-----Get Area-----
//------------------

//Gives the area of the balls of the pruned localized point cloud

//Just does Monte Carlo b/c it's fast and would extend to higher dimensions
//If this becomes a choke point, you could probably use a boost r tree to speed up the time, but determining the largest connected component seems much slower than this
//d is ub from previous fn
//minX and maxX are used in getVals function, they don't need to pass in any specific value
//numVertices is the number of vertices we ran the plainFlagAlgebra method on - in this case that is 6
//PFAMIndex is calculated automatically and is necessary for finding the correct PFAM file
//numIt is the number of iterations run by the Monte Carlo method
double area(const std::vector< std::pair<double, double> >  &coordinates, const double d, double &minX, double &maxX, const int numVertices, const int PFAMindex, const int numIt) {
    double a, b;
    int n = 0;
    maxX = 0.;
    minX = 1.;
    double minY = 1.;
    double maxY = 0.;

    //Determine # of data points and the range of the box we need to pull random points from 
    for(auto myPair: coordinates) {
        a = myPair.first;
        b = myPair.second;
        ++n;
        if(a > maxX) {
            maxX = a;
        }
        if(a < minX) {
            minX = a;
        }
        if(b > maxY) {
            maxY = b;
        }
        if(b < minY) {
            minY = b;
        }
    }

    maxX = std::min(maxX+d,1.);
    minX = std::max(minX-d,0.);
    maxY = maxY+d;
    minY = std::max(minY+d,0.);

    int numHits = 0;

    std::ifstream ub;
    std::ifstream lb;

    if(PFAMindex < 4) {
        ub.open("../PFAM/PFAM_"+std::to_string(numVertices)+"_3_"+std::to_string(PFAMindex)+".txt");
    }
    else {
        ub.open("../PFAM/PFAM_"+std::to_string(numVertices)+"_4_"+std::to_string(PFAMindex-4)+".txt");
    }

    if((PFAMindex == 0) || (PFAMindex == 3)) { //K3c or K3
        lb.open("../PFAM/PFAM_"+std::to_string(numVertices)+"_3_"+std::to_string(PFAMindex)+"_1.txt");
    }

    else if((PFAMindex == 4) || (PFAMindex == 14)) { //K4c or K4
        lb.open("../PFAM/PFAM_"+std::to_string(numVertices)+"_4_"+std::to_string(PFAMindex-4)+"_1.txt");
    }

    std::vector<double> ubVec;
    std::vector<double> lbVec;

    double val;
    while(ub >> val) {
        ubVec.push_back(val);
    }
    if(lb.is_open()) {
        while(lb >> val) {
            lbVec.push_back(val);
        }
    }
    int PFAMn = ubVec.size();

    //Set up generator 
    std::random_device rand_dev;
    std::minstd_rand generator(rand_dev());
    generator.seed(time(NULL)); //Can change
    std::uniform_real_distribution<double> distX(minX, maxX);
    std::uniform_real_distribution<double> distY(minY, maxY);
    
    for(int i = 0; i < numIt; ++i) {
        if((i % 100000) == 0) {
            std::cout << "Iteration " << i << "\n";
        }
        double x = distX(generator);
        double y = distY(generator);

        //Check if range from PFAM
        //Has the benefit of speeding up some calculations (specifically complete graphs)
        bool range = true;

        if(ubVec[floor(x*PFAMn)] < y) {
            range = false;
        } 

        if(lb.is_open()) {
            if(lbVec[floor(x*PFAMn)] > y) {
                range = false;
            }
        }

        if(range) {
            for(int j = 0; j < n; ++j) {
                a = coordinates[j].first;
                double dis1 = a - x;
                if (std::abs(dis1) < d) { //More function calls - but will reduce run time
                    b = coordinates[j].second;
                    double dis2 = b - y;
                    if(std::abs(dis2) < d) {
                        if( std::sqrt(dis1*dis1 + dis2*dis2) < d) {
                            ++numHits;
                            break;
                        }
                    }
                }
            }
        }
    }

    return ((double)numHits)/((double)numIt)*(maxX-minX)*(maxY-minY);
}

//---------------------------
//-----Get Relative Area-----
//---------------------------

//Returns x-value

//Output[i] = c, i = subgraph index, c = values
//NumVertices are from PFAM (6 for us)
//d is output by the prune function and is the radius of the balls such that ~95% of the points form a connected region
//n is size of PFAM files (100 for us)
//numIt is how accurate we want area calculations to be (used in the area function)
std::vector<double> getVals(const std::vector< std::vector< std::pair<double, double> > > &prunedCoordinates, const std::vector<double> & d, const int numVertices, const int n, const int numIt) {
    std::vector<double> output;

    for(int i = 0; i < prunedCoordinates.size(); ++i) {
        std::cout << "Calculating relative area, iteration " << i+1 << " out of " << prunedCoordinates.size() << "\n";
        std::ifstream lb;
        std::ifstream ub;
        double xMin;
        double xMax;

        double myArea = area(prunedCoordinates[i],d[i],xMin, xMax, numVertices, i, numIt);

        if(i <= 3) { //3 vtxs
            ub.open("../PFAM/PFAM_"+std::to_string(numVertices)+"_3_"+std::to_string(i)+".txt");
        }
        else { //4 vtxs
            ub.open("../PFAM/PFAM_"+std::to_string(numVertices)+"_4_"+std::to_string(i-4)+".txt");
        }

        if((i == 0) || (i == 3)) { //K3c or K3
            lb.open("../PFAM/PFAM_"+std::to_string(numVertices)+"_3_"+std::to_string(i)+"_1.txt");
        }

        else if((i == 4) || (i == 14)) { //K4 or K4c
            lb.open("../PFAM/PFAM_"+std::to_string(numVertices)+"_4_"+std::to_string(i-4)+"_1.txt");
        }

        double extremalArea = 0.;

        int index = 0;
        double val;
        while(ub >> val) {
            if((index > xMin*n - 1) && (index < xMax*n)) {
                extremalArea += val/n;
            }
            ++index; 
        }

        if(lb.is_open()) {
            index = 0;
            while( lb >> val) {
                if((index > xMin*n-1) && (index < xMax*n)) {
                    extremalArea -= val/n;
                }
                ++index; 
            }
        }

        output.push_back(myArea/extremalArea);
    }

    return output;
}


//-----------------------
//-----Main Function-----
//-----------------------

//This is the main function to call when calculating the x-value
//The directory you call this from should have a directory called "Data" which contains the file
//The data file is a a list of vertices in an edge either separated by a comma or a space
//The vertices must be labeled with integers, but the labels may skip values or not start at 0
//The directory must also contain a folder called PFAM which includes the results of bounds from the Plain Flag Algebra Method
//These files are labeled as "PFAM_a_b_c" where a is the number of vertices we run the method on, b is the number of vertices in our subgraph, and c is the index of the subgraph - we use the same ordering as in ESCAPE for our subgraphs. 
//The files contain the objective values, separated by line breaks, of the program (max H s.t i/n <= K2 <= (i+1)/n) for all i \in {0,...,n-1}
//n is the number of rows in file, and the values are ordered by i
//The results are stored in a ./Results/fileName (this directory will be created if it does not exist) with subdirectories Coordinates and Figures
//Coordinates.txt is a list of the coordinate values to create the localized scatter plots - the first column is edge density and the others are other subgraph densities in the same ordering as the ESCAPE ordering.
//PrunedCoordinate files are the pruned coordinates where the indexing is based on the indexing in Coordinates.txt
//Results.txt gives all x-values
void mainFn(std::string fileName, std::string fileExtension) {
    std::ifstream graphFile;
    int n;
    std::vector<int> indices;

    //If necessary, create file structure
    if (!std::filesystem::is_directory("../Results/") || !std::filesystem::exists("../Results/")) { // Check if src folder exists
        std::filesystem::create_directory("../Results/"); 
    }   
    if (!std::filesystem::is_directory("../Results/"+fileName) || !std::filesystem::exists("../Results/"+fileName)) { // Check if src folder exists
        std::filesystem::create_directory("../Results/"+fileName); 
    }
    if (!std::filesystem::is_directory("../Results/"+fileName+"/Coordinates") || !std::filesystem::exists("../Results/"+fileName+"/Coordinates")) { // Check if src folder exists
        std::filesystem::create_directory("../Results/"+fileName+"/Coordinates"); 
    }  
    if (!std::filesystem::is_directory("../Results/"+fileName+"/Figures") || !std::filesystem::exists("../Results/"+fileName+"/Figures")) { // Check if src folder exists
        std::filesystem::create_directory("../Results/"+fileName+"/Figures"); 
    }  

    //Open file with fileName from data folder
    graphFile.open("../Data/"+fileName+fileExtension);
    std::vector< std::pair<int, int> > edges = getEdgeList(graphFile,n,indices); //Read in edge list into a nice format
    std::vector< std::vector<double> > counts =  subgraphCounts(edges, n, indices); //Get coordinates for localized point clouds
    
    //Print coordinates.txt
    std::ofstream temp;
    temp.open("../Results/"+fileName+"/Coordinates/coordinates.txt");
    for(auto myVec: counts) {
        for(auto val: myVec) {
            temp << val << " ";
        }
        temp << "\n";
    }

    std::vector<double> d;
    std::vector< std::vector< std::pair<double, double> > > prunedCoordinates = prune(counts, d); //Prune Coordinates

    //Print pruned coordinates
    for(int i = 0; i < prunedCoordinates.size(); ++i) {
        std::ofstream tempFile;
        tempFile.open("../Results/"+fileName+"/Coordinates/prunedCoordinates"+ std::to_string(i) + ".txt");
        for(int j = 0; j < prunedCoordinates[i].size(); ++j) {
            tempFile << prunedCoordinates[i][j].first << " " << prunedCoordinates[i][j].second << "\n";
        }
        tempFile.close();
    }
    
    //Get and print x-vals
    std::ofstream areas;
    areas.open("../Results/"+fileName+"/results.txt");
    std::vector<double> results = getVals(prunedCoordinates, d, 6, 100, 10E5);
    areas << "K3c: " << results[0] << "\n";
    areas << "P3c: " << results[1] << "\n";
    areas << "P3: " << results[2] << "\n";
    areas << "K3: " << results[3] << "\n";
    areas << "K4c: " << results[4] << "\n";
    areas << "e+v+v: " << results[5] << "\n";
    areas << "e+e: " << results[6] << "\n";
    areas << "P3+v: " << results[7] << "\n";
    areas << "K3+v: " << results[8] << "\n";
    areas << "K1_3: " << results[9] << "\n";
    areas << "P4: " << results[10] << "\n";
    areas << "K3+pendant: " << results[11] << "\n";
    areas << "C4: " << results[12] << "\n";
    areas << "K4-e: " << results[13] << "\n";
    areas << "K4: " << results[14] << "\n";
    double avg = 0.;
    for(double val: results) {
        avg += val;
    }
    areas << "Average: " << avg/15. << "\n";
    //ESCAPE uses a temporary file that we make sure to close
    areas.close();
    std::remove("temp.edges");

    std::cout << "\n\nx-values:\n---------\n\n";
    std::cout << "K3c: " << results[0] << "\n";
    std::cout  << "P3c: " << results[1] << "\n";
    std::cout  << "P3: " << results[2] << "\n";
    std::cout  << "K3: " << results[3] << "\n";
    std::cout  << "K4c: " << results[4] << "\n";
    std::cout  << "e+v+v: " << results[5] << "\n";
    std::cout  << "e+e: " << results[6] << "\n";
    std::cout  << "P3+v: " << results[7] << "\n";
    std::cout << "K3+v: " << results[8] << "\n";
    std::cout  << "K1_3: " << results[9] << "\n";
    std::cout  << "P4: " << results[10] << "\n";
    std::cout  << "K3+pendant: " << results[11] << "\n";
    std::cout  << "C4: " << results[12] << "\n";
    std::cout  << "K4-e: " << results[13] << "\n";
    std::cout  << "K4: " << results[14] << "\n";
    std::cout << "Average: " << avg/15. << "\n";

    return;
}