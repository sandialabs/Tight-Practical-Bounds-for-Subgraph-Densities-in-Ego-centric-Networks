#include <iostream>
#include <chrono> //Time the function
#include <string>

#include "subgraphs.h" 

int main() {
    //The file we wish to calculate the x-value for must be found in a subdirectory called data
    std::cout << "\nThis program takes a file of edges in the Data folder and finds the x-value.\nSee the readme for more information.\n\n";
    std::cout << "Enter the file name, with extension (the extension must either be .txt or .csv): ";
    std::string totalFileName;
    std::cin >> totalFileName; //Name of file

    //Extract actual name and extension
    std::string fileName = totalFileName.substr(0, totalFileName.size() - 4);
    std::string fileExtension = totalFileName.substr(totalFileName.size() - 4, totalFileName.size());

    //More information about this function is included in subgraphs.h
    //The function will output the x-values and some other information in ./Results/fileName (this folder will be created if it does not exist)
    //The createGraphs.py code will create the figures (see that file for more information)
    auto start = std::chrono::high_resolution_clock::now();
    mainFn(fileName, fileExtension);
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    std::cout << "Total running time: " << duration.count() << " microseconds.\n";

    return 0;
}