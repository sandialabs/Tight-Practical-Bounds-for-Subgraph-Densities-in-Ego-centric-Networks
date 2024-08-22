Overview of Code
----------------

This code calculates the x-value of graphs for all subgraphs on 3 or 4 vertices, and generates figures of the localized point clouds. 


Installation
------------

(1) Download and make ESCAPE (https://bitbucket.org/seshadhri/escape/src/master/)
(2) If the ESCAPE library is not installed in a location the compiler will look for it, in CMakeLists.txt change FIND_LIBRARY(LIB_ESCAPE libescape) to FIND_LIBRARY(LIB_ESCAPE libescape PATH /path/to/escape)
(3) Navigate in terminal to build
(4) Run cmake ..
(5) Run make
(6) To calculate x-values run ./findXVal.out the code will prompt you for the file name, which must be located in Subgraphs_Code/Data. x-values will be stored in Subgraphs_Code/Results/file_name. If the file name is incorrect or if the file is incorrectly formatted, the code will seg fault. 
(7) To generate figures run python generateFigueres . It will prompt you for a file name. The figures will also be stored in Subgraphs_Code/Results/file_name 

Note: In the Data directory we include the file musae_facebook.csv (from Multi-scale Attributed Node Embedding by Rozemberczki, Allen, Sarkar). This can be used to test the build.

File Format
-----------

All files must be located in Subgraphs_Code/Data . The files are edge lists that are either space separated (and have a .txt extension) or are comma separated (and have a .csv extension). The files may not contain any preamble at the beginning, often with data from the SNAP database, one must delete the first few lines of the file. The labels of the vertices may be any integers (not causing overflow in C++), it does not even need to be a contiguous set of integers. Edges only need to appear once, e.g. the file does not need to include both i j and j i. However, repetitions of edges will not cause errors in the program. If the data comes from an directed graph, the code will just assume that the graph is undirected.

Plain Flag Algebra Method
-------------------------

Results from the plain flag algebra method have been included in the directory PFAM. These are the results of running the method on 6 vertices and this should be sufficient for most applications. If the user is interested in some plain flag algebra code, we direct them to https://github.com/clmatt/Plain-Flag-Algebra for example. 
