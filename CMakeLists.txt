cmake_minimum_required(VERSION 3.14)
project(CHMM)
set (CMAKE_CXX_STANDARD 20) #Sets C++ standard to C++20 

FIND_PACKAGE( Boost 1.40 COMPONENTS program_options REQUIRED ) #Need Boost for r-tree
INCLUDE_DIRECTORIES(include) #Find the Headers

#Find ESCAPE
FIND_LIBRARY(LIB_ESCAPE libescape
	PATH path/to/escape
) 

file(GLOB SOURCES "src/*.cpp") #Add .cpp files

add_executable(findXVal.out ${SOURCES})

TARGET_LINK_LIBRARIES(findXVal.out
	LINK_PUBLIC ${Boost_LIBRARIES} 
) 

target_link_libraries(findXVal.out
    	PRIVATE ${LIB_ESCAPE}
)
