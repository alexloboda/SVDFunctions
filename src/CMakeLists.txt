cmake_minimum_required(VERSION 3.13)
project(src)

find_package (Threads REQUIRED)

INCLUDE_DIRECTORIES( /usr/lib/R/library/BH/include )
INCLUDE_DIRECTORIES(/usr/lib/R/library/RcppEigen/include)
INCLUDE_DIRECTORIES( /usr/lib/R/library/Rcpp/include/ /usr/include/R )
LINK_DIRECTORIES( /usr/lib64/R/lib )
LINK_DIRECTORIES(/usr/lib/R/library/RcppEigen/libs)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wno-long-long -pedantic -pthread")
set(CMAKE_CXX_STANDARD 17)

set(LIBRARIES_TO_LINK profiler "${CMAKE_THREAD_LIBS_INIT}")

file(GLOB SVDF_SRC "*.cpp")

file(GLOB SVDF_TEST_SRC "test/*")

add_executable(test_run ${SVDF_TEST_SRC} ${SVDF_SRC} hw.cpp)

TARGET_LINK_LIBRARIES(test_run R z ${LIBRARIES_T_LINK})
