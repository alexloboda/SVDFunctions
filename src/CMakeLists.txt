cmake_minimum_required(VERSION 3.13)
project(src)

find_package (Threads REQUIRED)
find_package(Eigen3 3.3 REQUIRED NO_MODULE)

INCLUDE_DIRECTORIES( /usr/lib/R/library/BH/include )
INCLUDE_DIRECTORIES( /usr/lib/R/library/Rcpp/include/ /usr/include/R )
LINK_DIRECTORIES( /usr/lib64/R/lib )

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wno-long-long -pedantic -pthread")
set(CMAKE_CXX_STANDARD 17)

set(LIBRARIES_TO_LINK profiler "${CMAKE_THREAD_LIBS_INIT}" Eigen3::Eigen)

file(GLOB SVDF_SRC "*.cpp")

file(GLOB SVDF_TEST_SRC "test/*")

add_executable(test ${SVDF_TEST_SRC} ${SVDF_SRC} test/mvn_stat_test.cpp include/subsample.h subsample.cpp four_russians.cpp matching.cpp include/matching.h include/hw.h)

TARGET_LINK_LIBRARIES(test R z ${LIBRARIES_T_LINK})
