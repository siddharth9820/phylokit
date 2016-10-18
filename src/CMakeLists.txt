cmake_minimum_required (VERSION 2.6)
project (phylokit)


link_directories(/usr/local/lib )
file(GLOB SOURCES *.cpp util/*.cpp)
file(GLOB HEADERS *.hpp util/*.hpp)
add_definitions(-std=c++14 -g -O3 -Wall)

add_library(phylokit ${SOURCES})

if(UNIX)
  SET(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
  SET(BUILD_SHARED_LIBRARIES OFF)
  SET(CMAKE_EXE_LINKER_FLAGS "-static -static-libgcc")
endif()


install(FILES ${HEADERS} DESTINATION include/)
install(TARGETS phylokit DESTINATION lib/)

