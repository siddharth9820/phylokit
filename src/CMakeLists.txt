cmake_minimum_required (VERSION 2.6)
project (phylokit)


link_directories(/usr/local/lib )
file(GLOB SOURCES *.cpp util/*.cpp)
file(GLOB HEADERS *.hpp)
file(GLOB UHEADERS util/*.hpp)
add_definitions(-std=c++14 -g -O3 -Wall)


if(UNIX)
  SET(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
  SET(BUILD_SHARED_LIBRARIES OFF)
  SET(CMAKE_EXE_LINKER_FLAGS "-static -static-libgcc")
endif()


find_package(Boost)
if(Boost_FOUND)
#  INCLUDE_DIRECTORIES("${Boost_INCLUDE_DIRS}" /usr/include/python2.7)
  SET(Boost_USE_STATIC_LIBS OFF)
  SET(Boost_USE_MULTITHREADED ON)
  SET(Boost_USE_STATIC_RUNTIME OFF)
#  FIND_PACKAGE(Boost 1.45.0 COMPONENTS python)

#  ADD_LIBRARY(phylokit SHARED pyext/bindings.cxx ${SOURCES})
  ADD_LIBRARY(phylokit STATIC ${SOURCES})
  TARGET_LINK_LIBRARIES(phylokit ${Boost_LIBRARIES} )


  
elseif(NOT Boost_FOUND)

  add_library(phylokit STATIC ${SOURCES})

endif()

install(FILES ${HEADERS} DESTINATION include/)
install(FILES ${UHEADERS} DESTINATION include/util/)
install(TARGETS phylokit DESTINATION lib/)
