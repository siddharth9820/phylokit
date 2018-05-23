cmake_minimum_required (VERSION 2.6)
project (phylokit)


file(GLOB SOURCES *.cpp util/*.cpp)
file(GLOB HEADERS *.hpp)
file(GLOB UHEADERS util/*.hpp)

add_definitions(-std=c++14 -g -O3 -Wall)


find_package(Boost)
if(Boost_FOUND)
  INCLUDE_DIRECTORIES("${Boost_INCLUDE_DIRS}" /usr/include/python2.7)
  SET(Boost_USE_STATIC_LIBS OFF)
  SET(Boost_USE_MULTITHREADED ON)
  SET(Boost_USE_STATIC_RUNTIME OFF)
  FIND_PACKAGE(Boost 1.45.0 COMPONENTS python)

  ADD_LIBRARY(phylokit_python SHARED pyext/bindings.cxx ${SOURCES} )
  SET_TARGET_PROPERTIES(phylokit_python PROPERTIES OUTPUT_NAME pyphylokit)
  TARGET_LINK_LIBRARIES(phylokit_python ${Boost_LIBRARIES} )



elseif(NOT Boost_FOUND)
  INCLUDE_DIRECTORIES("boost/" /usr/include/python2.7)
endif()

ADD_LIBRARY(phylokit_shared SHARED ${SOURCES} )
SET_TARGET_PROPERTIES(phylokit_shared PROPERTIES OUTPUT_NAME phylokit)
TARGET_LINK_LIBRARIES(phylokit_shared ${Boost_LIBRARIES} )


ADD_LIBRARY(phylokit_static STATIC ${SOURCES})
SET_TARGET_PROPERTIES(phylokit_static PROPERTIES OUTPUT_NAME phylokit)
set_target_properties(phylokit_static PROPERTIES LINK_FLAGS "-static" )

install(FILES ${HEADERS} DESTINATION include/)
install(FILES ${UHEADERS} DESTINATION include/util/)
install(TARGETS phylokit_static DESTINATION lib/)
install(TARGETS phylokit_shared DESTINATION lib/)


ADD_EXECUTABLE(phylokit_test ${SOURCES} Tests/test.cpp Tests/runner.cpp)
