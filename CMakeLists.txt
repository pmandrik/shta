#----------------------------------------------------------------------------
# Setup the project
#
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
cmake_policy(SET CMP0011 NEW)

project(shta)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -ldl")

if(TURN_ON_WARNING)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wunused-variable -Wunused-parameter")
else(TURN_ON_WARNING)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unused-variable -Wno-unused-parameter")
endif(TURN_ON_WARNING)

#-- Extern

#-- ROOT
find_package(ROOT)

if(ROOT_FOUND)
  message(STATUS "found root v" ${ROOT_VERSION_MAJOR}.${ROOT_VERSION_MINOR}.${ROOT_VERSION_PATCH} " at: " ${ROOT_LIBRARY_DIR})
  add_definitions(-DUSE_ROOT=1)

  include(${ROOT_USE_FILE})
  include_directories(${CMAKE_SOURCE_DIR} ${ROOT_INCLUDE_DIRS})
  add_definitions(${ROOT_CXX_FLAGS})

  # ROOT_GENERATE_DICTIONARY(G__Event Event.h LINKDEF EventLinkDef.h)
  # add_library(Event SHARED Event.cxx G__Event.cxx)
  # target_link_libraries(Event ${ROOT_LIBRARIES})
  # target_link_libraries(Main Event)
else(ROOT_FOUND)
  message(STATUS "root not found")
endif(ROOT_FOUND)

#-- shta
set(SHTA_SOURCE_DIR "${CMAKE_SOURCE_DIR}")

include_directories( ${SHTA_SOURCE_DIR} )
include_directories( ${SHTA_SOURCE_DIR}/negative_scripts/ )
file(GLOB sources_shta  ${SHTA_SOURCE_DIR}/*.cpp)
file(GLOB headers_shta  ${SHTA_SOURCE_DIR}/*.hh)

#set(sources ${sources_toprex} ${sources})
#set(headers ${headers_toprex} ${headers})

#message(STATUS  ${sources_toprex} )
#message(STATUS  ${sources} )

#-- BIN
#--build-options
add_definitions(-DSHTA_SOURCE_DIR="${CMAKE_SOURCE_DIR}")
add_definitions(-DSHTA_BINARY_DIR="${CMAKE_BINARY_DIR}")

#-- BIN
if(USER_FUNCTION_NAME)
  message(STATUS "build shta user function " ${USER_FUNCTION_NAME} )
  get_filename_component(funcname ${USER_FUNCTION_NAME} NAME_WE)
  add_library( ${funcname} SHARED ${USER_FUNCTION_NAME} )
  target_link_libraries( ${funcname} ${ROOT_LIBRARIES} )
else(USER_FUNCTION_NAME)
  message(STATUS "build shta tests" )
  #file( GLOB SHTA_TESTS ${SHTA_SOURCE_DIR}/msu_13Tev/*.cpp ) # ${SHTA_SOURCE_DIR}/negative_scripts/*.cpp 
  file( GLOB SHTA_TESTS ${SHTA_SOURCE_DIR}/negative_scripts/*.cpp)
  foreach( test ${SHTA_TESTS} )
    message( STATUS ${test} )
    get_filename_component(testname ${test} NAME_WE)
    add_executable( ${testname} ${test} )
    target_link_libraries( ${testname} ${ROOT_LIBRARIES} )
  endforeach( test ${SHTA_TESTS} )
endif(USER_FUNCTION_NAME)

#install(TARGETS shta RUNTIME DESTINATION ${PROJECT_SOURCE_DIR})



#add_executable(TopReX TopRexMain.cpp ${sources} ${headers})
#target_link_libraries(TopReX ${ROOT_LIBRARIES})
#install(TARGETS TopReX RUNTIME DESTINATION ${PROJECT_SOURCE_DIR})





