cmake_minimum_required(VERSION 3.23)

project(
  2to2CollisionWidth
  VERSION 1.0.0
  LANGUAGES CXX
  DESCRIPTION
    "2 -> 2 Collision Widths"
)

#add_compile_options(-Ofast)

#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Ofast" )

include_directories(/opt/homebrew/include)
find_library(CUBA_LIB NAMES cuba PATHS /opt/homebrew/lib)
find_package(Python3 REQUIRED COMPONENTS Interpreter Development)
find_package(GSL REQUIRED)  



# Ask Python for NumPy include path dynamically
execute_process(
    COMMAND ${Python3_EXECUTABLE} -c "import numpy; print(numpy.get_include())"
    OUTPUT_VARIABLE NumPy_INCLUDE_DIR
    OUTPUT_STRIP_TRAILING_WHITESPACE
)

message(STATUS "Python3_EXECUTABLE = ${Python3_EXECUTABLE}")
message(STATUS "Python3_VERSION    = ${Python3_VERSION}")
message(STATUS "Python3_INCLUDE_DIRS = ${Python3_INCLUDE_DIRS}")
message(STATUS "Python3_LIBRARIES  = ${Python3_LIBRARIES}")

message(STATUS "Found NumPy include dir: ${NumPy_INCLUDE_DIR}")

add_executable(2to2collisionwidth 2to2.cpp 2to2collisionwidth.cpp 2to2collisionwidth.hpp constants.hpp)

target_compile_features(2to2collisionwidth PUBLIC cxx_std_17)
target_link_libraries(2to2collisionwidth GSL::gsl ${CUBA_LIB} ${Python3_LIBRARIES})
target_include_directories(2to2collisionwidth PRIVATE ${Python3_INCLUDE_DIRS} ${NumPy_INCLUDE_DIR})
set(CMAKE_BUILD_TYPE "Release")