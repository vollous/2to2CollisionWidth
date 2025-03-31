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

find_package(GSL REQUIRED)  

add_executable(2to2collisionwidth 2to2.cpp 2to2collisionwidth.cpp 2to2collisionwidth.hpp constants.hpp)

target_compile_features(2to2collisionwidth PUBLIC cxx_std_17)
target_link_libraries(2to2collisionwidth GSL::gsl ${CUBA_LIB})

set(CMAKE_BUILD_TYPE "Release")