# run cmake with " cmake -S [source dir of targerts] -B [build dir for makefile] "
# or enter build dir and just do " cmake [source dir] "

# set cmake version
cmake_minimum_required(VERSION "3.23.1")

# ass compiler flags
add_compile_options(-DNDEBUG -Wall -g -std=c++11)

# set project name and source directory of .cpp files
project("odesolver")
include_directories(${PROJECT_SOURCE_DIR})

# add executables to project (only one main() function!)
add_executable("${PROJECT_NAME}" "mainRK4.cpp"
"../../version2/src/superdrop_solver/superdrop.cpp"
"../../version2/src/superdrop_solver/common2allsuperdrops.cpp" 
"../../version2/src/superdrop_solver/readwritefuncs4SDM.cpp") 
#"../../version2/src/superdrop_solver/collisions4SDM.cpp" 
#"../../version2/src/superdrop_solver/condensationgrowth4SDM.cpp" 

