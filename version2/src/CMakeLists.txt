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
add_executable("${PROJECT_NAME}" "main.cpp" "./kinetics_solver/cvode_odesolver.cpp"  
"./kinetics_solver/differential_functions.cpp" "./superdrop_solver/collisions4SDM.cpp"
"./superdrop_solver/common2allsuperdrops.cpp" "./superdrop_solver/superdrop.cpp"
"./superdrop_solver/condensationgrowth4SDM.cpp" "./superdrop_solver/readwritefuncs4SDM.cpp") 
install(TARGETS "${PROJECT_NAME}" DESTINATION bin)        # ${CMAKE_INSTALL_PREFIX}/bin
install(FILES "testing.cpp" DESTINATION src)              # ${CMAKE_INSTALL_PREFIX}/src 

# # specify the C++ standard (or use flag in line 8 above)
# set(CMAKE_CXX_STANDARD 11)
# target_compile_features(${PROJECT_NAME} INTERFACE cxx_std_11)
# set_target_properties(${PROJECT_NAME} PROPERTIES CMAKE_CXX_STANDARD_REQUIRED ON)
# set_target_properties(${PROJECT_NAME} PROPERTIES CMAKE_CXX_EXTENSIONS ON)

# CVODE libaries stuff
# specify libraries to link against
set(LIBDIR /usr/local/sundials-6/lib)
find_library(CVODE_LIB NAMES sundials_cvode HINTS "${LIBDIR}")
find_library(NVECS_LIB NAMES sundials_nvecserial HINTS "${LIBDIR}")
#find_library(MANYVECS_LIB NAMES sundials_nvecmanyvector  HINTS "${LIBDIR}")

set(SUNDIALS_LIBS ${CVODE_LIB} ${NVECS_LIB})
target_link_libraries(${PROJECT_NAME} PRIVATE ${SUNDIALS_LIBS})


# CVODE include libraries stuff
target_include_directories(${PROJECT_NAME} PRIVATE /usr/local/sundials-6/include)