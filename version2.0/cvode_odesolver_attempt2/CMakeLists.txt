# run cmake with " cmake -S [source dir of targerts] -B [build dir for makefile] "
# or enter build dir and just do " cmake [source dir] "

# set cmake version
cmake_minimum_required(VERSION "3.23.1")

# ass compiler flags
add_compile_options(-DNDEBUG -Wall -g)

# set project name and source directory of .cpp files
project("Testing1")
include_directories(${PROJECT_SOURCE_DIR})

# add executables to project (only one main() function!)
add_executable("${PROJECT_NAME}" "testing.cpp" "cvode_differentials.cpp" "cvode_odesolver.cpp")

install(TARGETS "${PROJECT_NAME}" DESTINATION bin)        # ${CMAKE_INSTALL_PREFIX}/bin

install(FILES "testing.cpp" DESTINATION src)              # ${CMAKE_INSTALL_PREFIX}/bin 

# specify the C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_EXTENSIONS TRUE)

# CVODE libaries stuff
# specify libraries to link against
set(LIBDIR /usr/local/sundials-6/lib)
find_library(CVODE_LIB NAMES sundials_cvode HINTS "${LIBDIR}")
find_library(NVECS_LIB NAMES sundials_nvecserial HINTS "${LIBDIR}")
#find_library(MANYVECS_LIB NAMES sundials_nvecmanyvector  HINTS "${LIBDIR}")

set(SUNDIALS_LIBS ${CVODE_LIB} ${NVECS_LIB})
target_link_libraries(${PROJECT_NAME} PRIVATE ${SUNDIALS_LIBS})


# CVODE include libaries stuff
target_include_directories(${PROJECT_NAME} PRIVATE /usr/local/sundials-6/include)