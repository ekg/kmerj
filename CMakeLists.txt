# Specify the minimum version for CMake

cmake_minimum_required(VERSION 3.1)

# Project's name
project(kmerj)
# We build using c++14
set(CMAKE_CXX_STANDARD 14)

# Use all standard-compliant optimizations
set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -g")
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -g")

# Use openmp for parallelism, but it's configured differently on OSX
find_package(OpenMP)
if (OPENMP_FOUND)
  if (${CMAKE_SYSTEM_NAME} MATCHES "Linux")
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS} -fopenmp")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} -fopenmp")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
  elseif (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    # assumes clang build
    # we can't reliably detect when we're using clang, so for the time being we assume
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS} -Xpreprocessor -fopenmp")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} -Xpreprocessor -fopenmp")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS} -lomp")
  endif()
endif()

# Set the output folder where your program will be created
set(CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR}/bin)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})
set(LIBRARY_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/lib)

# The following folder will be included
include_directories("${PROJECT_SOURCE_DIR}")

# Add external projects
include(${CMAKE_ROOT}/Modules/ExternalProject.cmake)

# TODO: We're using INSTALL_DIR very wrong. We *should* be actually installing
# the external projects into their prefixes and working with the installed
# files. Instead we're building but not installing them and trying to work with
# the non-installed build trees.
# 
# Hence the blanked out INSTALL_COMMANDs to suppress the install step.
#
# We need to NOT blank out UPDATE_COMMAND or we can never change the Git revision we point to.
# The cost of this is that we have to re-configure on every build if we do update.

# sdsl-lite (full build using its cmake config)
ExternalProject_Add(sdsl-lite
  GIT_REPOSITORY "https://github.com/simongog/sdsl-lite.git"
  GIT_TAG "d52aa9a71513d132e30c09491b5899af449ebb94"
  CMAKE_ARGS -DCMAKE_INSTALL_PREFIX:PATH=${INSTALL_DIR} # TODO ADD static build flag
  UPDATE_COMMAND ""
  INSTALL_COMMAND "")
ExternalProject_Get_property(sdsl-lite INSTALL_DIR)
set(sdsl-lite_INCLUDE "${INSTALL_DIR}/src/sdsl-lite-build/include")
set(sdsl-lite-divsufsort_INCLUDE "${INSTALL_DIR}/src/sdsl-lite-build/external/libdivsufsort/include")
set(sdsl-lite_LIB "${INSTALL_DIR}/src/sdsl-lite-build/lib")
set(sdsl-lite-divsufsort_LIB "${INSTALL_DIR}/src/sdsl-lite-build/external/libdivsufsort/lib")

# DYNAMIC (full build using its cmake config)
ExternalProject_Add(dynamic
  GIT_REPOSITORY "https://github.com/vgteam/DYNAMIC.git"
  GIT_TAG "c8a87a12aec9d71a2e4151027dbcedd673e1e277"
  # we don't actually install dynamic... it's header only
  #CMAKE_ARGS -DCMAKE_INSTALL_PREFIX:PATH=${INSTALL_DIR}
  UPDATE_COMMAND ""
  INSTALL_COMMAND ""
  BUILD_COMMAND ""
  CONFIGURE_COMMAND "")
ExternalProject_Get_property(dynamic INSTALL_DIR)
set(dynamic_INCLUDE "${INSTALL_DIR}/src/dynamic/include")

# gzip file reader, header only
ExternalProject_Add(gzipreader
  GIT_REPOSITORY "https://github.com/gatoravi/gzip_reader.git"
  GIT_TAG "0ef26c0399e926087f9d6c4a56067a7bf1fc4f5e"
  UPDATE_COMMAND ""
  INSTALL_COMMAND "")
ExternalProject_Get_property(gzipreader SOURCE_DIR)
set(gzipreader_INCLUDE "${SOURCE_DIR}")

# taywee's C++ args library, header only
ExternalProject_Add(tayweeargs
  GIT_REPOSITORY "https://github.com/Taywee/args.git"
  GIT_TAG "3de44ec671db452cc0c4ef86399b108939768abb"
  UPDATE_COMMAND ""
  INSTALL_COMMAND "")
ExternalProject_Get_property(tayweeargs SOURCE_DIR)
set(tayweeargs_INCLUDE "${SOURCE_DIR}")

# mmmultimap (memory mapped multimap)
ExternalProject_Add(mmmultimap
  GIT_REPOSITORY "https://github.com/ekg/mmmultimap.git"
  GIT_TAG "ecbc55799e72076286e19e774d842168db14e8bb"
  BUILD_COMMAND ""
  UPDATE_COMMAND ""
  INSTALL_COMMAND "")
ExternalProject_Get_property(mmmultimap SOURCE_DIR)
set(mmmultimap_INCLUDE "${SOURCE_DIR}/src")

# In-place Parallel Super Scalar Samplesort (IPS⁴o), header only
ExternalProject_Add(ips4o
  GIT_REPOSITORY "https://github.com/SaschaWitt/ips4o.git"
  GIT_TAG "bff3ccf0bf349497f2bb10f825d160b792236367"
  INSTALL_COMMAND ""
  BUILD_COMMAND ""
  CONFIGURE_COMMAND "")
ExternalProject_Get_property(ips4o SOURCE_DIR)
set(ips4o_INCLUDE "${SOURCE_DIR}")

# BBHash perfect hasher
ExternalProject_Add(bbhash
  GIT_REPOSITORY "https://github.com/vgteam/BBHash.git"
  GIT_TAG "ccc22805eae48aca2a01a7ff5ef8e002857020d3"
  UPDATE_COMMAND ""
  INSTALL_COMMAND ""
  BUILD_COMMAND ""
  CONFIGURE_COMMAND "")
ExternalProject_Get_property(bbhash SOURCE_DIR)
set(bbhash_INCLUDE "${SOURCE_DIR}")

# ska
ExternalProject_Add(ska
  GIT_REPOSITORY "https://github.com/skarupke/flat_hash_map.git"
  GIT_TAG "2c4687431f978f02a3780e24b8b701d22aa32d9c"
  UPDATE_COMMAND ""
  INSTALL_COMMAND ""
  BUILD_COMMAND ""
  CONFIGURE_COMMAND "")
ExternalProject_Get_property(ska SOURCE_DIR)
set(ska_INCLUDE "${SOURCE_DIR}")

ExternalProject_Add(picosha256
  GIT_REPOSITORY "https://github.com/okdshin/PicoSHA2.git"
  GIT_TAG "b699e6c900be6e00152db5a3d123c1db42ea13d0"
  UPDATE_COMMAND ""
  INSTALL_COMMAND ""
  BUILD_COMMAND ""
  CONFIGURE_COMMAND "")
ExternalProject_Get_property(picosha256 SOURCE_DIR)
set(picosha256_INCLUDE "${SOURCE_DIR}")

set(CMAKE_BUILD_TYPE Release)

# set up our target executable and specify its dependencies and includes
add_executable(kmerj
  ${CMAKE_SOURCE_DIR}/src/kmer.cpp
  ${CMAKE_SOURCE_DIR}/src/main.cpp
  )
add_dependencies(kmerj sdsl-lite)
add_dependencies(kmerj dynamic)
add_dependencies(kmerj tayweeargs)
add_dependencies(kmerj bbhash)
add_dependencies(kmerj mmmultimap)
add_dependencies(kmerj ips4o)
add_dependencies(kmerj ska)
add_dependencies(kmerj picosha256)
target_include_directories(kmerj PUBLIC
  "${CMAKE_SOURCE_DIR}/src"
  "${sdsl-lite_INCLUDE}"
  "${sdsl-lite-divsufsort_INCLUDE}"
  "${dynamic_INCLUDE}"
  "${gzipreader_INCLUDE}"
  "${tayweeargs_INCLUDE}"
  "${ska_INCLUDE}"
  "${ips4o_INCLUDE}"
  "${mmmultimap_INCLUDE}"
  "${bbhash_INCLUDE}"
  "${picosha256_INCLUDE}")
target_link_libraries(kmerj
  "${sdsl-lite_LIB}/libsdsl.a"
  "${sdsl-lite-divsufsort_LIB}/libdivsufsort.a"
  "${sdsl-lite-divsufsort_LIB}/libdivsufsort64.a"
  "-latomic"
  z)

if (APPLE)
elseif (TRUE)
  set(CMAKE_EXE_LINKER_FLAGS "-static")
endif()
