# Minimum cmake version
cmake_minimum_required (VERSION 3.0.0)

# Name of project and that it is C++ only.
project (ED_LINK CXX)

# ----------------------------------------------------------------------------
# Dependencies
# ----------------------------------------------------------------------------

# Set static Boost
set(Boost_USE_STATIC_LIBS ON)

# Search for zlib as a dependency for SeqAn.
find_package (ZLIB)

# Load the SeqAn module and fail if not found.
find_package (SeqAn REQUIRED)
find_package (Boost COMPONENTS system filesystem REQUIRED)
find_package (OpenMP REQUIRED)

# ----------------------------------------------------------------------------
# Build Setup
# ----------------------------------------------------------------------------

# Add include directories.
include_directories (${SEQAN_INCLUDE_DIRS})
include_directories (${BOOST_INCLUDE_DIR})

# Add definitions set by find_package (SeqAn).
add_definitions (${SEQAN_DEFINITIONS})

# Add CXX flags found by find_package (SeqAn).
#set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -pthread ${SEQAN_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -pthread -static ${SEQAN_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")

# Add processor executable.
add_executable(processor ED_links.cpp)

target_link_libraries (processor ${Boost_FILESYSTEM_LIBRARY} ${Boost_SYSTEM_LIBRARY})

# Add example executable.
include_directories(${CMAKE_SOURCE_DIR}/sdsl-lite/include)
link_directories(${CMAKE_SOURCE_DIR}/sdsl-lite/lib)
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ffast-math -funroll-loops -msse4.2 -march=native -msse4.2 -march=native -DHAVE_CXA_DEMANGLE")

add_executable(example expl-21.cpp RefDB.cc)

target_link_libraries (example sdsl divsufsort divsufsort64)

# Add generator executable.
add_executable (generator generator.cpp RefDB.cc)

target_link_libraries (generator ${SEQAN_LIBRARIES} rt ${Boost_FILESYSTEM_LIBRARY} ${Boost_SYSTEM_LIBRARY} sdsl divsufsort divsufsort64)

# Add converter executable
add_executable (converter L_F_convert.cpp)

target_link_libraries (converter ${Boost_FILESYSTEM_LIBRARY} ${Boost_SYSTEM_LIBRARY})

# Add ETS_profiler executable.
add_executable (profiler ETS_DB.cpp)

target_link_libraries (profiler ${Boost_FILESYSTEM_LIBRARY} ${Boost_SYSTEM_LIBRARY} sdsl divsufsort divsufsort64)

# Add OSS executable.
add_executable (OSS optimalSolverLN.cc testOptimalSolverLN.cc RefDB.cc)

target_link_libraries (OSS ${Boost_FILESYSTEM_LIBRARY} ${Boost_SYSTEM_LIBRARY} sdsl divsufsort divsufsort64)

# Add ETSS executable.
add_executable (ETSS testEtsSolver.cc etsSolver.cc)

target_link_libraries (ETSS ${Boost_FILESYSTEM_LIBRARY} ${Boost_SYSTEM_LIBRARY} sdsl divsufsort divsufsort64)

