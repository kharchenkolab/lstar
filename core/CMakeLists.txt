cmake_minimum_required(VERSION 3.16)
project(lstar_core CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()
set(CMAKE_CXX_FLAGS_RELEASE "-O2")

# header-only core
add_library(lstar INTERFACE)
target_include_directories(lstar INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include)

# zlib enables reading gzip/zlib-compressed chunks (ubiquitous; optional).
# LIBZARR_HAS_ZLIB gates libzarr's gzip/zlib codec; LSTAR_HAVE_ZLIB is lstar's own guard.
find_package(ZLIB)
if(ZLIB_FOUND)
  target_link_libraries(lstar INTERFACE ZLIB::ZLIB)
  target_compile_definitions(lstar INTERFACE LSTAR_HAVE_ZLIB LIBZARR_HAS_ZLIB)
endif()

# zstd enables reading zstd-compressed chunks -- zarr-python 3's DEFAULT compressor, so it's what
# externally-produced v3 (and many v2) stores use. Optional: absent -> lstar reads gzip/none stores
# fine and throws a clear "compile with LIBZARR_HAS_ZSTD" error only on an actual zstd chunk.
find_library(ZSTD_LIBRARY NAMES zstd)
find_path(ZSTD_INCLUDE_DIR NAMES zstd.h)
if(ZSTD_LIBRARY AND ZSTD_INCLUDE_DIR)
  target_include_directories(lstar INTERFACE ${ZSTD_INCLUDE_DIR})
  target_link_libraries(lstar INTERFACE ${ZSTD_LIBRARY})
  target_compile_definitions(lstar INTERFACE LIBZARR_HAS_ZSTD)
  message(STATUS "lstar: zstd codec enabled (${ZSTD_LIBRARY})")
else()
  message(STATUS "lstar: zstd not found; zstd-compressed stores will be unreadable")
endif()

# OpenMP for the translation primitives.
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
  target_link_libraries(lstar INTERFACE OpenMP::OpenMP_CXX)
elseif(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  target_compile_options(lstar INTERFACE -fopenmp)
  target_link_libraries(lstar INTERFACE -fopenmp)
endif()

add_executable(test_crossimpl test/test_crossimpl.cpp)
target_link_libraries(test_crossimpl PRIVATE lstar)

# gcc < 9.1 needs an explicit link for <filesystem>
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 11)
  target_link_libraries(test_crossimpl PRIVATE stdc++fs)
endif()

add_executable(bench_colstats test/bench_colstats.cpp)
target_link_libraries(bench_colstats PRIVATE lstar)
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 11)
  target_link_libraries(bench_colstats PRIVATE stdc++fs)
endif()

add_executable(test_chunked test/test_chunked.cpp)
target_link_libraries(test_chunked PRIVATE lstar)
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 11)
  target_link_libraries(test_chunked PRIVATE stdc++fs)
endif()

add_executable(test_viewer test/test_viewer.cpp)
target_link_libraries(test_viewer PRIVATE lstar)
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 11)
  target_link_libraries(test_viewer PRIVATE stdc++fs)
endif()

add_executable(test_zip test/test_zip.cpp)
target_link_libraries(test_zip PRIVATE lstar)
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 11)
  target_link_libraries(test_zip PRIVATE stdc++fs)
endif()

add_executable(test_v3 test/test_v3.cpp)
target_link_libraries(test_v3 PRIVATE lstar)
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 11)
  target_link_libraries(test_v3 PRIVATE stdc++fs)
endif()

enable_testing()
add_test(NAME crossimpl COMMAND test_crossimpl)
add_test(NAME bench COMMAND bench_colstats)
add_test(NAME viewer COMMAND test_viewer)
