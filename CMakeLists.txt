cmake_minimum_required(VERSION 3.8)
project(thesis_project)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")

find_package(PkgConfig)

pkg_check_modules(HTSLIB REQUIRED htslib)
set(SOURCE_FILES main.cpp)
add_executable(thesis_project ${SOURCE_FILES})
target_link_libraries(thesis_project ${HTSLIB_LIBRARIES} pthread)
target_include_directories(thesis_project PUBLIC ${HTSLIB_INCLUDE_DIRS})
target_compile_options(thesis_project PUBLIC ${HTSLIB_CFLAGS_OTHER})