cmake_minimum_required(VERSION 3.20)
project(SuperLUExample LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# # check for SuperLU
# set(SUPERLU_DIR "${PROJECT_SOURCE_DIR}/superlu")
# include_directories("${SUPERLU_DIR}/SRC")
# if(NOT EXISTS "${SUPERLU_DIR}/CMakeLists.txt")
#   message(
#     FATAL_ERROR
#     "The submodules were not downloaded. You have to do it manually!"
#   )
# endif()

# # build the SuperLU (and CBLAS if needed) static libs!
# add_subdirectory(${SUPERLU_DIR})


# Find SuperLU using find_package
find_package(PkgConfig REQUIRED)
pkg_check_modules(SUPERLU REQUIRED superlu)

if(NOT SUPERLU_FOUND)
    message(FATAL_ERROR "SuperLU not found. Install it with: pacman -S mingw-w64-ucrt-x86_64-superlu")
endif()

# Find OpenBLAS (use pkg-config for consistency in MSYS2)
pkg_check_modules(OpenBLAS REQUIRED openblas)

if(NOT OpenBLAS_FOUND)
    message(FATAL_ERROR "OpenBLAS not found. Install it with: pacman -S mingw-w64-ucrt-x86_64-openblas")
endif()


# Add executable
add_executable(superlu_example src/main.cpp)

# Link SuperLU and include directories
target_include_directories(superlu_example PRIVATE ${SUPERLU_INCLUDE_DIRS} ${OpenBLAS_INCLUDE_DIRS})
target_link_libraries(superlu_example PRIVATE ${SUPERLU_LIBRARIES} ${OpenBLAS_LIBRARIES})
