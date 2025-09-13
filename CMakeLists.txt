cmake_minimum_required(VERSION 3.20)

project(
    Pom
    VERSION 1.1
    LANGUAGES CXX
)

# find_package(OpenMP REQUIRED)
add_compile_options(-fopenmp -O3)
add_link_options(-fopenmp -O3)

# Source code
add_subdirectory(src)