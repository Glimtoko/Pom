cmake_minimum_required(VERSION 3.30)

project(
    Pom
    VERSION 1.1
    LANGUAGES CXX
)

# find_package(OpenMP REQUIRED)
add_compile_options(-fopenmp)
add_link_options(-fopenmp)

# Source code
add_subdirectory(src)