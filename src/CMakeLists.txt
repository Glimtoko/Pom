set(SRCFILES    mesh/mesh2d.cpp
                hydro/hydro_MH2D.cpp
                hydro/hydro_flux.cpp
                hydro/hydro_timestep.cpp
                main.cpp)

add_executable(pom ${SRCFILES})

# if(OpenMP_CXX_FOUND)
#     target_link_directories(pom PUBLIC OpenMP::OpenMP_CXX)
# endif()

target_link_options(pom PUBLIC "-fopenmp")