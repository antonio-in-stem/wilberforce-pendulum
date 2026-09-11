# =============================================================================
#  Wilberforce Pendulum — CMake build
#
#     cmake -S . -B build
#     cmake --build build
#     ./build/wilberforce 14            # or build\Debug\wilberforce.exe on MSVC
#     ctest --test-dir build            # run the regression/physics tests
# =============================================================================
cmake_minimum_required(VERSION 3.16)
project(WilberforcePendulum
        VERSION 2.0.0
        DESCRIPTION "Numerical simulation of the Wilberforce pendulum (RK4, Fortran)"
        LANGUAGES Fortran)

set(CMAKE_Fortran_STANDARD 2008)
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_Fortran_FLAGS_RELEASE "-O2")
add_compile_options(-std=f2008 -Wall)

# --- library: the numerical core (CMake resolves Fortran module order) --------
add_library(wilberforce_core
    src/wilberforce_kinds.f90
    src/wilberforce_types.f90
    src/wilberforce_model.f90
    src/wilberforce_analytic.f90
    src/wilberforce_energy.f90
    src/wilberforce_io.f90
    src/wilberforce_integrator.f90
    src/wilberforce_presets.f90
    src/wilberforce_gnuplot.f90)

# --- driver executable --------------------------------------------------------
add_executable(wilberforce app/main.f90)
target_link_libraries(wilberforce PRIVATE wilberforce_core)

# --- test suite ---------------------------------------------------------------
enable_testing()
add_executable(test_wilberforce test/test_wilberforce.f90)
target_link_libraries(test_wilberforce PRIVATE wilberforce_core)
add_test(NAME physics_and_regression COMMAND test_wilberforce)
