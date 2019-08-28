CMAKE_MINIMUM_REQUIRED( VERSION 2.6)


PROJECT(RichSolver)


ENABLE_LANGUAGE(Fortran)


#add new cmake module files taken from
#https://github.com/SethMMorton/cmake_fortran_template
#
SET(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/Modules/")

INCLUDE(${CMAKE_MODULE_PATH}FindOpenMP_Fortran.cmake)
INCLUDE(${CMAKE_MODULE_PATH}GetGitRevisionDescription.cmake)
INCLUDE(${CMAKE_MODULE_PATH}SetCompileFlag.cmake)
INCLUDE(${CMAKE_MODULE_PATH}SetFortranFlags.cmake)
INCLUDE(${CMAKE_MODULE_PATH}SetParallelizationLibrary.cmake)
INCLUDE(${CMAKE_MODULE_PATH}SetUpLAPACK.cmake)

# SET_COMPILE_FLAGS(CMAKE_C_FLAGS "${CMAKE_C_FLAGS}" C REQUIRED
#                   "-Wall"     # GNU
#                   "-warn all" # Intel
#                  )


# make sure that the default is a RELEASE
if (NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE RELEASE CACHE STRING
      "Choose the type of build, options are: None Debug Release."
      FORCE)
endif (NOT CMAKE_BUILD_TYPE)

# default installation
get_filename_component (default_prefix ".." ABSOLUTE)
set (CMAKE_INSTALL_PREFIX ${default_prefix} CACHE STRING
      "Choose the installation directory; by default it installs in the NORMA directory."
      FORCE)

# FFLAGS depend on the compiler
get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)

if (Fortran_COMPILER_NAME MATCHES "gfortran.*")
  # gfortran
  set (CMAKE_Fortran_FLAGS_RELEASE "-funroll-all-loops -fno-f2c -O3")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-fno-f2c -O0 -g")
elseif (Fortran_COMPILER_NAME MATCHES "ifort.*")
  # ifort (untested)
  set (CMAKE_Fortran_FLAGS_RELEASE "-f77rtl -O3")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-f77rtl -O0 -g")
elseif (Fortran_COMPILER_NAME MATCHES "g77")
  # g77
  set (CMAKE_Fortran_FLAGS_RELEASE "-funroll-all-loops -fno-f2c -O3 -m32")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-fno-f2c -O0 -g -m32")
else (Fortran_COMPILER_NAME MATCHES "gfortran.*")
  message ("CMAKE_Fortran_COMPILER full path: " ${CMAKE_Fortran_COMPILER})
  message ("Fortran compiler: " ${Fortran_COMPILER_NAME})
  message ("No optimized Fortran compiler flags are known, we just try -O2...")
  set (CMAKE_Fortran_FLAGS_RELEASE "-O2")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g")
endif (Fortran_COMPILER_NAME MATCHES "gfortran.*")

add_library(core_model  src/soil_states/mass_soil_state.F90
                        src/soil_states/head_soil_state.F90
                        src/abstract_classes/abstract_run_options.F90
                        src/abstract_classes/abstract_sim_control.F90
                        src/abstract_classes/abstract_soil_parameter.F90
                        src/abstract_classes/abstract_soil_state.F90
                        src/utils/constants.F90
                        src/utils/matrix_utils.F90
                        src/soil_parameters/BrooksCorey_soil_patameters.F90
                        src/soil_parameters/ClappHorn_soil_patameters.F90 )

add_executable(mod_sim  src/mod_control/controller.F90
                        src/mod_control/run_options.F90
                        src/mod_control/run_mod.F90)
target_link_libraries(mod_sim core_model)


add_executable(pan_sim  src/pan_control/controller.F90
                        src/pan_control/run_options.F90
                        src/pan_control/run_pan.F90)
target_link_libraries(pan_sim core_model)


add_executable(all_sim  src/pan_control/controller.F90
                        src/pan_control/run_options.F90
                        src/mod_control/controller.F90
                        src/mod_control/run_options.F90
                        src/run_all_sims.F90)
target_link_libraries(all_sim core_model)



# install executables and scripts
install (TARGETS mod_sim pan_sim all_sim
         RUNTIME DESTINATION "bin")

