## sources
set(SOURCES_F77 dgqvec.f dqhqr.f dqrinc.f dqrinr.f dqrot.f dqrqh.f dqrtv1.f)
set(SOURCES_F90 lassofun.f90 wrappers.f90)
set(SOURCES_Fortran_Main main.f90)
set(SOURCES_Fortran ${SOURCES_F77} ${SOURCES_F90})
set(SOURCES ${SOURCES_Fortran} ${SOURCES_Fortran_Main})

## gfortran-specific compilation flags
if ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "GNU")
  set_source_files_properties(${SOURCES} PROPERTIES COMPILE_FLAGS -fexternal-blas)
  set_source_files_properties(${SOURCES_F77} PROPERTIES COMPILE_FLAGS "-fallow-argument-mismatch -std=legacy")
endif()

## executable and library
add_library(lasso OBJECT ${SOURCES_Fortran})
add_executable(lasso_main ${SOURCES})
add_dependencies(lasso_main lasso)

## lapack
set_property(TARGET lasso APPEND PROPERTY LINK_FLAGS "${LAPACK_LINKER_FLAGS}")
target_link_libraries(lasso PRIVATE ${LAPACK_LIBRARIES})
set_property(TARGET lasso_main APPEND PROPERTY LINK_FLAGS "${LAPACK_LINKER_FLAGS}")
target_link_libraries(lasso_main PRIVATE ${LAPACK_LIBRARIES})

## install the binary
install(TARGETS lasso lasso_main RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

