# the following are directories where stuff will be installed to
set(INCLUDE_INSTALL_DIR      "${CMAKE_INSTALL_PREFIX}/include"             CACHE PATH "The subdirectory to the header prefix" FORCE)

set(eigen_header_SRCS 
  linearsolver.h
  linearsolverbase.h
  ludecomposition.h
  ludecompositionbase.h
  matrix.h
  matrixbase.h
  projective.h
  util.h
  vector.h
  vectorbase.h
  regressioninternal.h
  regression.h
  )

install(FILES 
  ${eigen_header_SRCS}
  DESTINATION ${INCLUDE_INSTALL_DIR}/eigen
  )

