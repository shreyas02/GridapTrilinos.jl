#ifndef _TRILINOS_TYPES_
#define _TRILINOS_TYPES_

#include "headers_and_helpers.hpp"

struct SolverResult {
  std::string name;
  int num_iters = 0;
  double residual = 0.0;
  double solve_time = 0.0;
  int verbose = 0;
  int depth = 0;
};

struct TpetraMatrixData {
  RCP<const Tpetra_map> rowMap;
  RCP<const Tpetra_map> colMap;
  RCP<crs_matrix_type> matrix;
  std::vector<size_t> rowptr;
  std::vector<local_ordinal_type> colind;
  std::vector<double> values;
};

struct TpetraVectorData {
  RCP<vec_type> vector;
};

struct TrilinosSolveData {
  RCP<vec_type> x;
  SolverResult result;
};

TrilinosSolveData TrilinosSolve(
  const RCP<crs_matrix_type>& A,
  const RCP<vec_type>& b,
  const std::string& parameterFilePath,
  bool verbose);

#endif
