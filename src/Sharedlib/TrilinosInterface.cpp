#include "TrilinosTypes.hpp"

template<> struct jlcxx::IsMirroredType<SolverResult> : std::false_type { };

template<> struct jlcxx::IsMirroredType<TpetraMatrixData> : std::false_type { };

template<> struct jlcxx::IsMirroredType<TpetraVectorData> : std::false_type { };

template<> struct jlcxx::IsMirroredType<TrilinosSolveData> : std::false_type { };

// Function to construct Tpetra matrix.

TpetraMatrixData ConstructTpetraMatrix(
  jlcxx::ArrayRef<double> A_nzval,
  jlcxx::ArrayRef<int64_t> A_colind,
  jlcxx::ArrayRef<int64_t> A_rowptr,
  jlcxx::ArrayRef<int64_t> RowPartition,
  jlcxx::ArrayRef<int64_t> ColPartition,
  int64_t LinSysSize,
  int64_t LocRowSize,
  jlcxx::ArrayRef<int32_t> OwnToValRow,
  const RCP<const Comm<int>>& comm) {

  TpetraMatrixData data;
  const Tpetra::global_size_t numGblIndices = LinSysSize;
  const global_ordinal_type indexBase = 0;

  std::vector<global_ordinal_type> ownedRowIndices(static_cast<size_t>(LocRowSize));
  for(int64_t i = 0; i < LocRowSize; ++i){
    ownedRowIndices[static_cast<size_t>(i)] =
      static_cast<global_ordinal_type>(RowPartition[OwnToValRow[i] - 1] - 1);
  }

  std::vector<global_ordinal_type> colPartition(static_cast<size_t>(ColPartition.size()));
  std::transform(ColPartition.begin(), ColPartition.end(), colPartition.begin(),
    [](int64_t value) { return static_cast<global_ordinal_type>(value - 1); });

  Teuchos::ArrayView<const global_ordinal_type> rowIndices(
    ownedRowIndices.data(),
    ownedRowIndices.size());
    
  Teuchos::ArrayView<const global_ordinal_type> colIndices(
    colPartition.data(),
    colPartition.size());

  data.rowMap = rcp(new Tpetra_map(numGblIndices, rowIndices, indexBase, comm));
  data.colMap = rcp(new Tpetra_map(numGblIndices, colIndices, indexBase, comm));

  const size_t rowptrSize = static_cast<size_t>(LocRowSize + 1);
  const size_t nnzOwnedRows = static_cast<size_t>(A_rowptr[LocRowSize]);

  data.rowptr.resize(rowptrSize);
  std::transform(A_rowptr.begin(), A_rowptr.begin() + rowptrSize, data.rowptr.begin(),
    [](int64_t value) { return static_cast<size_t>(value); });

  data.colind.resize(nnzOwnedRows);
  std::transform(A_colind.begin(), A_colind.begin() + nnzOwnedRows, data.colind.begin(),
    [](int64_t value) { return static_cast<local_ordinal_type>(value); });

  data.values.resize(nnzOwnedRows);
  std::copy(A_nzval.begin(), A_nzval.begin() + nnzOwnedRows, data.values.begin());

  Teuchos::ArrayRCP<size_t> rowptrView(data.rowptr.data(), 0, data.rowptr.size(), false);
  Teuchos::ArrayRCP<local_ordinal_type> colindView(data.colind.data(), 0, data.colind.size(), false);
  Teuchos::ArrayRCP<double> valuesView(data.values.data(), 0, data.values.size(), false);

  data.matrix = rcp(new crs_matrix_type(data.rowMap, data.colMap, rowptrView, colindView, valuesView));
  data.matrix->fillComplete();

  return data;
}

// Julia-facing Tpetra matrix construction entry point.

TpetraMatrixData ConstructTpetraMatrixWrapper(
  jlcxx::ArrayRef<double> A_nzval,
  jlcxx::ArrayRef<int64_t> A_colind,
  jlcxx::ArrayRef<int64_t> A_rowptr,
  jlcxx::ArrayRef<int64_t> RowPartition,
  jlcxx::ArrayRef<int64_t> ColPartition,
  int64_t LinSysSize,
  int64_t LocRowSize,
  jlcxx::ArrayRef<int32_t> OwnToValRow) {

  MPI_Comm yourComm = MPI_COMM_WORLD;
  RCP<const Comm<int> > comm (new MpiComm<int> (yourComm));
  return ConstructTpetraMatrix(
    A_nzval,
    A_colind,
    A_rowptr,
    RowPartition,
    ColPartition,
    LinSysSize,
    LocRowSize,
    OwnToValRow,
    comm);
}

// Construct Tpetra vector.

TpetraVectorData ConstructTpetraVectorWrapper(
  jlcxx::ArrayRef<double> RhsVec,
  int64_t LocRowSize,
  jlcxx::ArrayRef<int32_t> OwnToValRow,
  const TpetraMatrixData& matrixData) {

  TpetraVectorData data;
  data.vector = rcp(new vec_type(matrixData.rowMap));
  for(int64_t i = 0; i < LocRowSize; ++i){
    local_ordinal_type row = static_cast<local_ordinal_type>(i);
    double value = RhsVec[OwnToValRow[i] - 1];
    data.vector->sumIntoLocalValue(row, value);
  }
  return data;
}

// Copy the solution back to the Gridap vector.

SolverResult CopySolutionWrapper(
  const TrilinosSolveData& solveData,
  jlcxx::ArrayRef<double> LocSoln,
  int64_t LocRowSize,
  jlcxx::ArrayRef<int32_t> OwnToValSol) {

  auto x_data_host = solveData.x->getLocalViewHost(Tpetra::Access::ReadOnly);
  for(size_t i = 0; i< LocRowSize; ++i){
    LocSoln[OwnToValSol[i] - 1] = x_data_host(i, 0);
  }
  return solveData.result;
}

// Julia-facing solve entry point.

TrilinosSolveData TrilinosSolveWrapper(
  const TpetraMatrixData& matrixData,
  const TpetraVectorData& vectorData,
  std::string parameterFilePath) {

  MPI_Comm yourComm = MPI_COMM_WORLD;
  RCP<const Comm<int> > comm (new MpiComm<int> (yourComm));
  const int myRank = comm->getRank ();
  const bool verbose = (myRank == 0); // Only print on rank 0

  return TrilinosSolve(
    matrixData.matrix,
    vectorData.vector,
    parameterFilePath,
    verbose);
}

void KokkosInitialize() {
  if (!Kokkos::is_initialized()) {
    Kokkos::initialize();
  }
}

void KokkosFinalize() {
  if (Kokkos::is_initialized() && !Kokkos::is_finalized()) {
    Kokkos::finalize();
  }
}

JLCXX_MODULE define_julia_module(jlcxx::Module& mod) {
    mod.add_type<SolverResult>("SolverResult")
        .method("num_iters", [](const SolverResult& r) { return r.num_iters;})
        .method("residual", [](const SolverResult& r) { return r.residual;})
        .method("solve_time", [](const SolverResult& r) { return r.solve_time;})
        .method("name", [](const SolverResult& r) { return r.name;})
        .method("verbose", [](const SolverResult& r) { return r.verbose;})
        .method("depth", [](const SolverResult& r) { return r.depth;});
    mod.add_type<TpetraMatrixData>("TpetraMatrixData");
    mod.add_type<TpetraVectorData>("TpetraVectorData");
    mod.add_type<TrilinosSolveData>("TrilinosSolveData");
    mod.method("ConstructTpetraMatrixWrapper", static_cast<TpetraMatrixData (*)(
        jlcxx::ArrayRef<double>,
        jlcxx::ArrayRef<int64_t>,
        jlcxx::ArrayRef<int64_t>,
        jlcxx::ArrayRef<int64_t>,
        jlcxx::ArrayRef<int64_t>,
        int64_t,
        int64_t,
        jlcxx::ArrayRef<int32_t>)>(&ConstructTpetraMatrixWrapper));
    mod.method("ConstructTpetraVectorWrapper", &ConstructTpetraVectorWrapper);
    mod.method("TrilinosSolveWrapper", static_cast<TrilinosSolveData (*)(
        const TpetraMatrixData&,
        const TpetraVectorData&,
        std::string)>(&TrilinosSolveWrapper));
    mod.method("CopySolutionWrapper", static_cast<SolverResult (*)(
        const TrilinosSolveData&,
        jlcxx::ArrayRef<double>,
        int64_t,
        jlcxx::ArrayRef<int32_t>)>(&CopySolutionWrapper));
    mod.method("KokkosInitialize", &KokkosInitialize);
    mod.method("KokkosFinalize", &KokkosFinalize);
}
