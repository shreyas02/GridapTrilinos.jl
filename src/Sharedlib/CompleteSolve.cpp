#include "headers_and_helpers.hpp"

struct SolverResult {
  std::string name;
  int num_iters = 0;
  double residual = 0.0;
  double solve_time = 0.0;
  std::vector<double> residuals{};
  int verbose = 0;
  int depth = 0;
};

template<> struct jlcxx::IsMirroredType<SolverResult> : std::false_type { };

SolverResult TrilinosParallel(jlcxx::ArrayRef<double> A_nzval, jlcxx::ArrayRef<int64_t> A_colind, 
  jlcxx::ArrayRef<int64_t> A_rowptr, jlcxx::ArrayRef<double> RhsVec, jlcxx::ArrayRef<double> LocSoln, jlcxx::ArrayRef<int64_t> RowPartition,
  jlcxx::ArrayRef<int64_t> ColPartition, int64_t LinSysSize, int64_t LocRowSize, int64_t LocColSize,
  jlcxx::ArrayRef<int32_t> OwnToValSol, jlcxx::ArrayRef<int32_t> OwnToValRow, int64_t maxNumEntriesPerRowInput,
  std::string parameterFilePath) {

  int argc = 0; char** argv = nullptr;

  // Setting the Global system size 
  Tpetra::global_size_t numGblIndices = LinSysSize;
  SolverResult result;

  MPI_Comm yourComm = MPI_COMM_WORLD;
  {
    // Passing MPI comm to Trilinos 
    RCP<const Comm<int> > comm (new MpiComm<int> (yourComm));
    const int myRank = comm->getRank ();
    const int numProcs = comm->getSize ();

    // Initialise FancyOStream
    RCP<Teuchos::FancyOStream> out = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
    const bool verbose = (myRank == 0); // Only print on rank 0 

    //// Seting up Teuchos timer
    Teuchos::Time timer("solve_time");

    //// Begin custom Tpetra map
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

    // Setting the IndexBase
    const global_ordinal_type indexBase = 0;

    // Row Map
    RCP<const Tpetra_map> rowMap =
            rcp(new Tpetra_map(numGblIndices, rowIndices, indexBase, comm));

    // Col Map
    RCP<const Tpetra_map> colMap =
            rcp(new Tpetra_map(numGblIndices, colIndices, indexBase, comm));
    //// End custom Tpetra map

    //// Begin Tpetra Matrix Assembly
    const size_t rowptrSize = static_cast<size_t>(LocRowSize + 1);
    const size_t nnzOwnedRows = static_cast<size_t>(A_rowptr[LocRowSize]);

    std::vector<size_t> rowptr(rowptrSize);
    std::transform(A_rowptr.begin(), A_rowptr.begin() + rowptrSize, rowptr.begin(),
      [](int64_t value) { return static_cast<size_t>(value); });

    std::vector<local_ordinal_type> colind(nnzOwnedRows);
    std::transform(A_colind.begin(), A_colind.begin() + nnzOwnedRows, colind.begin(),
      [](int64_t value) { return static_cast<local_ordinal_type>(value); });

    Teuchos::ArrayRCP<size_t> rowptrView(rowptr.data(), 0, rowptr.size(), false);
    Teuchos::ArrayRCP<local_ordinal_type> colindView(colind.data(), 0, colind.size(), false);
    Teuchos::ArrayRCP<double> valuesView(A_nzval.data(), 0, nnzOwnedRows, false);

    RCP<crs_matrix_type> A = rcp(new crs_matrix_type(rowMap, colMap, rowptrView, colindView, valuesView));
    A->fillComplete();
    size_t locSize = A->getLocalNumRows();

    //// End Tpetra Matrix assembly
    
    ///// Begin Initializing Solution Vector 
    RCP<vec_type> x(new vec_type(A->getDomainMap())); x->putScalar(Teuchos::ScalarTraits<scalar_type>::zero());
    ///// End Initializing Solution Vector

    //// Begin Right Hand side vector
    RCP<Tpetra::Vector<>> b = rcp(new Tpetra::Vector<>(rowMap));
    for(int64_t i = 0; i < LocRowSize; ++i){
      local_ordinal_type row = static_cast<local_ordinal_type>(i);
      double value = RhsVec[OwnToValRow[i] - 1];
      b->sumIntoLocalValue(row, value);
    }
    //// End Right Hand side vector

    //// Begin converting Tpetra objects to Xpetra objects
    RCP<Xpetra::Matrix<scalar_type, local_ordinal_type, global_ordinal_type, NO>> xpetra_A = MueLu::TpetraCrs_To_XpetraMatrix<scalar_type, local_ordinal_type, global_ordinal_type, NO>(A);
    RCP<Xpetra::Vector<scalar_type, local_ordinal_type, global_ordinal_type, NO>> xpetra_b = Xpetra::toXpetra(b);
    RCP<Xpetra::Vector<scalar_type, local_ordinal_type, global_ordinal_type, NO>> xpetra_x = Xpetra::toXpetra(x);
    //// End converting Tpetra objects to Xpetra objects

    //// Begin Reading parameters list from input XML file
    RCP<ParameterList> parameterList = getParametersFromXmlFile(parameterFilePath);
    RCP<ParameterList> belosList = sublist(parameterList,"Belos List");
    RCP<ParameterList> precList = sublist(parameterList,"Preconditioner List");
    //// End Reading parameters list from input XML file

    //// Begin constructing the preconditioner 
    RCP<onelevelpreconditioner_type> prec(new onelevelpreconditioner_type(xpetra_A,precList));
    prec->initialize(false);
    prec->compute(); 
    RCP<operatort_type> belosPrec = rcp(new xpetraop_type(prec));
    //// End constructing the preconditioner 

    //// Begin solving the system using Belos
    // Constructing the linear problem 
    RCP<operatort_type> belosA = rcp(new xpetraop_type(xpetra_A));
    RCP<linear_problem_type> linear_problem (new linear_problem_type(belosA,xpetra_x,xpetra_b));
    linear_problem->setProblem(xpetra_x,xpetra_b);
    if(locSize != 0){
    linear_problem->setRightPrec(belosPrec); // Specify the preconditioner
    } 
    // Sending the parameters to the solver
    solverfactory_type solverfactory;
    RCP<solver_type> solver = solverfactory.create(parameterList->get("Belos Solver Type","GMRES"),belosList);
    solver->setProblem(linear_problem);
    // Solve
    timer.start(); // Start timer for solve
    Belos::ReturnType ret = solver->solve();
    bool success = false; success = (ret==Belos::Converged);
    if (success) {
      if (verbose)
        std::cout << "\nEnd Result: Problem Solved!" << std::endl;
    } else {
      if (verbose)
        std::cout << "\nEnd Result: Error in solving" << std::endl;
    }
    timer.stop(); // End timer for solve
    //// End solving the system using Belos

    //// Begin to Capture results from the solver
    const std::string belos_solver = parameterList->get("Belos Solver Type", "GMRES");
    const std::string overlap_type = precList->get("OverlappingOperator Type", "UnknownPreconditioner");
    result.name = belos_solver + " + " + overlap_type + " Preconditioner";
    result.verbose = belosList->get("Verbosity", 0);
    result.depth   = belosList->get("Output Frequency", 0);
    result.num_iters = solver->getNumIters();
    result.residual = solver->achievedTol();
    result.solve_time = timer.totalElapsedTime();
    //// End to Capture results from the solver

    //// Begin copying the solution
    auto x_data_host = x->getLocalViewHost(Tpetra::Access::ReadOnly);
    for(size_t i = 0; i< LocRowSize; ++i){
      LocSoln[OwnToValSol[i] - 1] = x_data_host(i, 0);
    }
    //// End copying the solution
  }
  return result;
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
        .method("_residualsCxx", [](const SolverResult& r) { return r.residuals;})
        .method("name", [](const SolverResult& r) { return r.name;})
        .method("verbose", [](const SolverResult& r) { return r.verbose;})
        .method("depth", [](const SolverResult& r) { return r.depth;});
    mod.method("TrilinosParallel", &TrilinosParallel);
    mod.method("KokkosInitialize", &KokkosInitialize);
    mod.method("KokkosFinalize", &KokkosFinalize);
}
