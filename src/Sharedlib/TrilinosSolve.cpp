#include "TrilinosTypes.hpp"

// Matrix-only setup path. This is the expensive part for direct solvers and
// preconditioners: Thyra/Stratimikos builds the LinearOpWithSolve and stores
// any symbolic/numeric setup state inside it.
TrilinosSolverCache TrilinosSolve(
  const RCP<crs_matrix_type>& A,
  const std::string& parameterFilePath,
  bool verbose) {

  TrilinosSolverCache solverCache;

  // Accept both plain Stratimikos XML files and files where the solver list is
  // nested under "Thyra Solver List".
  RCP<ParameterList> parameterList = getParametersFromXmlFile(parameterFilePath);
  RCP<ParameterList> thyraList = parameterList;
  if (parameterList->isSublist("Thyra Solver List")) {
    thyraList = sublist(parameterList, "Thyra Solver List");
  }
  solverCache.linearSolverName = thyraList->get<std::string>("Linear Solver Type");
  solverCache.preconditionerName = thyraList->get("Preconditioner Type", "None");

  // Wrap the existing Tpetra matrix; the LOWS object owns the setup/factorization state.
  auto thyraA =
    Thyra::createConstLinearOp<scalar_type, local_ordinal_type, global_ordinal_type, node_type>(A);

  // Build the generic Stratimikos solver factory from the Thyra XML list.
  Stratimikos::DefaultLinearSolverBuilder builder;
  Stratimikos::enableFROSch<scalar_type, local_ordinal_type, global_ordinal_type, node_type>(builder);
  builder.setParameterList(thyraList);
  auto solverFactory = Thyra::createLinearSolveStrategy(builder);
  solverCache.solver = Thyra::linearOpWithSolve<scalar_type>(*solverFactory, thyraA);
  (void)verbose;

  return solverCache;
}

// RHS-only solve path. The cached LOWS object still points at the matrix used
// during setup, so repeated calls can reuse any factorization/preconditioner
// state while only rebuilding Thyra wrappers for b and x.
TrilinosSolveData TrilinosSolve(
  const TrilinosSolverCache& solverCache,
  const RCP<vec_type>& b,
  bool verbose) {

  TrilinosSolveData data;
  Teuchos::Time timer("solve_time");

  // Allocate the solution vector in the same Tpetra space as b.
  RCP<vec_type> x(new vec_type(b->getMap()));
  x->putScalar(Teuchos::ScalarTraits<scalar_type>::zero());

  auto thyraB =
    Thyra::createConstMultiVector<scalar_type, local_ordinal_type, global_ordinal_type, node_type>(
      rcp_implicit_cast<const multivec_type>(b));
  auto thyraX =
    Thyra::createMultiVector<scalar_type, local_ordinal_type, global_ordinal_type, node_type>(
      rcp_implicit_cast<multivec_type>(x));

  // Solve A*x = b through Thyra, reusing any setup held by the LOWS object.
  timer.start();
  Thyra::SolveStatus<scalar_type> solveStatus =
    Thyra::solve<scalar_type>(*solverCache.solver, Thyra::NOTRANS, *thyraB, thyraX.ptr());
  const bool success = (solveStatus.solveStatus == Thyra::SOLVE_STATUS_CONVERGED);
  if (success) {
    if (verbose)
      std::cout << "\nEnd Result: Problem Solved!" << std::endl;
  } else {
    if (verbose)
      std::cout << "\nEnd Result: Error in solving" << std::endl;
  }
  timer.stop();

  // Store a small Julia-facing solve log.
  data.result.name =
    solverCache.linearSolverName + " + " + solverCache.preconditionerName + " Preconditioner";
  if (!solveStatus.extraParameters.is_null()) {
    data.result.num_iters = solveStatus.extraParameters->get("Iteration Count", 0);
  }
  data.result.residual =
    solveStatus.achievedTol == Thyra::SolveStatus<scalar_type>::unknownTolerance()
      ? 0.0
      : solveStatus.achievedTol;
  data.result.solve_time = timer.totalElapsedTime();
  data.x = x;

  return data;
}

// Compatibility path for callers that do not keep a numerical setup cache.
TrilinosSolveData TrilinosSolve(
  const RCP<crs_matrix_type>& A,
  const RCP<vec_type>& b,
  const std::string& parameterFilePath,
  bool verbose) {

  TrilinosSolverCache solverCache = TrilinosSolve(A, parameterFilePath, verbose);
  return TrilinosSolve(solverCache, b, verbose);
}
