#include "TrilinosTypes.hpp"

TrilinosSolveData TrilinosSolve(
  const RCP<crs_matrix_type>& A,
  const RCP<vec_type>& b,
  const std::string& parameterFilePath,
  bool verbose) {

  TrilinosSolveData data;
  Teuchos::Time timer("solve_time");

  // Allocate the solution vector in the same Tpetra space as A.
  RCP<vec_type> x(new vec_type(A->getDomainMap()));
  x->putScalar(Teuchos::ScalarTraits<scalar_type>::zero());

  RCP<ParameterList> parameterList = getParametersFromXmlFile(parameterFilePath);
  RCP<ParameterList> thyraList = parameterList;
  if (parameterList->isSublist("Thyra Solver List")) {
    thyraList = sublist(parameterList,"Thyra Solver List");
  }
  const std::string linearSolverName = thyraList->get<std::string>("Linear Solver Type");
  const std::string preconditionerName = thyraList->get("Preconditioner Type", "None");

  // Wrap the existing Tpetra objects; Thyra does not own new matrix/vector data here.
  auto thyraA =
    Thyra::createConstLinearOp<scalar_type, local_ordinal_type, global_ordinal_type, node_type>(A);
  auto thyraB =
    Thyra::createConstMultiVector<scalar_type, local_ordinal_type, global_ordinal_type, node_type>(
      rcp_implicit_cast<const multivec_type>(b));
  auto thyraX =
    Thyra::createMultiVector<scalar_type, local_ordinal_type, global_ordinal_type, node_type>(
      rcp_implicit_cast<multivec_type>(x));

  // Build the generic Stratimikos solver factory from the Thyra XML list.
  Stratimikos::DefaultLinearSolverBuilder builder;
  Stratimikos::enableFROSch<scalar_type, local_ordinal_type, global_ordinal_type, node_type>(builder);
  builder.setParameterList(thyraList);
  auto solverFactory = Thyra::createLinearSolveStrategy(builder);
  auto solver = Thyra::linearOpWithSolve<scalar_type>(*solverFactory, thyraA);

  // Solve A*x = b through Thyra.
  timer.start();
  Thyra::SolveStatus<scalar_type> solveStatus =
    Thyra::solve<scalar_type>(*solver, Thyra::NOTRANS, *thyraB, thyraX.ptr());
  bool success = (solveStatus.solveStatus == Thyra::SOLVE_STATUS_CONVERGED);
  if (success) {
    if (verbose)
      std::cout << "\nEnd Result: Problem Solved!" << std::endl;
  } else {
    if (verbose)
      std::cout << "\nEnd Result: Error in solving" << std::endl;
  }
  timer.stop();

  // Store a small Julia-facing solve log.
  data.result.name = linearSolverName + " + " + preconditionerName + " Preconditioner";
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
