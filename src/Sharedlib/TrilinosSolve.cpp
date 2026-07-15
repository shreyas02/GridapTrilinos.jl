#include "TrilinosTypes.hpp"

TrilinosSolveData TrilinosSolve(
  const RCP<crs_matrix_type>& A,
  const RCP<vec_type>& b,
  const std::string& parameterFilePath,
  bool verbose) {

  TrilinosSolveData data;
  Teuchos::Time timer("solve_time");

  RCP<vec_type> x(new vec_type(A->getDomainMap()));
  x->putScalar(Teuchos::ScalarTraits<scalar_type>::zero());

  const size_t locSize = A->getLocalNumRows();

  RCP<Xpetra::Matrix<scalar_type, local_ordinal_type, global_ordinal_type, NO>> xpetra_A = MueLu::TpetraCrs_To_XpetraMatrix<scalar_type, local_ordinal_type, global_ordinal_type, NO>(A);
  RCP<Xpetra::Vector<scalar_type, local_ordinal_type, global_ordinal_type, NO>> xpetra_b = Xpetra::toXpetra(b);
  RCP<Xpetra::Vector<scalar_type, local_ordinal_type, global_ordinal_type, NO>> xpetra_x = Xpetra::toXpetra(x);

  RCP<ParameterList> parameterList = getParametersFromXmlFile(parameterFilePath);
  RCP<ParameterList> belosList = sublist(parameterList,"Belos List");
  RCP<ParameterList> precList = sublist(parameterList,"Preconditioner List");

  RCP<onelevelpreconditioner_type> prec(new onelevelpreconditioner_type(xpetra_A,precList));
  prec->initialize(false);
  prec->compute();
  RCP<operatort_type> belosPrec = rcp(new xpetraop_type(prec));

  RCP<operatort_type> belosA = rcp(new xpetraop_type(xpetra_A));
  RCP<linear_problem_type> linear_problem (new linear_problem_type(belosA,xpetra_x,xpetra_b));
  linear_problem->setProblem(xpetra_x,xpetra_b);
  if(locSize != 0){
    linear_problem->setRightPrec(belosPrec);
  }

  solverfactory_type solverfactory;
  RCP<solver_type> solver = solverfactory.create(parameterList->get("Belos Solver Type","GMRES"),belosList);
  solver->setProblem(linear_problem);

  timer.start();
  Belos::ReturnType ret = solver->solve();
  bool success = (ret==Belos::Converged);
  if (success) {
    if (verbose)
      std::cout << "\nEnd Result: Problem Solved!" << std::endl;
  } else {
    if (verbose)
      std::cout << "\nEnd Result: Error in solving" << std::endl;
  }
  timer.stop();

  const std::string belos_solver = parameterList->get("Belos Solver Type", "GMRES");
  const std::string overlap_type = precList->get("OverlappingOperator Type", "UnknownPreconditioner");
  data.result.name = belos_solver + " + " + overlap_type + " Preconditioner";
  data.result.verbose = belosList->get("Verbosity", 0);
  data.result.depth = belosList->get("Output Frequency", 0);
  data.result.num_iters = solver->getNumIters();
  data.result.residual = solver->achievedTol();
  data.result.solve_time = timer.totalElapsedTime();
  data.x = x;

  return data;
}
