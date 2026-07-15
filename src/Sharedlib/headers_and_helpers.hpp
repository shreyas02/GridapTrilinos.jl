#ifndef _HEADERS_AND_HELPERS_
#define _HEADERS_AND_HELPERS_

// std
#include <algorithm>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdlib.h>
#include <exception>
#include <filesystem>
#include <vector>

// MPI
#include <mpi.h>

// Solvers Amesos2 and Belos  
#include <Amesos2.hpp>
#include <BelosConfigDefs.hpp>
#include <BelosLinearProblem.hpp>
#include <BelosTpetraAdapter.hpp>
#include <BelosBlockGmresSolMgr.hpp>
#include <BelosBlockCGSolMgr.hpp>
#include <BelosBiCGStabSolMgr.hpp>

// Teuchos
#include <Teuchos_Array.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_ScalarTraits.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_StandardCatchMacros.hpp>
#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_StackedTimer.hpp>
#include <Teuchos_Tuple.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>
#include <Teuchos_AbstractFactoryStd.hpp>

// Thyra 
#include <Thyra_LinearOpWithSolveBase.hpp>
#include <Thyra_VectorBase.hpp>
#include <Thyra_SolveSupportTypes.hpp>
#include <Thyra_LinearOpWithSolveBase.hpp>
#include <Thyra_LinearOpWithSolveFactoryHelpers.hpp>
#include <Thyra_BelosLinearOpWithSolveFactory.hpp>
#include <Thyra_TpetraLinearOp.hpp>
#include <Thyra_TpetraMultiVector.hpp>
#include <Thyra_TpetraVector.hpp>
#include <Thyra_TpetraThyraWrappers.hpp>
#include <Thyra_VectorBase.hpp>
#include <Thyra_VectorStdOps.hpp>
#include <Thyra_VectorSpaceBase_def.hpp>
#include <Thyra_VectorSpaceBase_decl.hpp>
#include <Thyra_Ifpack2PreconditionerFactory_def.hpp>

// Stratimikos
#include <Stratimikos_DefaultLinearSolverBuilder.hpp>
#include <Stratimikos_FROSch_def.hpp>
#include <Stratimikos_LinearSolverBuilder.hpp>

// Tpetra 
#include <Tpetra_Core.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_Map.hpp>
#include <Tpetra_MultiVector.hpp>
#include <Tpetra_Vector.hpp>
#include <Tpetra_Version.hpp>

// Kokkos 
#include <Kokkos_Core.hpp> 

// jlcxx
#include <jlcxx/jlcxx.hpp>
#include <jlcxx/array.hpp>
#include <jlcxx/stl.hpp>

// namespaces
using namespace Belos;
using namespace std;
using namespace Teuchos;

// typeDefs

// Basic Tpetra types with default template parameters
typedef Tpetra::CrsMatrix<> crs_matrix_type;
typedef Tpetra::Map<> Tpetra_map;
typedef Tpetra::MultiVector<> multivec_type;
typedef Tpetra::Vector<> vec_type;
typedef multivec_type::scalar_type scalar_type;
typedef multivec_type::local_ordinal_type local_ordinal_type;
typedef multivec_type::global_ordinal_type global_ordinal_type; 
typedef crs_matrix_type::node_type node_type;

// Define the scalar traits and magnitude type
typedef Teuchos::ScalarTraits<scalar_type> scalar_traits;
typedef typename scalar_traits::magnitudeType magnitude_type;

// Define Tpetra operator and multivector types explicitly using scalar_type
typedef Tpetra::Operator<scalar_type> tpetra_operator;
typedef Tpetra::MultiVector<scalar_type> multivec;

inline void enableIfpack2Preconditioner(Stratimikos::DefaultLinearSolverBuilder& builder)
{
  const RCP<const ParameterList> precValidParams =
    sublist(builder.getValidParameters(), "Preconditioner Types");

  if (!precValidParams->isParameter("Ifpack2")) {
    using Base = Thyra::PreconditionerFactoryBase<scalar_type>;
    using Impl = Thyra::Ifpack2PreconditionerFactory<crs_matrix_type>;
    builder.setPreconditioningStrategyFactory(abstractFactoryStd<Base, Impl>(), "Ifpack2");
  }
}

#endif
