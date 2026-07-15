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

// Thyra 
#include <Thyra_LinearOpWithSolveBase.hpp>
#include <Thyra_VectorBase.hpp>
#include <Thyra_SolveSupportTypes.hpp>
#include <Thyra_LinearOpWithSolveBase.hpp>
#include <Thyra_LinearOpWithSolveFactoryHelpers.hpp>
#include <Thyra_TpetraLinearOp.hpp>
#include <Thyra_TpetraMultiVector.hpp>
#include <Thyra_TpetraVector.hpp>
#include <Thyra_TpetraThyraWrappers.hpp>
#include <Thyra_VectorBase.hpp>
#include <Thyra_VectorStdOps.hpp>
#ifdef HAVE_SHYLU_DDFROSCH_EPETRA
#include <Thyra_EpetraLinearOp.hpp>
#endif
#include <Thyra_VectorSpaceBase_def.hpp>
#include <Thyra_VectorSpaceBase_decl.hpp>

// Xpetra
#include <Xpetra_Map.hpp>
#include <Xpetra_CrsMatrixWrap.hpp>
#include <Xpetra_DefaultPlatform.hpp>
#ifdef HAVE_SHYLU_DDFROSCH_EPETRA
#include <Xpetra_EpetraCrsMatrix.hpp>
#endif
#include <Xpetra_Parameters.hpp>

// Tpetra 
#include <Tpetra_Core.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_Map.hpp>
#include <Tpetra_MultiVector.hpp>
#include <Tpetra_Vector.hpp>
#include <Tpetra_Version.hpp>

// Galeri::Xpetra
#include "Galeri_XpetraProblemFactory.hpp"
#include "Galeri_XpetraMatrixTypes.hpp"
#include "Galeri_XpetraParameters.hpp"
#include "Galeri_XpetraUtils.hpp"
#include "Galeri_XpetraMaps.hpp"

// FROSch
#include <ShyLU_DDFROSch_config.h>
#include <FROSch_Tools_def.hpp>
#include <FROSch_SchwarzPreconditioners_fwd.hpp>
#include <FROSch_OneLevelPreconditioner_def.hpp>

// MueLU
#include <MueLu_Utilities.hpp>  
#include <MueLu_Utilities_decl.hpp>

// Kokkos 
#include <Kokkos_Core.hpp> 

// jlcxx
#include <jlcxx/jlcxx.hpp>
#include <jlcxx/array.hpp>
#include <jlcxx/stl.hpp>

// namespaces
using namespace Belos;
using namespace FROSch;
using namespace std;
using namespace Teuchos;
using namespace Xpetra;

// typeDefs

// Basic Tpetra types with default template parameters
typedef Tpetra::CrsMatrix<> crs_matrix_type;
typedef Tpetra::Map<> Tpetra_map;
typedef Tpetra::MultiVector<> multivec_type;
typedef Tpetra::Vector<> vec_type;
typedef multivec_type::scalar_type scalar_type;
typedef multivec_type::local_ordinal_type local_ordinal_type;
typedef multivec_type::global_ordinal_type global_ordinal_type; 
typedef MueLu::DefaultNode NO;

typedef MultiVector<double,int,FROSch::DefaultGlobalOrdinal,Tpetra::KokkosClassic::DefaultNode::DefaultNodeType> multivector_type;
// typedef multivector_type::scalar_type scalar_type;
// typedef multivector_type::local_ordinal_type local_ordinal_type;
// typedef multivector_type::global_ordinal_type global_ordinal_type;
typedef multivector_type::node_type node_type;
typedef MultiVectorFactory<scalar_type,local_ordinal_type,global_ordinal_type,node_type> multivectorfactory_type;
typedef Map<local_ordinal_type,global_ordinal_type,node_type> map_type;
typedef Matrix<scalar_type,local_ordinal_type,global_ordinal_type,node_type> matrix_type;
typedef CrsMatrixWrap<scalar_type,local_ordinal_type,global_ordinal_type,node_type> crsmatrixwrap_type;

typedef Galeri::Xpetra::Problem<Map<local_ordinal_type,global_ordinal_type,node_type>,crsmatrixwrap_type,multivector_type> problem_type;

typedef Belos::OperatorT<multivector_type> operatort_type;
typedef Belos::LinearProblem<scalar_type,multivector_type,operatort_type> linear_problem_type;
typedef Belos::SolverFactory<scalar_type,multivector_type,operatort_type> solverfactory_type;
typedef Belos::SolverManager<scalar_type,multivector_type,operatort_type> solver_type;
typedef XpetraOp<scalar_type,local_ordinal_type,global_ordinal_type,node_type> xpetraop_type;

typedef FROSch::OneLevelPreconditioner<scalar_type,local_ordinal_type,global_ordinal_type,node_type> onelevelpreconditioner_type;
typedef FROSch::TwoLevelPreconditioner<scalar_type,local_ordinal_type,global_ordinal_type,node_type> twolevelpreconditioner_type;

// Define the scalar traits and magnitude type
typedef Teuchos::ScalarTraits<scalar_type> scalar_traits;
typedef typename scalar_traits::magnitudeType magnitude_type;

// Define Tpetra operator and multivector types explicitly using scalar_type
typedef Tpetra::Operator<scalar_type> tpetra_operator;
typedef Tpetra::MultiVector<scalar_type> multivec;

#endif
