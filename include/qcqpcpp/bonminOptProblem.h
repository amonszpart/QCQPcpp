#ifndef BONMINOPT_H
#define BONMINOPT_H

#include "coin/BonTMINLP.hpp"
#include "Eigen/Dense"
#include "Eigen/Sparse"
#include "qcqpcpp/optProblem.h"

#warning ":LSJDF:LSDJ:LSDJFL:SJD"
using namespace  Ipopt;
using namespace Bonmin;

namespace qcqpcpp
{

class BonminOptException : public std::runtime_error
{
        using std::runtime_error::runtime_error;
};

template <typename _Scalar>
class BonminOpt : public qcqpcpp::OptProblem<_Scalar,int>
{
        typedef Eigen::Matrix<_Scalar,1,-1> VectorX;
        enum { NeedsToAlign = (sizeof(VectorX)%16)==0 };
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW_IF(NeedsToAlign)

        typedef typename qcqpcpp::OptProblem<_Scalar,int>               ParentType;
        typedef typename ParentType::SparseMatrix                       SparseMatrix;

        //! \brief BonminOpt    Default constructor.
        BonminOpt()
            : _printSol(false)
            , _debug   (false) {}

        //! \brief ~BonminOpt   virtual destructor.
        virtual ~BonminOpt() {}

        //! \brief BonminOpt    Copy constructor.
//        BonminOpt(const BonminOpt &other)
//            : _printSol(other._printSol)
//            , _debug(false) {}

        //! \brief operator=    Assignment operator. No data = nothing to assign.
        //MyTMINLP& operator=(const MyTMINLP&) {}

        /** \name Functions from OptProblem. */
        //@{
        //! \brief update               main Sets up problem tailored to solver. Bonminopt needs jacobian and hessian computation.
        //! \param verbose              Controls logging to console.
        virtual int update( bool verbose /* = false */ );

        //! \brief optimize             Runs optimization. Please call update before.
        //! \param x_out                Output values.
        //! \param objective_sense      Minimization or Maximization.
        virtual int optimize( std::vector<_Scalar> *x_out /* = NULL */, typename ParentType::OBJ_SENSE objecitve_sense /* = OBJ_SENSE::MINIMIZE */ );

        virtual _Scalar getINF() const override { return std::numeric_limits<_Scalar>::max(); } //DBL_MAX
        //@}

        static inline TMINLP::VariableType getVarTypeCustom( typename ParentType::VAR_TYPE var_type );
        inline SparseMatrix const&  getJacobian            ()        const { return _jacobian; }
        inline SparseMatrix const&  getHessian             ()        const { return _hessian; }
        inline VectorX      const&  getCachedqo            ()        const { return _qo; }
        inline SparseMatrix const&  getCachedQo            ()        const { return _Qo; }
        inline SparseMatrix const&  getCachedA             ()        const { return _A; }

        inline void printSolutionAtEndOfAlgorithm          ()              { _printSol = true; }
        inline bool isDebug                                ()        const { return _debug; }
        inline bool isPrintSol                             ()        const { return _printSol; }

    protected:
        SparseMatrix                 _jacobian; //!< \brief Cached Jacobian of linear constraints.
        SparseMatrix                 _hessian;  //!< \brief Cached Objective+quadConstriants Hessian, lower triangle only!

        VectorX                      _qo;       //!< \brief Cached full linear objective vector.
        SparseMatrix                 _Qo;       //!< \brief Cached full quadratic objective matrix.
        SparseMatrix                 _A;        //!< \brief Cached full linear constraint matrix.

    private:
        bool                          _printSol;  //!< \brief Flag, to print x in the end.
        bool                          _debug;

}; //...class BonminOpt

template <typename _Scalar>
class BonminTMINLP : public TMINLP
{
    public:
        typedef typename qcqpcpp::OptProblem<_Scalar,int>               ParentType;
        typedef typename ParentType::VectorX                            VectorX;
        typedef          Eigen::Map<const Eigen::Matrix<_Scalar,-1,1> > MatrixMapT;

        BonminTMINLP( BonminOpt<_Scalar> &delegate ) : _delegate( delegate ) {}

        //! \brief ~BonminOpt   virtual destructor.
        virtual ~BonminTMINLP() {}

        /** \name Overloaded functions specific to a TMINLP.*/
        //@{
        //! \brief get_variables_types  Pass the type of the variables (INTEGER, BINARY, CONTINUOUS) to the optimizer.
        //! \param n                    size of var_types (has to be equal to the number of variables in the problem)
        //! \param var_types            types of the variables (has to be filled by function).
        virtual bool get_variables_types( Index n, VariableType* var_types );

        //! \brief get_variables_linearity      Pass info about linear and nonlinear variables.
        virtual bool get_variables_linearity(Index n, Ipopt::TNLP::LinearityType* var_types);

        //! \brief get_constraints_linearity    Pass the type of the constraints (LINEAR, NON_LINEAR) to the optimizer.
        //! \param m                    Size of const_types (has to be equal to the number of constraints in the problem)
        //! \param const_types          Types of the constraints (has to be filled by function).
        virtual bool get_constraints_linearity(Index m, Ipopt::TNLP::LinearityType* const_types);
        //@}

        /** \name Overloaded functions defining a TNLP.
         * This group of function implement the various elements needed to define and solve a TNLP.
         * They are the same as those in a standard Ipopt NLP problem*/
        //@{
        //! \brief get_nlp_info         Method to pass the main dimensions of the problem to Ipopt.
        //! \param n                    Number of variables in problem.
        //! \param m                    Number of constraints.
        //! \param nnz_jac_g            Number of nonzeroes in Jacobian of constraints system.
        //! \param nnz_h_lag            Number of nonzeroes in Hessian of the Lagrangean.
        //! \param index_style          Indicate wether arrays are numbered from 0 (C-style) or from 1 (Fortran).
        //! \return true                in case of success.
        virtual bool get_nlp_info( Index& n, Index&m, Index& nnz_jac_g, Index& nnz_h_lag, TNLP::IndexStyleEnum& index_style );

        //! \brief get_bounds_info      Method to pass the bounds on variables and constraints to Ipopt.
        //! \param n                    Size of x_l and x_u (has to be equal to the number of variables in the problem)
        //! \param x_l                  Lower bounds on variables (function should fill it).
        //! \param x_u                  Upper bounds on the variables (function should fill it).
        //! \param m                    Size of g_l and g_u (has to be equal to the number of constraints in the problem).
        //! \param g_l                  Lower bounds of the constraints (function should fill it).
        //! \param g_u                  Upper bounds of the constraints (function should fill it).
        //! \return true                in case of success.
        virtual bool get_bounds_info( Index n, Number* x_l, Number* x_u, Index m, Number* g_l, Number* g_u );

        //! \brief get_starting_point   Method to to pass the starting point for optimization to Ipopt.
        //! \param init_x               Do we initialize primals?
        //! \param x                    Pass starting primal points (function should fill it if init_x is 1).
        //! \param m                    Size of lambda (has to be equal to the number of constraints in the problem).
        //! \param init_lambda          Do we initialize duals of constraints?
        //! \param lambda               Lower bounds of the constraints (function should fill it).
        //! \return true                in case of success.
        virtual bool get_starting_point( Index n, bool init_x, Number* x,
                                         bool init_z, Number* z_L, Number* z_U,
                                         Index m, bool init_lambda,
                                         Number* lambda );

        //! \brief eval_f               Method which compute the value of the objective function at point x.
        //! \param n                    size of array x (has to be the number of variables in the problem).
        //! \param x                    Point where to evaluate.
        //! \param new_x                Is this the first time we evaluate functions at this point? (in the present context we don't care).
        //! \param obj_value            Value of objective in x (has to be computed by the function).
        //! \return true                in case of success.
        virtual bool eval_f( Index n, const Number* x, bool new_x, Number& obj_value );

        //! \brief eval_grad_f          Method which compute the gradient of the objective at a point x.
        //! \param n                    Size of array x (has to be the number of variables in the problem).
        //! \param x                    Point where to evaluate.
        //! \param new_x                Is this the first time we evaluate functions at this point? (in the present context we don't care).
        //! \param grad_f               Gradient of objective taken in x (function has to fill it).
        //! \return true                in case of success.
        virtual bool eval_grad_f( Index n, const Number* x, bool new_x, Number* grad_f);

        //! \brief eval_g               Method which compute the value of the functions defining the constraints at a point x.
        //! \param n                    Size of array x (has to be the number of variables in the problem).
        //! \param x                    Point where to evaluate.
        //! \param new_x                Is this the first time we evaluate functions at this point? (in the present context we don't care).
        //! \param m                    Size of array g (has to be equal to the number of constraints in the problem)
        //! \param grad_f               Values of the constraints (function has to fill it).
        //! \return true                In case of success.
        virtual bool eval_g( Index n, const Number* x, bool new_x, Index m, Number* g);

        //! \brief eval_jac_g           Method to compute the Jacobian of the functions defining the constraints.
        //! \brief                      If the parameter values==NULL fill the arrays iCol and jRow which store the position of
        //! \brief                      the non-zero element of the Jacobian.
        //! \brief                      If the paramenter values!=NULL fill values with the non-zero elements of the Jacobian.
        //! \param n                    Size of array x (has to be the number of variables in the problem).
        //! \param x                    Point where to evaluate.
        //! \param new_x                Is this the first time we evaluate functions at this point? (in the present context we don't care).
        //! \param m                    Size of array g (has to be equal to the number of constraints in the problem)
        //! \param grad_f               Values of the constraints (function has to fill it).
        //! \return true                in case of success.
        virtual bool eval_jac_g( Index n, const Number* x, bool new_x,
                                 Index m, Index nele_jac, Index* iRow, Index *jCol,
                                 Number* values );

        //! \brief eval_h               Method to compute the Jacobian of the functions defining the constraints.
        //! \brief                      If the parameter values==NULL fill the arrays iCol and jRow which store the position of
        //! \brief                      the non-zero element of the Jacobian.
        //! \brief                      If the paramenter values!=NULL fill values with the non-zero elements of the Jacobian.
        //! \param n                    Size of array x (has to be the number of variables in the problem).
        //! \param x                    Point where to evaluate.
        //! \param new_x                Is this the first time we evaluate functions at this point? (in the present context we don't care).
        //! \param m                    Size of array g (has to be equal to the number of constraints in the problem)
        //! \param grad_f               Values of the constraints (function has to fill it).
        //! \return true                in case of success.
        virtual bool eval_h( Index n, const Number* x, bool new_x,
                             Number obj_factor, Index m, const Number* lambda,
                             bool new_lambda, Index nele_hess, Index* iRow,
                             Index* jCol, Number* values );


        //! \brief finalize_solution    Method called by Ipopt at the end of optimization.
        virtual void finalize_solution( TMINLP::SolverReturn status,
                                        Index n, const Number* x, Number obj_value );
        //@}

        virtual const SosInfo      * sosConstraints() const { return NULL; }
        virtual const BranchingInfo* branchingInfo () const { return NULL; }

    protected:
        BonminOpt<_Scalar>  &_delegate;
}; // ...class BonminTMINLP

} //...ns qcqpcpp

//_____________________________________________________________________________________________________________________
// HPP
#include <random>
#include "coin/BonBonminSetup.hpp"
#include "coin/BonCbc.hpp"          // Bab

#define MYDEBUG 1
namespace qcqpcpp
{

template <typename _Scalar> int
BonminOpt<_Scalar>::update( bool verbose /* = false */ )
{
    _jacobian = this->estimateJacobianOfConstraints();
    if ( isDebug() )
    {
        std::cout<<"[" << __func__ << "]: " << "jacobian ok" << _jacobian << std::endl; fflush(stdout);
    }

    _hessian = this->estimateHessianOfLagrangian();
    if ( isDebug() )
    {
        std::cout<<"[" << __func__ << "]: " << "hessian ok" << _hessian << std::endl; fflush(stdout);
    }

    // NOTE: _qo is a columnvector for convenient multiplication reasons, whilst sparsematrices are rowVectors by convention
    _qo = this->getLinObjectivesMatrix().transpose();
    _Qo = this->getQuadraticObjectivesMatrix();
    _A  = this->getLinConstraintsMatrix();

    this->_updated = true;

    return EXIT_SUCCESS;
}

template <typename _Scalar> int
BonminOpt<_Scalar>::optimize( std::vector<_Scalar> *x_out /* = NULL */, typename ParentType::OBJ_SENSE objecitve_sense /* = MINIMIZE */ )
{
    if ( !this->_updated )
    {
        std::cerr << "[" << __func__ << "]: " << "Please call update() first!" << std::endl;
        return EXIT_FAILURE;
    }

    // Now initialize from tminlp
    BonminSetup bonmin2;
    bonmin2.initializeOptionsAndJournalist();
    bonmin2.readOptionsString("bonmin.algorithm B-BB\n");
    // need this relay, otherwise, we'll end up with a double free corruption thing
    SmartPtr<BonminTMINLP<_Scalar> > problem = new BonminTMINLP<_Scalar>(*this);
    bonmin2.initialize( GetRawPtr(problem) );

    // Set up done, now let's branch and bound
    try
    {
        Bab bb;
        bb( bonmin2 ); // process parameter file using Ipopt and do branch and bound using Cbc
    }
    catch ( TNLPSolver::UnsolvedError *E )
    {
        // There has been a failure to solve a problem with Ipopt.
        std::cerr << "Ipopt has failed to solve a problem" << std::endl;
    }
    catch ( OsiTMINLPInterface::SimpleError &E )
    {
        std::cerr << E.className() << "::" << E.methodName()
                  << std::endl
                  << E.message()
                  << std::endl;
    }
    catch ( CoinError &E )
    {
        std::cerr << E.className() << "::" << E.methodName()
                  << std::endl
                  << E.message()
                  << std::endl;
    }
    catch ( std::exception &ex )
    {
        std::cerr << "[" << __func__ << "]: " << "exception: " << ex.what() << std::endl;
    }

    // output answer
    if ( x_out )
    {
        x_out->reserve( this->_x.size() );
        for ( int i = 0; i != this->_x.size(); ++i )
            x_out->push_back( this->_x[i] );
    }

    return EXIT_SUCCESS;
}

template <typename _Scalar> TMINLP::VariableType
BonminOpt<_Scalar>::getVarTypeCustom( typename ParentType::VAR_TYPE var_type )
{
    switch ( var_type )
    {
        case ParentType::VAR_TYPE::CONTINUOUS:
            return TMINLP::VariableType::CONTINUOUS;
            break;
        case ParentType::VAR_TYPE::INTEGER:
            return TMINLP::VariableType::INTEGER;
            break;
        case ParentType::VAR_TYPE::BINARY:
            return TMINLP::VariableType::BINARY;
            break;
        default:
            std::cerr << "[" << __func__ << "]: " << "Unrecognized file type, returning continuous" << std::endl;
            return TMINLP::VariableType::CONTINUOUS;
            break;
    } //... switch

} //...BonminOpt::getVarTypeCustom()

template <typename _Scalar> bool
BonminTMINLP<_Scalar>::get_variables_types( Index n, VariableType* var_types )
{
    if ( _delegate.isDebug() )
    {
        std::cout << "[" << __func__ << "]: " << "call(n = " << n << ")" << std::endl; fflush( stdout );
    }

    if ( n != _delegate.getVarCount() )
        throw new BonminOptException( "[BonminOpt::get_variables_types] n != getVarCount()" );

    for ( int j = 0; j != _delegate.getVarCount(); ++j )
        var_types[ j ] = BonminOpt<_Scalar>::getVarTypeCustom( _delegate.getVarType(j) );

    if ( _delegate.isDebug() )
    {
        std::cout << "[" << __func__ << "]: " << "finish" << std::endl; fflush( stdout );
    }

    return true;
} //...BonminOpt::get_variables_type()

template <typename _Scalar> bool
BonminTMINLP<_Scalar>::get_variables_linearity( Index n, Ipopt::TNLP::LinearityType* var_types )
{
    if ( _delegate.isDebug() )
    {
        std::cout << "[" << __func__ << "]: " << "call(n = " << n << ")" << std::endl; fflush( stdout );
    }

    if ( n != _delegate.getVarCount() )
        throw new BonminOptException( "[BonminOpt::get_variables_types] n != getVarCount()" );

    // NOTE: this is hard coded for now, non-lin cases not tested
    for ( int j = 0; j != _delegate.getVarCount(); ++j )
        var_types[j] = Ipopt::TNLP::LINEAR;

    if ( _delegate.isDebug() )
    {
        std::cout << "[" << __func__ << "]: " << "finish" << std::endl; fflush( stdout );
    }

    return true;
} //...BonminOpt::get_variables_linearity()


template <typename _Scalar> bool
BonminTMINLP<_Scalar>::get_constraints_linearity( Index m, Ipopt::TNLP::LinearityType* const_types )
{
    if ( _delegate.isDebug() )
    {
        std::cout << "[" << __func__ << "]: " << "call(m = " << m << ")" << std::endl; fflush( stdout );
    }

    if ( m != _delegate.getConstraintCount() )
        throw new BonminOptException( "[BonminOpt::get_variables_types] m != getConstraintCount()" );

    // NOTE: this is hard coded for now, non-lin cases not tested
    for ( int i = 0; i != _delegate.getConstraintCount(); ++i )
        const_types[ i ] = Ipopt::TNLP::LINEAR;

    if ( _delegate.isDebug() )
    {
        std::cout << "[" << __func__ << "]: " << "finish" << std::endl; fflush( stdout );
    }

    return true;
} //...BonminOpt::get_constraints_linearity()

template <typename _Scalar> bool
BonminTMINLP<_Scalar>::get_nlp_info( Index                & n
                                , Index                & m
                                , Index                & nnz_jac_g
                                , Index                & nnz_h_lag
                                , TNLP::IndexStyleEnum & index_style )
{
    if ( _delegate.isDebug() )
    {
        std::cout << "[" << __func__ << "]: " << "call(n = " << n
                  << ", m = " << m
                  << ")" << std::endl; fflush( stdout );
    }

    n           = _delegate.getVarCount();          // number of variable
    m           = _delegate.getConstraintCount();   // number of constraints
    nnz_jac_g   = _delegate.getJacobian().nonZeros();        // number of non zeroes in Jacobian
    nnz_h_lag   = _delegate.getHessian().nonZeros();          // number of non zeroes in Hessian of Lagrangean
    index_style = TNLP::C_STYLE;                // zero-indexed

    if ( _delegate.isDebug() )
    {
        std::cout << "[" << __func__ << "]: " << "finish" << std::endl; fflush( stdout );
    }

    return true;
} //...BonminOpt::get_nlp_info()

template <typename _Scalar> bool
BonminTMINLP<_Scalar>::get_bounds_info( Index    n
                                   , Number * x_l
                                   , Number * x_u
                                   , Index    m
                                   , Number * g_l
                                   , Number * g_u )
{
    if ( _delegate.isDebug() )
    {
        std::cout << "[" << __func__ << "]: " << "call(n = " << n
              << ", m = " << m
              << ")" << std::endl; fflush( stdout );
    }

    if ( n != _delegate.getVarCount() )
        throw new BonminOptException( "[BonminOpt::get_bounds_info] n != getVarCount()" );
    if ( m != _delegate.getConstraintCount() )
        throw new BonminOptException( "[BonminOpt::get_bounds_info] m != getConstraintCount()" );

    for ( int j = 0; j != _delegate.getVarCount(); ++j )
    {
        x_l[ j ] = _delegate.getVarLowerBound( j );
        x_u[ j ] = _delegate.getVarUpperBound( j );
    } //...for vars

    for ( int i = 0; i != _delegate.getConstraintCount(); ++i )
    {
        g_l[ i ] = _delegate.getConstraintLowerBound( i );
        g_u[ i ] = _delegate.getConstraintUpperBound( i );
    } //...for constraints

//    x_l[0] = 0;
//    x_u[0] = 1;

//    x_l[1] = 0;
//    x_u[1] = 1;

//    x_l[2] = 0;
//    x_u[2] = 1;

//    x_l[3] = 0;
//    x_u[3] = 1;

//    g_l[0] = 1.;
//    g_u[0] = DBL_MAX;

//    g_l[1] = 1.;
//    g_u[1] = DBL_MAX;

    if ( _delegate.isDebug() )
    {
        std::cout << "[" << __func__ << "]: " << "finish" << std::endl; fflush( stdout );
    }

    return true;
} //...BonminOpt::get_bounds_info()

template <typename _Scalar> bool
BonminTMINLP<_Scalar>::get_starting_point( Index    n
                                      , bool     init_x
                                      , Number * x
                                      , bool     init_z
                                      , Number * z_L
                                      , Number * z_U
                                      , Index    m
                                      , bool     init_lambda
                                      , Number * lambda )
{
    if ( _delegate.isDebug() )
    {
        std::cout << "[" << __func__ << "]: " << "call(n = " << n
                  << ", m = " << m
                  << ")" << std::endl; fflush( stdout );
    }

    if ( n != _delegate.getVarCount() )
        throw new BonminOptException( "[BonminOpt::get_starting_point] n != getVarCount()" );
    if ( m != _delegate.getConstraintCount() )
        throw new BonminOptException( "[BonminOpt::get_starting_point] m != getConstraintCount()" );

    if ( !init_x )
        throw new BonminOptException( "[BonminOpt::get_starting_point] !init_x..." );
    //assert( init_x );

    if ( init_lambda )
        throw new BonminOptException( "[BonminOpt::get_starting_point] init_lambda..." );
    //assert( !init_lambda );

    // random for now...
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<>         *int_distribution  = NULL;
    std::uniform_real_distribution<_Scalar> *real_distribution = NULL;
    for ( int j = 0; j != _delegate.getVarCount(); ++j )
    {
        if (    (_delegate.getVarType(j) == ParentType::VAR_TYPE::INTEGER)
             || (_delegate.getVarType(j) == ParentType::VAR_TYPE::BINARY ) )
        {
            if ( !int_distribution )
                int_distribution = new std::uniform_int_distribution<>( _delegate.getVarLowerBound(j), _delegate.getVarUpperBound(j) );
            x[ j ] = (*int_distribution)( gen );
        }
        else
        {
            if ( !real_distribution )
                real_distribution = new std::uniform_real_distribution<>( _delegate.getVarLowerBound(j), _delegate.getVarUpperBound(j) );

            x[ j ] = (*real_distribution)( gen );
        }
    } //...for variables

//    x[0] = 1;
//    x[1] = 1;
//    x[2] = 1;
//    x[3] = 1;

    // cleanup
    if ( int_distribution ) { delete int_distribution; int_distribution = NULL; }
    if ( real_distribution ) { delete real_distribution; real_distribution = NULL; }

    if ( _delegate.isDebug() )
    {
        std::cout << "[" << __func__ << "]: " << "finish" << std::endl;
        fflush( stdout );
    }

    return true;
} //...BonminOpt::get_starting_point()

template <typename _Scalar> bool
BonminTMINLP<_Scalar>::eval_f( Index n, const Number* x, bool new_x, Number& obj_value )
{
    if ( _delegate.isDebug() )
    {
        std::cout << "[" << __func__ << "]: " << "call(n = " << n
                  << ")" << std::endl;
        fflush( stdout );
    }

    if ( n != _delegate.getVarCount() )
        throw new BonminOptException( "[BonminOpt::eval_f] n != getVarCount()" );


    MatrixMapT x_eig( x, _delegate.getVarCount() );
    obj_value = (x_eig.transpose() * _delegate.getCachedQo() * x_eig + _delegate.getCachedqo() * x_eig).coeff( 0 );

    if ( _delegate.isDebug() )
    {
        std::cout << "[" << __func__ << "]: " << "finish" << std::endl;
        fflush( stdout );
    }

    return true;
} //...BonminOpt::eval_f()

template <typename _Scalar> bool
BonminTMINLP<_Scalar>::eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f)
{
    if ( _delegate.isDebug() )
    {
        std::cout << "[" << __func__ << "]: " << "call(n = " << n
                  << ")" << std::endl;
        fflush( stdout );
    }

    if ( n != _delegate.getVarCount() )
        throw new BonminOptException( "[BonminOpt::eval_grad_f] n != getVarCount()" );

  for ( int i = 0; i != _delegate.getVarCount(); ++i )
      grad_f[i] = _delegate.getCachedqo().coeff(i);

//  grad_f[0] = -1.;
//  grad_f[1] = -1.;
//  grad_f[2] = -1.;
//  grad_f[3] = 0.;

  if ( _delegate.isDebug() )
  {
      std::cout << "[" << __func__ << "]: " << "finish" << std::endl;
      fflush( stdout );
  }

  return true;
} //...BonminOpt::eval_grad_f()

template <typename _Scalar> bool
BonminTMINLP<_Scalar>::eval_g( Index n, const Number* x, bool new_x, Index m, Number* g )
{
    if ( _delegate.isDebug() )
    {
        std::cout << "[" << __func__ << "]: " << "call(n = " << n
                  << ", m = " << m
                  << ")" << std::endl;
        fflush( stdout );
    }

    if ( n != _delegate.getVarCount() )
        throw new BonminOptException( "[BonminOpt::eval_g] n != getVarCount()" );
    if ( m != _delegate.getConstraintCount() )
        throw new BonminOptException( "[BonminOpt::eval_g] m != getConstraintCount()" );

    MatrixMapT x_eig( x, _delegate.getVarCount() );
    Eigen::Matrix<_Scalar,-1,1> c = _delegate.getCachedA() * x_eig;
    for ( int i = 0; i != _delegate.getConstraintCount(); ++i )
    {
        g[i] = c(i);
    }

    //  g[0] = (x[1] - 1./2.)*(x[1] - 1./2.) + (x[2] - 1./2.)*(x[2] - 1./2.);
    //  g[1] = x[0] - x[1];
    //  g[2] = x[0] + x[2] + x[3];

    if ( _delegate.isDebug() )
    {
        std::cout << "[" << __func__ << "]: " << "finish" << std::endl;
        fflush( stdout );
    }

    return true;
} //...BonminOpt::eval_g()

template <typename _Scalar> bool
BonminTMINLP<_Scalar>::eval_jac_g( Index         n
                              , Number const* x
                              , bool          new_x
                              , Index         m
                              , Index         nnz_jac
                              , Index       * iRow
                              , Index       * jCol
                              , Number      * values )
{
    bool ret_val = false;

    if ( _delegate.isDebug() )
    {
        std::cout << "[" << __func__ << "]: " << "call(n = " << n
                  << ", m = " << m
                  << ")" << std::endl;
        fflush( stdout );
    }

    if ( n != _delegate.getVarCount() )
        throw new BonminOptException( "[BonminOpt::eval_jac_g] n != getVarCount()" );
    if ( nnz_jac != _delegate.getJacobian().nonZeros() )
        throw new BonminOptException( "[BonminOpt::eval_jac_g] nnz_jac != _jacobian.nonZeros()" );
//        assert( nnz_jac == 4 );

    int entry_id = 0;
    if ( values == NULL )
    {
        for ( int row = 0; row != _delegate.getJacobian().outerSize(); ++row )
        {
            for ( typename Eigen::SparseMatrix<_Scalar,Eigen::RowMajor>::InnerIterator it(_delegate.getJacobian(),row);
                  it; ++it, ++entry_id )
            {
                iRow[ entry_id ] = it.row();
                jCol[ entry_id ] = it.col();
            } // ...for col
        } // ...for row

        //    iRow[0] = 2;
        //    jCol[0] = 1;

        //    iRow[1] = 3;
        //    jCol[1] = 1;

        //    iRow[2] = 1;
        //    jCol[2] = 2;

        //    iRow[3] = 2;
        //    jCol[3] = 2;

        //    iRow[4] = 1;
        //    jCol[4] = 3;

        //    iRow[5] = 3;
        //    jCol[5] = 3;

        //    iRow[6] = 3;
        //    jCol[6] = 4;

        ret_val = true;
    } // ... if values == NULL
    else
    {
        // NOTE: ONLY for linear constraints
        for ( int row = 0; row != _delegate.getJacobian().outerSize(); ++row )
        {
            for ( typename Eigen::SparseMatrix<_Scalar,Eigen::RowMajor>::InnerIterator it(_delegate.getJacobian(),row);
                  it; ++it, ++entry_id )
            {
                values[ entry_id ] = it.value();
            } // ...for col
        } // ...for row

        //    values[0] = 1.;
        //    values[1] = 1;

        //    values[2] = 2*x[1] - 1;
        //    values[3] = -1.;

        //    values[4] = 2*x[2] - 1;
        //    values[5] = 1.;

        //    values[6] = 1.;

        ret_val = true;
    } // ... else values != NULL

    if ( _delegate.isDebug() )
    {
        std::cout << "[" << __func__ << "]: " << "finish" << std::endl;
        fflush( stdout );
    }

    return ret_val;
} // ... BonminOpt::eval_jac_g()

template <typename _Scalar> bool
BonminTMINLP<_Scalar>::eval_h( Index          n
                          , Number const * x
                          , bool           new_x
                          , Number         obj_factor
                          , Index          m
                          , Number const * lambda
                          , bool           new_lambda
                          , Index          nele_hess
                          , Index        * iRow
                          , Index        * jCol
                          , Number       * values )
{
    bool ret_val = false;

    if ( _delegate.isDebug() )
    {
        std::cout << "[" << __func__ << "]: " << "call(n = " << n
              << ", m = " << m
              << ")" << std::endl;
        fflush( stdout );
    }

    if ( n != _delegate.getVarCount() )
        throw new BonminOptException( "[BonminOpt::eval_h] n != getVarCount()" );
    if ( m != _delegate.getConstraintCount() )
        throw new BonminOptException( "[BonminOpt::eval_h] m != getConstraintCount()" );
    if ( nele_hess != _delegate.getHessian().nonZeros() )
        throw new BonminOptException( "[BonminOpt::eval_h] nele_hess != _delegate.getHessian().nonZeros()" );
//    assert( nele_hess == 5 );

    int entry_id = 0;
    if ( values == NULL )
    {
        for ( int row = 0; row != _delegate.getHessian().outerSize(); ++row )
        {
            for ( typename Eigen::SparseMatrix<_Scalar,Eigen::RowMajor>::InnerIterator it(_delegate.getHessian(),row);
                  it; ++it, ++entry_id )
            {
                iRow[ entry_id ] = it.row();
                jCol[ entry_id ] = it.col();
            } // ...for col
        } // ...for row

        ret_val = true;
    }
    else {
        // NOTE: lower triangular only please!
        for ( int row = 0; row != _delegate.getHessian().outerSize(); ++row )
        {
            for ( typename Eigen::SparseMatrix<_Scalar,Eigen::RowMajor>::InnerIterator it(_delegate.getHessian(),row);
                  it; ++it, ++entry_id )
            {
                values[ entry_id ] = it.value();
            } // ...for col
        } // ...for row
        ret_val = true;

        //    values[0] = 2*lambda[0];
        //    values[1] = 2*lambda[0];

        //        2,1,0.571389
        //        3,2,0.571389
        //        4,1,0.91724
        //        4,2,0.717526
        //        4,3,0.91724
    }

    if ( _delegate.isDebug() )
    {
        std::cout << "[" << __func__ << "]: " << "finished" << std::endl;
        fflush( stdout );
    }

    return ret_val;
} // ...BonminOpt::eval_h()

template <typename _Scalar> void
BonminTMINLP<_Scalar>::finalize_solution( TMINLP::SolverReturn   status
                                        , Index                  n
                                        , Number const         * x
                                        , Number                 obj_value )
{
    std::cout << "Problem status: "  << status    << std::endl;
    std::cout << "Objective value: " << obj_value << std::endl;
    if ( _delegate.isPrintSol() && (x != NULL) )
    {
        std::cout << "Solution:" << std::endl;
        for ( int i = 0 ; i < n; ++i )
        {
            std::cout << "x[" << i << "] = " << x[i];
            if ( i < n-1 ) std::cout << ", ";
        }
        std::cout << std::endl;
    } // ... printSol_

    // save solution
    VectorX sol( n );
    for ( int i = 0 ; i < n; ++i )
    {
        sol(i) = x[i];
    }
    _delegate.setSolution( sol );
} // ... BonminOpt::finalize_solution()

} //...ns GF2

#endif // BONMINOPT_H