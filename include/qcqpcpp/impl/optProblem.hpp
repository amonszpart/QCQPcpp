#ifndef QCQPCPP_SGOPTPROBLEM_HPP
#define QCQPCPP_SGOPTPROBLEM_HPP

#include <exception>

namespace qcqpcpp
{

template <typename _Scalar, typename _ReturnType> int
OptProblem<_Scalar,_ReturnType>::addVariable( BOUND bound_type, Scalar lower_bound, Scalar upper_bound, VAR_TYPE var_type )
{
    // var bound type { MSK_BK_FX = fixed, MSK_BK_FR = free, MSK_BK_LO = blx[j] .. INF, MSK_BK_RA = blx[j] .. bux[j], MSK_BK_UP = -INF .. bux[j] }
    _bkx.push_back( bound_type  );
    // lower bounds
    _blx.push_back( lower_bound );
    // upper bounds
    _bux.push_back( upper_bound );
    // linear objective coeff
    _linObjs  .push_back( Scalar(0)   );
    // variable type
    _type_x.push_back( var_type );

    //return EXIT_SUCCESS;
    return _bkx.size()-1;
} // ...OptProblem::addVariable

//______________________________________________________________________________

template <typename _Scalar, typename _ReturnType> int
OptProblem<_Scalar,_ReturnType>::setLinObjective( int j, Scalar coeff )
{
    if ( static_cast<int>(_linObjs.size()) <= j )
    {
        std::cerr << "[" << __func__ << "]: " << "please add var " << j << " before setting any coeffs" << std::endl;
        return EXIT_FAILURE;
    }

    _linObjs[ j ] = coeff;

    return EXIT_SUCCESS;
} // ...OptProblem::setVarLinCoeff

template <typename _Scalar, typename _ReturnType> int
OptProblem<_Scalar,_ReturnType>::addLinObjective( int j, Scalar coeff )
{
    if ( static_cast<int>(_linObjs.size()) <= j )
    {
        std::cerr << "[" << __func__ << "]: " << "please add var " << j << " before adding any coeffs" << std::endl;
        return EXIT_FAILURE;
    }

    _linObjs[ j ] += coeff;

    return EXIT_SUCCESS;
} // ...OptProblem::addVarLinCoeff

template <typename _Scalar, typename _ReturnType> int
OptProblem<_Scalar,_ReturnType>::addLinObjectives( SparseMatrix const& mx )
{
    int err = EXIT_SUCCESS;

    // mx should be a vector, one dim should be varCount long, the other 1 long.
    if (  ( (mx.outerSize() != getVarCount()) && (mx.innerSize() != getVarCount()) )
       || ( (mx.outerSize() != 1            ) && (mx.innerSize() != 1            ) )
       )
    {
        std::cerr << "[" << __func__ << "]: " << "getVarCount " << getVarCount() << " != " << std::max(mx.outerSize(),mx.innerSize()) << " variables in mx" << std::endl;
        return EXIT_FAILURE;
    }

    // set quadratic objective
    for ( int k = 0; k < mx.outerSize(); ++k )
        for ( typename OptProblem<_Scalar,_ReturnType>::SparseMatrix::InnerIterator it(mx, k); it; ++it )
        {
            //std::cout << "[" << __func__ << "]: " << "adding at " << std::max(it.row(),it.col()) << " = " << it.value() << std::endl;
            err += addLinObjective( std::max(it.row(),it.col()), it.value() );
        }

    return err;
} // ...OptProblem::addVarLinCoeff

template <typename _Scalar, typename _ReturnType> typename OptProblem<_Scalar,_ReturnType>::SparseMatrix
OptProblem<_Scalar,_ReturnType>::getLinObjectivesMatrix() const
{
    SparseMatrix smx( this->getVarCount(), 1 );
    for ( int i = 0; i != this->_linObjs.size(); ++i )
    {
        smx.insert( i, 0 ) = this->_linObjs[i];
    }

    return smx;
} // ...OptProblem::getQuadraticObjectivesMatrix()

template <typename _Scalar, typename _ReturnType> int
OptProblem<_Scalar,_ReturnType>::addQObjective( int i, int j, Scalar coeff )
{
    if ( (i >= this->getVarCount()) || j >= (this->getVarCount()) )
    {
        std::cerr << "[" << __func__ << "]: " << "i " << i << " or j " << j << " > " << this->getVarCount() << ", please call addVariable() first! returning..." << std::endl;
        return EXIT_FAILURE;
    }

    // only lower triangle
    if ( j > i )
        std::swap( i, j );

    _quadObjList.push_back( Eigen::Triplet<Scalar>(i,j,coeff) );

    return EXIT_SUCCESS;
} // ...OptProblem::addQObjectives()

template <typename _Scalar, typename _ReturnType> int
OptProblem<_Scalar,_ReturnType>::addQObjectives( SparseMatrix const& mx )
{
    int err = EXIT_SUCCESS;

    if ( (mx.outerSize() != getVarCount()) || (mx.innerSize() != getVarCount()) )
    {
        std::cerr << "[" << __func__ << "]: " << "getVarCount " << getVarCount() << " != " << mx.outerSize() << " || " << mx.innerSize() << " variables in mx" << std::endl;
        return EXIT_FAILURE;
    }

    // set quadratic objective
    _quadObjList.reserve( _quadObjList.size() + mx.nonZeros() );
    for ( int k = 0; k < mx.outerSize(); ++k )
        for ( typename OptProblem<_Scalar,_ReturnType>::SparseMatrix::InnerIterator it(mx, k); it; ++it )
            err += addQObjective( it.row(), it.col(), it.value() );

    return err;
} // ...OptProblem::addQObjectives()

template <typename _Scalar, typename _ReturnType> int
OptProblem<_Scalar,_ReturnType>::setQObjectives( SparseMatrix const& mx )
{
    _quadObjList.clear();
    return addQObjectives( mx );
} // ...OptProblem::setQObjectives()

template <typename _Scalar, typename _ReturnType> typename OptProblem<_Scalar,_ReturnType>::SparseMatrix
OptProblem<_Scalar,_ReturnType>::getQuadraticObjectivesMatrix() const
{
    SparseMatrix mx( this->getVarCount(), this->getVarCount() );
    mx.setFromTriplets( this->_quadObjList.begin(), this->_quadObjList.end() );

    return mx;
} // ...OptProblem::getQuadraticObjectivesMatrix()

template <typename _Scalar, typename _ReturnType> typename OptProblem<_Scalar,_ReturnType>::SparseMatrix
OptProblem<_Scalar,_ReturnType>::estimateHessianOfLagrangian() const
{
    if ( _quadConstrList.size() )
    {
        std::cerr << "[" << __func__ << "]: " << "Hessian NOT implemented for quadratic constraints yet..." << std::endl;
        throw new std::runtime_error( "[OptProblem::estimateHessianOfLagrangian] Hessian NOT implemented for quadratic constraints yet..." );
    }

    SparseMatrix hessian( getVarCount(), getVarCount() );
    hessian.reserve( _linConstrList.size() );
    for ( int i = 0; i != _linConstrList.size(); ++i )
    {
        // (c * x^2)'' == 2c. Second derivative has a 2 multiplier if it's a squared variable.
        if ( _linConstrList[i].row() != _linConstrList[i].col() )
            hessian.insert( _linConstrList[i].row(), _linConstrList[i].col() ) = _linConstrList[i].value();
        else
            hessian.insert( _linConstrList[i].row(), _linConstrList[i].col() ) = _Scalar(2) * _linConstrList[i].value();
    } // for linConstrList

    return hessian;
} // ...OptProblem::estimateHessianOfLagrangian()

//______________________________________________________________________________

template <typename _Scalar, typename _ReturnType> int
OptProblem<_Scalar,_ReturnType>::addLinConstraint( BOUND bound_type, Scalar lower_bound, Scalar upper_bound, std::vector<Scalar> coeffs  )
{
    // usage:
    //    coeffs[0] * x_0 + coeffs[1] * x_1 ... + coeffs[n] * x_n >= lower_bound              , bound_type == BOUND::GREATER_EQ
    // OR coeffs[0] * x_0 + coeffs[1] * x_1 ... + coeffs[n] * x_n <= upper_bound              , bound_type == BOUND::LESS_EQ
    // OR coeffs[0] * x_0 + coeffs[1] * x_1 ... + coeffs[n] * x_n =  lower_bound = upper_bound, bound_type == BOUND::FIXED
    if ( coeffs.size() != getVarCount() )
    {
        std::cerr << "[" << __func__ << "]: " << "A line in the constraints matrix A has to be varCount " << getVarCount() << " long, not " << coeffs.size() << std::endl;
        return EXIT_FAILURE;
    }

    // prepare
    SparseMatrix row_vector( 1, coeffs.size() );
    // add coeffs from new line
    for ( size_t col = 0; col != coeffs.size(); ++col )
    {
        // add non-zero elements to sparse representation
        if ( coeffs[col] != Scalar(0) )
        {
            row_vector.insert( 0, col ) = coeffs[col];
        }
    }

    // work
    return this->addLinConstraint( bound_type, lower_bound, upper_bound, &row_vector );
}

template <typename _Scalar, typename _ReturnType> int
OptProblem<_Scalar,_ReturnType>::addLinConstraint( BOUND bound_type, Scalar lower_bound, Scalar upper_bound, SparseMatrix *row_vector /* = NULL */ )
{
    // check bounds
    if      ( bound_type == BOUND::EQUAL )
    {
        if ( lower_bound != upper_bound )
        {
            std::cerr << "[" << __func__ << "]: " << "If bound_type is FIXED, upper_bound " << upper_bound << " == " << lower_bound << " lower_bound must hold" << std::endl;
            return EXIT_FAILURE;
        }
    }
    else if ( bound_type == BOUND::GREATER_EQ )
    {
        if ( upper_bound != getINF() )
        {
            std::cerr << "[" << __func__ << "]: " << "If bound_type is LOwer, upper_bound should probably be infinity, and not " << upper_bound << std::endl;
            return EXIT_FAILURE;
        }
    }
    else if ( bound_type == BOUND::LESS_EQ )
    {
        if ( lower_bound != -getINF() )
        {
            std::cerr << "[" << __func__ << "]: " << "If bound_type is UPper, lower_bound should probably be -infinity, and not " << lower_bound << std::endl;
            return EXIT_FAILURE;
        }
    }
    else if ( bound_type == BOUND::RANGE )
    {
        if ( (lower_bound == -getINF()) || (upper_bound == getINF()) )
        {
            std::cerr << "[" << __func__ << "]: " << "If bound_type is RAnge, bounds should not be infinity: " << lower_bound << ", " << upper_bound << std::endl;
            return EXIT_FAILURE;
        }
    }
    else
    {
        std::cerr << "[" << __func__ << "]: " << "bound_type is not defined used for constraints..." << std::endl;
        return EXIT_FAILURE;
    }

    // contraint matrix row
    const int row = _bkc.size();
    // var bound type
    _bkc.push_back( bound_type  );
    // lower bounds
    _blc.push_back( lower_bound );
    // upper bounds
    _buc.push_back( upper_bound );

    // constraint coeffs
    if ( row_vector )
    {
        if ( (row_vector->cols() != this->getVarCount()) || (row_vector->rows() != 1) )
        {
            std::cerr << "[" << __func__ << "]: " << "row_vector->cols()(" << row_vector->cols() << ") != this->getVarCount() (" << this->getVarCount() << ") || row_vector->rows()( " << row_vector->rows() << " != 1, returning" << std::endl;
            return EXIT_FAILURE;
        }

        // add coeffs from new line
        for ( int k = 0; k != row_vector->outerSize(); ++k )
        {
            for ( typename SparseMatrix::InnerIterator it(*row_vector,k); it; ++it )
            {
                _linConstrList.push_back( SparseEntry(row, it.col(), it.value()) );
            }
        }
    }

    return EXIT_SUCCESS;
} // ...OptProblem::addConstraint

template <typename _Scalar, typename _ReturnType> int
OptProblem<_Scalar,_ReturnType>::addLinConstraints( SparseMatrix const& mx )
{
    if ( mx.rows() != this->getConstraintCount() || mx.cols() != this->getVarCount() )
    {
        std::cerr << "[" << __func__ << "]: " << "mx has to be m x n, m == constrCount(" << this->getConstraintCount() << "), n == varCount(" << this->getVarCount() << ")...returning" << std::endl;
        return EXIT_FAILURE;
    }

    this->_linConstrList.reserve( this->_linConstrList.size() + mx.nonZeros() );
    for ( int row = 0; row != mx.outerSize(); ++row )
        for ( typename SparseMatrix::InnerIterator it(mx,row); it; ++it )
        {
            this->_linConstrList.push_back( SparseEntry(it.row(), it.col(), it.value()) );
        } // ... for cols

    return EXIT_SUCCESS;
} // ...OptProblem::setLinConstraints()

template <typename _Scalar, typename _ReturnType> int
OptProblem<_Scalar,_ReturnType>::setLinConstraints( SparseMatrix const& mx )
{
    _linConstrList.clear();
    return this->addLinConstraints( mx );;
} // ...OptProblem::setLinConstraints()

template <typename _Scalar, typename _ReturnType> int
OptProblem<_Scalar,_ReturnType>::addQConstraint( int constr_id, int i, int j, Scalar coeff )
{
    if ( constr_id >= getConstraintCount() )
    {
        std::cerr << "constr_id <= getConstraintCount, please use addLinConstraint to initialize constr_id first! exiting" << std::endl;
        return EXIT_FAILURE;
    }

    if ( _quadConstrList.size() <= constr_id )
        _quadConstrList.resize( constr_id + 1 );

    // only lower triangle
    if ( j > i )
        std::swap( i, j );

    _quadConstrList[constr_id].push_back( SparseEntry(i,j,coeff) );
    std::cout << "[" << __func__ << "]: " << "qconstrlist is now " << _quadConstrList.size() << "( " << _quadConstrList[0].size() << ")" << " long" << std::endl;

    return EXIT_SUCCESS;
} // ...OptProblem::addQConstraint

template <typename _Scalar, typename _ReturnType> typename OptProblem<_Scalar,_ReturnType>::SparseMatrix
OptProblem<_Scalar,_ReturnType>::getLinConstraintsMatrix() const
{
    SparseMatrix mx( this->getConstraintCount(), this->getVarCount() );
    mx.setFromTriplets( this->_linConstrList.begin(), this->_linConstrList.end() );

    return mx;
} // ...OptProblem::getLinConstraintsMatrix()

template <typename _Scalar, typename _ReturnType> typename OptProblem<_Scalar,_ReturnType>::SparseMatrix
OptProblem<_Scalar,_ReturnType>::getQuadraticConstraintsMatrix( int i ) const
{
    SparseMatrix mx( this->getVarCount(), this->getVarCount() );
    mx.setFromTriplets( this->getQuadraticConstraints(i).begin(), this->getQuadraticConstraints(i).end() );

    return mx;
} // ...OptProblem::getLinConstraintsMatrix()

template <typename _Scalar, typename _ReturnType> typename OptProblem<_Scalar,_ReturnType>::SparseMatrix
OptProblem<_Scalar,_ReturnType>::estimateJacobianOfConstraints() const
{
    if ( _quadConstrList.size() )
    {
        std::cerr << "[" << __func__ << "]: " << "Jacobian NOT implemented for quadratic constraints yet..." << std::endl;
        throw new std::runtime_error( "[OptProblem::estimateJacobianOfConstraints] Jacobian NOT implemented for quadratic constraints yet..." );
    }

    SparseMatrix jacobian( getConstraintCount(), getVarCount() );
    jacobian.setFromTriplets( _linConstrList.begin(), _linConstrList.end() );

    return jacobian;
} // ...OptProblem::estimateJacobianOfConstraints()

//______________________________________________________________________________

template <typename _Scalar, typename _ReturnType> int
OptProblem<_Scalar,_ReturnType>::printProblem() const
{
    // linear constraints
    SparseMatrix A = this->getLinConstraintsMatrix();
    std::cout << "A: ";
    for ( int row = 0; row != A.outerSize(); ++row )
    {
        for ( typename SparseMatrix::InnerIterator it(A,row); it; ++it )
        {
            std::cout << "(" << it.col() << "," << it.value() << "), ";
        }
        std::cout << std::endl;
    }

    // linear objectives
    std::cout<<"_linObjs:";for(size_t vi=0;vi!=_linObjs.size();++vi)std::cout<<_linObjs[vi]<<" ";std::cout << "\n";

    return EXIT_SUCCESS;
} // ...OptProblem::printProblem()

} // ...namespace qcqpp

#endif // QCQPCPP_SGOPTPROBLEM_HPP
