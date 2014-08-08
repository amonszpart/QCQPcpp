#ifndef QCQPCPP_SGOPTPROBLEM_HPP
#define QCQPCPP_SGOPTPROBLEM_HPP

namespace qcqpcpp
{

template <typename _Scalar, typename _ReturnType> int
OptProblem<_Scalar,_ReturnType>::addVariable( SG_BOUND bound_type, Scalar lower_bound, Scalar upper_bound, SG_VAR_TYPE var_type )
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

    return EXIT_SUCCESS;
} // ...MosekOpt::addVariable

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
} // ...MosekOpt::setVarLinCoeff

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
} // ...MosekOpt::addVarLinCoeff

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
} // ...MosekOpt::addVarLinCoeff

template <typename _Scalar, typename _ReturnType> int
OptProblem<_Scalar,_ReturnType>::addQObjective( int i, int j, Scalar coeff )
{
    if ( (i >= this->getVarCount()) || j >= (this->getVarCount()) )
    {
        std::cerr << "[" << __func__ << "]: " << "i " << i << " or j " << j << " > " << this->getVarCount() << ", please call addVariable() first! returning..." << std::endl;
        return EXIT_FAILURE;
    }

    _quadObjList.push_back( Eigen::Triplet<Scalar>(i,j,coeff) );

    return EXIT_SUCCESS;
} // ...MosekOpt::addQObjectives()

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
} // ...MosekOpt::addQObjectives()

template <typename _Scalar, typename _ReturnType> int
OptProblem<_Scalar,_ReturnType>::setQObjectives( SparseMatrix const& mx )
{
    _quadObjList.clear();
    return addQObjectives( mx );
} // ...MosekOpt::setQObjectives()

template <typename _Scalar, typename _ReturnType> int
OptProblem<_Scalar,_ReturnType>::addLinConstraint( std::vector<Scalar> coeffs, SG_BOUND bound_type, Scalar lower_bound, Scalar upper_bound )
{
    // usage:
    //    coeffs[0] * x_0 + coeffs[1] * x_1 ... + coeffs[n] * x_n >= lower_bound              , bound_type == MSK_BK_LO
    // OR coeffs[0] * x_0 + coeffs[1] * x_1 ... + coeffs[n] * x_n <= upper_bound              , bound_type == MSK_BK_UP
    // OR coeffs[0] * x_0 + coeffs[1] * x_1 ... + coeffs[n] * x_n =  lower_bound = upper_bound, bound_type == MSK_BK_FX
    if ( coeffs.size() != getVarCount() )
    {
        std::cerr << "[" << __func__ << "]: " << "A line in the constraints matrix A has to be varCount " << getVarCount() << " long, not " << coeffs.size() << std::endl;
        return EXIT_FAILURE;
    }

    if      ( bound_type == SG_BOUND::EQUAL )
    {
        if ( lower_bound != upper_bound )
        {
            std::cerr << "[" << __func__ << "]: " << "If bound_type is FIXED, upper_bound " << upper_bound << " == " << lower_bound << " lower_bound must hold" << std::endl;
            return EXIT_FAILURE;
        }
    }
    else if ( bound_type == SG_BOUND::GREATER_EQ )
    {
        if ( upper_bound != MSK_INFINITY )
        {
            std::cerr << "[" << __func__ << "]: " << "If bound_type is LOwer, upper_bound should probably be infinity, and not " << upper_bound << std::endl;
            return EXIT_FAILURE;
        }
    }
    else if ( bound_type == SG_BOUND::LESS_EQ )
    {
        if ( lower_bound != -MSK_INFINITY )
        {
            std::cerr << "[" << __func__ << "]: " << "If bound_type is UPper, lower_bound should probably be -infinity, and not " << lower_bound << std::endl;
            return EXIT_FAILURE;
        }
    }
    else if ( bound_type == SG_BOUND::RANGE )
    {
        if ( (lower_bound == -MSK_INFINITY) || (upper_bound == MSK_INFINITY) )
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
    {
        // add coeffs from new line
        for ( size_t i = 0; i != coeffs.size(); ++i )
        {
            // add non-zero elements to sparse representation
            if ( coeffs[i] != Scalar(0) )
            {
                _linConstrList.push_back( SparseEntry(row, i, coeffs[i]) );
            }
        }
    }

    return EXIT_SUCCESS;
} // ...MosekOpt::addConstraint

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

    _quadConstrList[constr_id].push_back( SparseEntry(i,j,coeff) );
    std::cout << "[" << __func__ << "]: " << "qconstrlist is now " << _quadConstrList.size() << "( " << _quadConstrList[0].size() << ")" << " long" << std::endl;

    return EXIT_SUCCESS;
} // ...MosekOpt::addQConstraint

template <typename _Scalar, typename _ReturnType> int
OptProblem<_Scalar,_ReturnType>::printProblem() const
{
    std::cout<<"LinObjective:";for(size_t vi=0;vi!=_linObjs.size();++vi)std::cout<<_linObjs[vi]<<" ";std::cout << "\n";

//    std::cout<<"_aptrb:";for(size_t vi=0;vi!=_aptrb.size();++vi)std::cout<<_aptrb[vi]<<" ";std::cout << "\n";
//    std::cout<<"_aptre:";for(size_t vi=0;vi!=_aptre.size();++vi)std::cout<<_aptre[vi]<<" ";std::cout << "\n";
//    std::cout<<"_asub:";for(size_t vi=0;vi!=_asub.size();++vi)std::cout<<_asub[vi]<<" ";std::cout << "\n";
//    std::cout<<"_aval:";for(size_t vi=0;vi!=_aval.size();++vi)std::cout<<_aval[vi]<<" ";std::cout << "\n";
    return EXIT_SUCCESS;
} // ...MosekOpt::printProblem()

//template <typename _Scalar, typename _ReturnType> size_t
//SGOpt<_Scalar,_ReturnType>::_countNonZeros( std::vector<Scalar> coeffs )
//{
//    int cnt = 0;
//    for ( size_t i = 0; i != coeffs.size(); ++i )
//        if ( coeffs[i] != Scalar(0) )
//            ++cnt;

//    return cnt;
//} // ... Mosek::_countNonZeros

} // ...namespace qcqpp

#endif // QCQPCPP_SGOPTPROBLEM_HPP
