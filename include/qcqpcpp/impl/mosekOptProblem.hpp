#ifndef QCQPCPP_MOSEKOPT_HPP
#define QCQPCPP_MOSEKOPT_HPP

#include <iomanip>
#include "mosek.h"              /* Include the MOSEK definition file. */

namespace qcqpcpp
{

template <typename _Scalar>
MosekOpt<_Scalar>::MosekOpt( MSKenv_t *env )
    : _r      ( MSK_RES_OK  )
    , _env    ( env         )
    , _ownsEnv( !env        ) // is true, if env == NULL
    , _task   ( NULL        )
{
    if ( !_env )
    {
        _r = MSK_makeenv( &_env, NULL );
        if ( MSK_RES_OK != _r )
        {
            std::cerr << "[" << __func__ << "]: " << "cannot create Mosek environment, nothing will work from now on..." << std::endl;
            return;
        }
    }
} // ... MosekOpt::MosekOpt()

template <typename _Scalar>
MosekOpt<_Scalar>::~MosekOpt()
{
    if ( _task )
    {
        MSK_deletetask( &_task );
        _task = NULL;
    }

    if ( _env && _ownsEnv )
    {
        MSK_deleteenv( &_env );
        _env = NULL;
    }
}

template <typename _Scalar> typename MosekOpt<_Scalar>::ReturnType
MosekOpt<_Scalar>::
update( bool verbose )
{
    if ( _task != NULL )
    {
        std::cerr << "[" << __func__ << "]: " << "update can only be called once! returning." << std::endl;
        return MSK_RES_ERR_UNKNOWN;
    }

    /* Create the optimization task. */
    if ( MSK_RES_OK == _r )
    {
        _r = MSK_maketask( _env, this->getConstraintCount(), this->getVarCount(), &_task );
        if ( MSK_RES_OK != _r )
            std::cerr << "[" << __func__ << "]: " << "could not create task with " << this->getVarCount() << " vars, and " << this->getConstraintCount() << " constraints" << std::endl;
    }

    // redirect output
    if ( MSK_RES_OK == _r )
    {
        _r = MSK_linkfunctotaskstream( _task, MSK_STREAM_LOG, NULL, mosekPrintStr );
        if ( MSK_RES_OK != _r )
            std::cerr << "[" << __func__ << "]: " << "could not create rewire output to mosekPrintStr(), continuing though..." << std::endl;
    }

    // Append _numCon empty constraints. The constraints will initially have no bounds.
    if ( MSK_RES_OK == _r )
    {
        if ( verbose ) std::cout << "my: MSK_appendcons(_task,"<< this->getConstraintCount() <<");" << std::endl;
        _r = MSK_appendcons( _task, this->getConstraintCount() );
        if ( MSK_RES_OK != _r )
            std::cerr << "[" << __func__ << "]: " << "could not append " << this->getConstraintCount() << " constraints" << std::endl;
    }

    // Append _numVar variables. The variables will initially be fixed at zero (x=0).
    if ( MSK_RES_OK == _r )
    {
        if ( verbose ) std::cout << "my: MSK_appendvars(_task," << this->getVarCount() <<");" << std::endl;
        _r = MSK_appendvars( _task, this->getVarCount() );
        if ( MSK_RES_OK != _r )
            std::cerr << "[" << __func__ << "]: " << "could not append " << this->getVarCount() << " variables" << std::endl;
    }

    // Optionally add a constant term to the objective.
    if ( MSK_RES_OK == _r )
    {
        if ( verbose ) std::cout << "my: MSK_putcfix(_task," << this->getObjectiveBias() << ");" << std::endl;
        _r = MSK_putcfix( _task, this->getObjectiveBias() );
        if ( MSK_RES_OK != _r )
            std::cerr << "[" << __func__ << "]: " << "could not add constant " << this->getObjectiveBias() << " to objective function" << std::endl;
    }

    // set Variables
    for ( size_t j = 0; (j < this->getVarCount()) && (MSK_RES_OK == _r); ++j )
    {
        // set Variable j's Bounds // blx[j] <= x_j <= bux[j]
        if ( MSK_RES_OK == _r )
        {
            _r = MSK_putvarbound( _task,
                                  j,                                                     /* Index of variable.*/
                                  MosekOpt<Scalar>::ToMosek( this->getVarBoundType(j) ), /* Bound key.*/
                                  this->getVarLowerBound(j),                             /* Numerical value of lower bound.*/
                                  this->getVarUpperBound(j) );                           /* Numerical value of upper bound.*/

            if ( verbose ) std::cout << "my: MSK_putvarbound(_task," << j << "," << this->getVarBoundType(j) << "," << this->getVarLowerBound(j) << "," << this->getVarUpperBound(j) << ");" << std::endl;
        }

        // set Variable j's Type
        if ( MSK_RES_OK == _r )
        {
            _r = MSK_putvartype( _task, j, MosekOpt<Scalar>::ToMosek(this->getVarType(j)) );
        }

        // set Variable j's linear coefficient in the objective function
        if ( MSK_RES_OK == _r )
        {
            if ( verbose ) std::cout << "my: putcj(_task," << j << "," << this->getLinObjectives()[j] << ")" << std::endl;
            _r = MSK_putcj( _task, j, this->getLinObjectives()[j] );
        }
    }

    // set Quadratic Objectives
    if ( MSK_RES_OK == _r )
    {
        const int numNonZeros = this->getQuadraticObjectives().size();
        MSKint32t *qsubi = new MSKint32t[numNonZeros],
                  *qsubj = new MSKint32t[numNonZeros];
        double    *qval  = new double[numNonZeros];

        for ( size_t qi = 0; qi != this->getQuadraticObjectives().size(); ++qi )
        {
            qsubi[qi] = this->getQuadraticObjectives()[qi].row();
            qsubj[qi] = this->getQuadraticObjectives()[qi].col();
            qval [qi] = this->getQuadraticObjectives()[qi].value();
        }

        if ( verbose ) std::cout<<"my: putqobj( _task, " << numNonZeros << ",\n";
        for ( size_t vi = 0; vi != numNonZeros; ++vi )
        {
            if ( verbose ) std::cout << qsubi[vi] << "," << qsubj[vi] << ", " << qval[vi] << std::endl;
        }
        if ( verbose ) std::cout << ");" << std::endl;

        _r = MSK_putqobj( _task, numNonZeros, qsubi, qsubj, qval );

        if ( qsubi ) { delete[] qsubi; qsubi = NULL; }
        if ( qsubj ) { delete[] qsubj; qsubj = NULL; }
        if ( qval  ) { delete[] qval ; qval  = NULL; }

        if ( MSK_RES_OK != _r )
            std::cerr << "[" << __func__ << "]: " << "Setting Quadratic Objectives caused error code " << (int)_r << std::endl;
    } // ...Quadratic objective

    // set Linear Constraints
    {
        typename ParentType::SparseMatrix A( this->getLinConstraintsMatrix() );
//        ( this->getConstraintCount(), this->getVarCount() );
//        A.setFromTriplets( this->getLinConstraints().begin(), this->getLinConstraints().end() );
        std::vector<Scalar>         aval;                //!< \brief Linear constraints coeff matrix (sparse)
        std::vector<int>            asub;                //!< \brief Linear constraints coeff matrix indices
        std::vector<int>            aptrb, aptre;
        for ( int row = 0; (row < A.outerSize()) && (MSK_RES_OK == _r); ++row )
        {
            // set Constraint Bounds for row
            if ( MSK_RES_OK == _r )
            {
                if ( verbose ) std::cout << "my: MSK_putconbound( _task, " << row << ", "
                                         << MosekOpt<Scalar>::ToMosek( this->getConstraintBoundType(row) ) << ", "
                                         << this->getConstraintLowerBound( row ) << ", "
                                         << this->getConstraintUpperBound( row ) << ")"
                                         << std::endl; fflush( stdout );

                _r = MSK_putconbound( _task,
                                      row,                                                          /* Index of constraint.*/
                                      MosekOpt<Scalar>::ToMosek(this->getConstraintBoundType(row)), /* Bound key.*/
                                      this->getConstraintLowerBound(row),                           /* Numerical value of lower bound.*/
                                      this->getConstraintUpperBound(row) );                         /* Numerical value of upper bound.*/
            }

            // set Linear Constraint row
            if ( MSK_RES_OK == _r )
            {
                // new line starts at index == current size
                aptrb.push_back( aval.size() );
                // add coeffs from new line
                for ( typename ParentType::SparseMatrix::InnerIterator it(A,row); it; ++it )
                {
                    if ( row != it.row() ) std::cerr << "[" << __func__ << "]: " << "this shouldn't happen" << std::endl;
                    // coeff value
                    aval.push_back( it.value() );  // TODO: A should be a matrix, not a vector...
                    // coeff subscript
                    asub.push_back( it.col() );
                }
                // new line ends at index == new size
                aptre.push_back( aval.size() );

                if ( verbose ) {
                    std::cout << "my: MSK_putarow( _task, "
                              << row << ", "
                              << aptre[row] - aptrb[row] << ", "
                              << *(asub.data() + aptrb[row]) << ", "
                              << *(aval.data() + aptrb[row]) << ");"
                              << std::endl; fflush( stdout );
                }

                _r = MSK_putarow( _task,
                                  row,                 /* Row index.*/
                                  aptre[row] - aptrb[row], /* Number of non-zeros in row i.*/
                                  asub.data() + aptrb[row],     /* Pointer to column indexes of row i.*/
                                  aval.data() + aptrb[row]);    /* Pointer to values of row i.*/
            }
        } // ... for A.rows

        // report error
        if ( MSK_RES_OK != _r )
            std::cerr << "[" << __func__ << "]: " << "Setting Lin constraints caused error code " << (int)_r << std::endl;
    } // ...set Linear Constraints

    // set Quadratic constraints
    if ( verbose ) std::cout << "[" << __func__ << "]: " << "adding q constraints" << std::endl;
    for ( size_t constr_id = 0; (constr_id != this->getQuadraticConstraints().size()) && (MSK_RES_OK == _r); ++constr_id )
    {
        const int numNonZeros = this->getQuadraticConstraints(constr_id).size();

        MSKint32t *qsubi = new MSKint32t[numNonZeros],
                  *qsubj = new MSKint32t[numNonZeros];
        double    *qval  = new double[numNonZeros];

        for ( size_t qi = 0; qi != this->getQuadraticConstraints(constr_id).size(); ++qi )
        {
            qsubi[qi] = this->getQuadraticConstraints(constr_id)[qi].row();
            qsubj[qi] = this->getQuadraticConstraints(constr_id)[qi].col();
            qval [qi] = this->getQuadraticConstraints(constr_id)[qi].value();
        }

        if ( verbose ) std::cout<<"my: MSK_putqonk( _task, " << constr_id << ", " << numNonZeros << ",\n";
        for(size_t vi=0;vi!=numNonZeros;++vi)
        {
            if ( verbose ) std::cout << qsubi[vi] << "," << qsubj[vi] << ", " << qval[vi] << std::endl;
        }
        if ( verbose ) std::cout << "); " << std::endl;

        _r = MSK_putqconk(_task,
                          constr_id,
                          numNonZeros,
                          qsubi,
                          qsubj,
                          qval);

        if ( qsubi ) { delete[] qsubi; qsubi = NULL; }
        if ( qsubj ) { delete[] qsubj; qsubj = NULL; }
        if ( qval  ) { delete[] qval ; qval  = NULL; }

        if ( MSK_RES_OK != _r )
            std::cerr << "[" << __func__ << "]: " << "Setting Quad constraints caused error code " << (int)_r << std::endl;
    } // ...set Quadratic Constraints

    // save to file
    {
        if ( _r == MSK_RES_OK )
        {
            _r = MSK_putintparam( _task, MSK_IPAR_WRITE_DATA_FORMAT, MSK_DATA_FORMAT_LP );
            if ( _r == MSK_RES_OK )
            {
                _r = MSK_writedata( _task, "mosek.lp" );
                if ( _r == MSK_RES_OK )
                {
                    std::cerr << "[" << __func__ << "]: " << "Writedata did not work" << std::endl;
                }
            }
        }
    }

    // return error code
    return _r;
} // ...MosekOpt::update()

template <typename _Scalar> typename MosekOpt<_Scalar>::ReturnType
MosekOpt<_Scalar>::optimize( std::vector<_Scalar> *x_out, OBJ_SENSE objective_sense )
{
    // cache problem size
    const int numvar = this->getVarCount();

    // determine problem type
    MSKobjsense_enum objsense = (objective_sense == OBJ_SENSE::MINIMIZE) ? MSK_OBJECTIVE_SENSE_MINIMIZE
                                                                         : MSK_OBJECTIVE_SENSE_MAXIMIZE;
    if ( MSK_RES_OK == _r )
        _r = MSK_putobjsense( _task, objsense );

    if ( MSK_RES_OK == _r  )
    {
        // set termination sensitivity
        MSKrescodee trmcode;
        MSK_putdouparam( _task, MSK_DPAR_MIO_TOL_REL_GAP, 1e-10f );

        if (_r == MSK_RES_OK)
        {
            //_r = MSK_putintparam(_task, MSK_IPAR_OPTIMIZER, MSK_OPTIMIZER_MIXED_INT_CONIC );
            if ( _r != MSK_RES_OK )
            {
                std::cerr << "[" << __func__ << "]: " << "setting MSK_OPTIMIZER_MIXED_INT_CONIC did not work!" << std::endl;
            }
        }

        if ( _r == MSK_RES_OK )
        {
            _r = MSK_putintparam( _task, MSK_IPAR_MIO_PRESOLVE_USE, MSK_OFF );
            if ( _r != MSK_RES_OK )
            {
                std::cerr << "[" << __func__ << "]: " << "setting MSK_IPAR_MIO_PRESOLVE_USE did not work!" << std::endl;
            }
        }

        if ( _r == MSK_RES_OK )
        {
            _r = MSK_putintparam( _task, MSK_IPAR_MIO_HEURISTIC_LEVEL, 5 );
            if ( _r != MSK_RES_OK )
            {
                std::cerr << "[" << __func__ << "]: " << "setting MSK_IPAR_MIO_HEURISTIC_LEVEL did not work!" << std::endl;
            }
        }

        // Run optimizer
        _r = MSK_optimizetrm( _task, &trmcode );

        // Print a summary containing information about the solution for debugging purposes.
        MSK_solutionsummary( _task, MSK_STREAM_LOG );

        // save solution
        double *xx = (double*) calloc(numvar,sizeof(double));
        if ( _r == MSK_RES_OK )
        {
            MSKsolstae solsta;

            if ( _r == MSK_RES_OK )
            {
                _r = MSK_getsolsta( _task, MSK_SOL_ITR, &solsta );
                if ( _r != MSK_RES_OK )
                {
                    _r = MSK_getsolsta( _task, MSK_SOL_ITG, &solsta );
                }
                if ( _r != MSK_RES_OK )
                {
                    std::cerr << "[" << __func__ << "]: " << "neithter MSK_SOL_ITR, nor MSK_SOL_ITR worked" << std::endl;
                }
            }

            switch ( solsta )
            {
                case MSK_SOL_STA_OPTIMAL:
                case MSK_SOL_STA_NEAR_OPTIMAL:
                {
                    if ( xx )
                    {
                        MSK_getxx(_task,
                                  MSK_SOL_ITR,    /* Request the basic solution. */
                                  xx);

                        printf("Optimal primal solution\n");
                        if ( x_out ) { x_out->clear(); x_out->reserve(numvar); }
                        for( int j=0; j<numvar; ++j)
                        {
                            if ( x_out )
                                x_out->push_back( xx[j] );
                            //printf("x[%d]: %e\n",j,xx[j]);
                        }
                    }
                    else
                    {
                        _r = MSK_RES_ERR_SPACE;
                    }
                    break;
                }

                case MSK_SOL_STA_DUAL_INFEAS_CER:
                case MSK_SOL_STA_PRIM_INFEAS_CER:
                case MSK_SOL_STA_NEAR_DUAL_INFEAS_CER:
                case MSK_SOL_STA_NEAR_PRIM_INFEAS_CER:
                    printf("Primal or dual infeasibility certificate found.\n");
                    break;
                case MSK_SOL_STA_UNKNOWN:
                {
                    MSKprostae prosta;
                    MSK_getprosta(_task,MSK_SOL_ITG,&prosta);
                    switch (prosta)
                    {
                        case MSK_PRO_STA_PRIM_INFEAS_OR_UNBOUNDED:
                            printf("Problem status Infeasible or unbounded\n");
                            break;
                        case MSK_PRO_STA_PRIM_INFEAS:
                            printf("Problem status Infeasible.\n");
                            break;
                        case MSK_PRO_STA_UNKNOWN:
                            printf("Problem status unknown.\n");
                            break;
                        default:
                            printf("Other problem status.");
                            break;
                    }
                    char symname[MSK_MAX_STR_LEN];
                    char desc[MSK_MAX_STR_LEN];

                    /* If the solutions status is unknown, print the termination code
               indicating why the optimizer terminated prematurely. */

                    MSK_getcodedesc(trmcode,
                                    symname,
                                    desc);

                    printf("The solutuion status is unknown.\n");
                    printf("The optimizer terminitated with code: %s\n",symname);
                    break;
                }
                    // ITG
                    //asdf todo: consolidate this last part:
                case MSK_SOL_STA_INTEGER_OPTIMAL:
                case MSK_SOL_STA_NEAR_INTEGER_OPTIMAL :
                    MSK_getxx(_task,
                              MSK_SOL_ITG,    /* Request the integer solution. */
                              xx);

                    printf("Optimal integer solution.\n");
                    if ( x_out ) { x_out->clear(); x_out->reserve(numvar); }
                    for( int j=0; j<numvar; ++j)
                    {
                        if ( x_out )
                            x_out->push_back( xx[j] );
                        //printf("x[%d]: %e\n",j,xx[j]);
                    }
                    break;

                case MSK_SOL_STA_PRIM_FEAS:
                    /* A feasible but not necessarily optimal solution was located. */
                    MSK_getxx(_task,MSK_SOL_ITG,xx);

                    printf("Feasible solution.\n");
                    if ( x_out ) { x_out->clear(); x_out->reserve(numvar); }
                    for( int j=0; j<numvar; ++j)
                    {
                        if ( x_out )
                            x_out->push_back( xx[j] );
                        //printf("x[%d]: %e\n",j,xx[j]);
                    }
                    break;

                default:
                    std::cerr << "[" << __func__ << "]: " << "unknown code " << (int)solsta << std::endl;
                    break;
            }

            if ( xx ) { free(xx); xx = NULL; }
        }
    }

    if ( MSK_RES_OK != _r )
    {
        /* In case of an error print error code and description. */
        char symname[MSK_MAX_STR_LEN];
        char desc[MSK_MAX_STR_LEN];

        printf("An error occurred while optimizing.\n");
        MSK_getcodedesc( _r,
                         symname,
                         desc);
        printf("Error %s - '%s'\n",symname,desc);
    }

    return _r;
} // ...MosekOpt::optimize()

template <typename _Scalar> Eigen::Matrix<_Scalar,3,1>
MosekOpt<_Scalar>::checkSolution( std::vector<_Scalar> x, Eigen::Matrix<_Scalar,3,1> weights ) const // TODO: probably move to optproblem.h
{
    Eigen::Matrix<_Scalar,3,1> energy; energy.setZero();

    SparseMatrix complexity( x.size(), 1 );
    for ( int row = 0; row != x.size(); ++row )
        complexity.insert( row, 0 ) = weights(2);

    // X
    SparseMatrix mx( x.size(), 1 );
    for ( size_t i = 0; i != x.size(); ++i )
        mx.insert( i, 0 ) = x[i];

    SparseMatrix linObj = this->getLinObjectivesMatrix();
    SparseMatrix data   = linObj - complexity;

    // qo
    SparseMatrix e02 = mx.transpose() * linObj;
    std::cout << "[" << __func__ << "]: " << "qo * x = " << e02.coeffRef(0,0) << std::endl; fflush(stdout);

    // datacost
    SparseMatrix e0 = (mx.transpose() * data);
    energy(0) = e0.coeffRef(0,0);
    //std::cout << "[" << __func__ << "]: " << "data: " << energy(0) << std::endl; fflush(stdout);

    // Qo
    SparseMatrix e1 = mx.transpose() * this->getQuadraticObjectivesMatrix() * mx;
    energy(1) = e1.coeffRef(0,0);
    //std::cout << "[" << __func__ << "]: " << "x' * Qo * x = pw = " << energy(1) << std::endl; fflush(stdout);

    // complexity
    SparseMatrix e2 = mx.transpose() * complexity;
    energy(2) = e2.coeffRef(0,0);
    //std::cout << "[" << __func__ << "]: " << "complx = " << energy(2) << std::endl; fflush(stdout);
    std::cout << "[" << __func__ << "]: " << std::setprecision(9) << energy(0) << " + " << energy(1) << " + " << energy(2) << " = " << energy.sum() << std::endl;

    return energy;
}

template <typename _Scalar> MSKboundkeye
MosekOpt<_Scalar>::ToMosek( typename MosekOpt<_Scalar>::BOUND bound )
{
    switch( bound )
    {
        case BOUND::GREATER_EQ:
            return MSK_BK_LO;
            break;
        case BOUND::LESS_EQ:
            return MSK_BK_UP;
            break;
        case BOUND::EQUAL:
            return MSK_BK_FX;
            break;
        case BOUND::FREE:
            return MSK_BK_FR;
            break;
        case BOUND::RANGE:
            return MSK_BK_RA;
            break;
        default:
            std::cerr << "[" << __func__ << "]: " << "Invalid bound type, cannot convert" << std::endl;
            return MSK_BK_FR;
            break;
    }
} // ...MosekOpt::ToMosek( bound )

template <typename _Scalar> MSKvariabletypee
MosekOpt<_Scalar>::ToMosek( MosekOpt<_Scalar>::VAR_TYPE var_type )
{
    switch( var_type  )
    {
        case VAR_TYPE::CONTINUOUS:
            return MSK_VAR_TYPE_CONT;
            break;
        case VAR_TYPE::INTEGER:
            return MSK_VAR_TYPE_INT;
            break;
        default:
            std::cerr << "[" << __func__ << "]: " << "Invalid var type, cannot convert" << std::endl;
            return MSK_VAR_TYPE_CONT;
            break;
    }
} // ...MosekOpt::ToMosek( var_type )

} //...namespace qcqpcpp

#endif // QCQPCPP_MOSEKOPT_HPP
