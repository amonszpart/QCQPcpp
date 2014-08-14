#ifndef QCQPCPP_MOSEKOPT_H
#define QCQPCPP_MOSEKOPT_H

#include "Eigen/Sparse"
#include "mosek.h"              /* Include the MOSEK definition file. */
#include "qcqpcpp/optProblem.h"

namespace qcqpcpp
{

static void MSKAPI mosekPrintStr( void *handle, MSKCONST char str[] )
{
    printf( "%s", str );
    fflush( stdout );
} // ... mosekPrintStr()

template <typename _Scalar>
class MosekOpt : public OptProblem<_Scalar,MSKrescodee>
{
    public:
        typedef          OptProblem<_Scalar,MSKrescodee> ParentType;
        typedef typename ParentType::Scalar         Scalar;
        typedef typename ParentType::ReturnType     ReturnType;
        typedef typename ParentType::OBJ_SENSE      OBJ_SENSE;
        typedef typename ParentType::BOUND          BOUND;
        typedef typename ParentType::VAR_TYPE       VAR_TYPE;
        typedef typename ParentType::SparseMatrix   SparseMatrix;

        virtual ReturnType update  ( bool verbose = false ) override;
        virtual ReturnType optimize( std::vector<_Scalar> *x_out = NULL
                                     , OBJ_SENSE objecitve_sense = OBJ_SENSE::MINIMIZE ) override;
        virtual Scalar getINF() const override { return 1.0e30; }

        MosekOpt( MSKenv_t *env = NULL );
        virtual ~MosekOpt();

        static inline MSKboundkeye     ToMosek( typename ParentType::BOUND bound );
        static inline MSKvariabletypee ToMosek( typename ParentType::VAR_TYPE var_type );

        Eigen::Matrix<_Scalar,3,1> checkSolution( std::vector<_Scalar> x, Eigen::Matrix<_Scalar,3,1> weights ) const;

    protected:
        ReturnType                _r;                   //!< \brief Mosek error code for operation calls
        MSKenv_t                  _env;
        bool                      _ownsEnv;
        MSKtask_t                 _task;                //!< \brief One "Problem" in Mosek
    private:
//        MosekOpt();

}; // ...class MosekOpt

} //... namespace qcqpcpp

#ifndef QCQPCPP_INC_MOSEKOPT_HPP
#   define QCQPCPP_INC_MOSEKOPT_HPP
#   include "qcqpcpp/impl/mosekOptProblem.hpp"
#endif // QCQPCPP_INC_MOSEKOPT_HPP

#endif // QCQPCPP_INC_MOSEKOPT_HPP

