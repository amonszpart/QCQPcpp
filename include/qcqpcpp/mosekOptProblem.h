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
        typedef typename ParentType::SG_OBJ_SENSE   SG_OBJ_SENSE;
        typedef typename ParentType::SG_BOUND       SG_BOUND;
        typedef typename ParentType::SG_VAR_TYPE    SG_VAR_TYPE;

        virtual ReturnType update  ( bool verbose = false ) override;
        virtual ReturnType optimize( std::vector<_Scalar> *x_out = NULL
                                     , SG_OBJ_SENSE objecitve_sense = SG_OBJ_SENSE::MINIMIZE ) override;

        MosekOpt( MSKenv_t *env = NULL );
        virtual ~MosekOpt();

        static inline MSKboundkeye     ToMosek( typename ParentType::SG_BOUND bound );
        static inline MSKvariabletypee ToMosek( typename ParentType::SG_VAR_TYPE var_type );

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

