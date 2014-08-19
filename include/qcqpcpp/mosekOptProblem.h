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

//! \brief Specialization of optimisation problem to solve using the Mosek solver
//! \tparam _Scalar Coefficient storage type (floating point).
template <typename _Scalar>
class MosekOpt : public OptProblem<_Scalar>
{
    public:
        typedef          OptProblem<_Scalar>             ParentType;   //!< \brief Parent type typedef.
        typedef typename ParentType::Scalar              Scalar;       //!< \brief Coefficient storage type inherited from parent.
        //typedef MSKrescodee                              ReturnType;   //!< \brief Implementation return type inherited from parent (aka MSKrescodee in this case).
        typedef typename ParentType::ReturnType          ReturnType;   //!< \brief Implementation return type inherited from parent (aka MSKrescodee in this case).
        typedef typename ParentType::OBJ_SENSE           OBJ_SENSE;    //!< \brief Problem optimization sense type inherited from parent.
        typedef typename ParentType::BOUND               BOUND;        //!< \brief Problem bound type inherited from parent.
        typedef typename ParentType::VAR_TYPE            VAR_TYPE;     //!< \brief Problem variable continuity type inherited from parent.
        typedef typename ParentType::SparseMatrix        SparseMatrix; //!< \brief SparseMatrix type inherited from parent.

        //! \brief                  Prepare for optimization. Reads parent content and sets up specialized problem. Has to be called before #optimize().
        virtual ReturnType update  ( bool verbose = false ) override;

        //! \brief                  Run optimization using the current solver. #update() has to be called before.
        //! \param x_out            Returns the output, if requested. The output is also stored in OptProblem::_x
        //! \param objective_sense  Specifies, whether it's a minimization or a maximization problem.
        virtual ReturnType optimize( std::vector<_Scalar> *x_out = NULL
                                     , OBJ_SENSE objecitve_sense = OBJ_SENSE::MINIMIZE ) override;

        //! \brief Overrides infinity (default: std::numeric_limits<Scalar>::max()) to make implementation specific (i.e. here MSK_INFINITY == 1e30).
        virtual Scalar getINF() const override { return MSK_INFINITY; }

        //! \brief      Custom constructor taking optional Mosek environment as initialization.
        //! \pram env   Optional parameter to pass already created Mosek environment. It is advised by Mosek to only have one environment, so if you have to instances, make sure you share the env.
        MosekOpt( MSKenv_t *env = NULL );
        //! \brief Virtual destructor. Frees the Mosek environment, if _ownsEnv is true.
        virtual ~MosekOpt();

        //! \brief      Copy constructor.
        MosekOpt( ParentType const& other );

        static inline MSKboundkeye     getBoundTypeCustom( typename ParentType::BOUND bound );       //!< \brief Converts OptProblem::BOUND to mosek bound type (MSK_BK_*).
        static inline MSKvariabletypee getVarTypeCustom  ( typename ParentType::VAR_TYPE var_type ); //!< \brief Converts OptProblem::VAR_TYPE to mosek var type (MSK_VAR_TYPE_CONT, MSK_VAR_TYPE_INT).

               inline MSKenv_t         getMosekEnv       () const { return _env; } //!< \brief Mosek environment variable getter.
    protected:
        MSKrescodee               _r;                   //!< \brief Mosek error code for operation calls
        MSKenv_t                  _env;                 //!< \brief Mosek environment variable containing all problems. Make sure to only have one of these, and share between instances of MosekOpt.
        bool                      _ownsEnv;             //!< \brief Shows, if this instance is the one that has to free the Mosek environment in the end.
        MSKtask_t                 _task;                //!< \brief One "Problem" in Mosek.

        inline void _storeSolution( double* xx, int n ); //!< \brief #optimize() calls this to store output in OptProblem::_x.
        inline void _init         ( MSKenv_t *env );    //!< \brief Instance init. Should be called from all constructors.
    private:
//        MosekOpt();

}; // ...class MosekOpt

} //... namespace qcqpcpp

#ifndef QCQPCPP_INC_MOSEKOPT_HPP
#   define QCQPCPP_INC_MOSEKOPT_HPP
#   include "qcqpcpp/impl/mosekOptProblem.hpp"
#endif // QCQPCPP_INC_MOSEKOPT_HPP

#endif // QCQPCPP_INC_MOSEKOPT_HPP

