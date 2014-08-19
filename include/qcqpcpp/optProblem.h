#ifndef QCQPCPP_INC_SGOPTPROBLEM_H
#define QCQPCPP_INC_SGOPTPROBLEM_H

#include <iostream> // cout, cerr, endl
#include "Eigen/Sparse"

namespace qcqpcpp
{

template <typename _Scalar>
class OptProblem
{
    public:
        typedef _Scalar                                     Scalar;         //!< \brief Some implementations require this to be double (i.e. Mosek does).
        //typedef _ReturnType                                 ReturnType;     //!< \brief General returntype of the library used, int or MSKrescodee, etc..
        typedef int                                         ReturnType;     //!< \brief Error code storage type.
        typedef Eigen::Matrix<_Scalar,-1,1>                 VectorX;        //!< \brief RowVector of unkown dimensions.
        typedef Eigen::Triplet<Scalar>                      SparseEntry;    //!< \brief This type is an element of a list of entries before construction of a SparseMatrix.
        typedef std::vector<SparseEntry>                    SparseEntries;  //!< \brief A list of entries before construction. SparseMatrices are created in a lazy fashion.
        typedef Eigen::SparseVector<Scalar,Eigen::RowMajor> SparseVector;   //!< \brief
        typedef Eigen::SparseMatrix<Scalar,Eigen::RowMajor> SparseMatrix;   //!< \brief To store quadratic objective (Q_o) and quadratic constraints (Qi). SparseMatrices are created in a lazy fashion.

        enum OBJ_SENSE   { MINIMIZE       //!< \brief Minimize objective function in optimize()
                         , MAXIMIZE       //!< \brief Maximize objective function in optimize()
                         };
        enum VAR_TYPE    { CONTINUOUS = 0 //!< \brief Default
                         , INTEGER    = 1 //!< \brief Mixed Integer Programming
                         , BINARY     = 2
                         };
        enum LINEARITY   { LINEAR     = 0
                         , NON_LINEAR = 1 //!< \brief NOT USED, here for completeness sake
                         };
        enum BOUND       { GREATER_EQ = 0 //!< \brief Means: lower_bound < ... < +INF
                         , LESS_EQ    = 1 //!< \brief Means: -INF        < ... < upper_bound
                         , EQUAL      = 2 //!< \brief Means: Const       = ... = Const
                         , FREE       = 3 //!< \brief Means: -INF        = ... = +INF
                         , RANGE      = 4 //!< \brief Means: lower_bound = ... = upper_bound
                         };

        //! \brief Internal structure to store a variable
        struct Variable
        {
                BOUND          _bkx;                 //!< \brief Variable bound key
                Scalar         _blx;                 //!< \brief Variable numerical value of lower bound
                Scalar         _bux;                 //!< \brief Variable numerical value of upper bound
                VAR_TYPE       _type_x;              //!< \brief Variable type
                LINEARITY      _lin_x;               //!< \brief Vairable linearity.
        };

        //static constexpr Scalar INF = 1.0e30; //!< \brief Used, when a bound direction is "unbounded"
        virtual Scalar getINF() const { return std::numeric_limits<Scalar>::max(); }

        //! \brief Make sure, to call before optimize is called, but AFTER the problem details are added using the setters below.
        virtual ReturnType                       update                 ( bool verbose = false )                                          { std::cerr << "[" << __func__ << "]: " << "Please override with library specific set-up logic. See MosekOpt.h for concept." << std::endl; return ReturnType(EXIT_FAILURE); }
        //! \brief Call to run optimization, AFTER the problem details were added, and <i>update()</i> was called.
        virtual ReturnType                       optimize               ( std::vector<Scalar> *x_out           = NULL,
                                                                          OBJ_SENSE            objecitve_sense = OBJ_SENSE::MINIMIZE ) { std::cerr << "[" << __func__ << "]: " << "Please override with library specific set-up logic. See MosekOpt.h for concept." << std::endl; return ReturnType(EXIT_FAILURE); }

        //! \brief Constructor unused for the moment. Start with <i>addVariable()</i>.
                                                  OptProblem             () : _updated(false) {}
        //! \brief Destructor unused for the moment. Declared virtual for inheritence.
        virtual                                  ~OptProblem            () { std::cout << "[" << __func__ << "][INFO]: " << "Empty destructor" << std::endl; }

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //// Variables ////////////////////////////////////////////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////

        //! \brief addVariable  Start setting up the problem by defining the X vector.
        //! \param bound_type   A variable can be FREE (unbounded), EQUAL (fixed), lower bounded (GREATER_EQ) or upper bounded (LESS_EQ), or bounded (RANGE).
        //! \param lower_bound  Specify lower bound (set to -INF), if not lower bounded (FREE, LESS_EQ).
        //! \param upper_bound  Specify upper bound (set to +INF), if not upper bounded (FREE, GREATER_EQ).
        //! \param var_type     Continuous (default) or integer (for mixed integer programming).
        //! \return             0-indexed id of the just added variable.
        inline int                               addVariable            ( BOUND     bound_type     = BOUND::FREE
                                                                        , Scalar    lower_bound    = -getINF()
                                                                        , Scalar    upper_bound    = +getINF()
                                                                        , VAR_TYPE  var_type       = VAR_TYPE::CONTINUOUS
                                                                        , LINEARITY var_lin        = LINEARITY::LINEAR );
        inline int                               addVariable            ( Variable const& var ) { return this->addVariable( var._bkx, var._blx, var._bux, var._type_x, var._lin_x ); }

        inline VAR_TYPE                   const  getVarType             ( int j )         const { return _type_x[j]; }
        inline BOUND                      const  getVarBoundType        ( int j )         const { return _bkx[j]; }
        inline Scalar                     const  getVarLowerBound       ( int j )         const { return _blx[j]; }
        inline Scalar                     const  getVarUpperBound       ( int j )         const { return _bux[j]; }
        inline LINEARITY                  const  getVarLinearity        ( int j )         const { return _lin_x[j]; }
        inline size_t                            getVarCount            ()                const { return _bkx.size(); }     //!< \brief Returns the number of variables currently in the system

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //// Objective function ///////////////////////////////////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////

        /// Formulation: minimize/maximize X' * 1/2 * Q_o * X + q_o * X + c

        // Const objective function (c)
        inline int                               setObjectiveBias               ( Scalar cfix ) { _cfix = cfix; return EXIT_SUCCESS; } //!< \brief Is always added to the objective function as bias.
        inline Scalar                            getObjectiveBias               ()        const { return _cfix; }                      //!< \brief Read-only getter for fixed objective function bias.
        // Linear objective function (q_o)
        inline int                               setLinObjective                ( int j, Scalar coeff );                 //!< \brief Sets the j-th linear coefficient to coeff, regardless previous value.
        inline int                               setLinObjectives               ( Eigen::Matrix<_Scalar,-1,1> const& obj ) { for ( int i = 0; i != obj.rows(); ++i ) setLinObjective(i, obj(i) ); return EXIT_SUCCESS; }
        inline int                               addLinObjective                ( int j, Scalar coeff );                 //!< \brief Adds coeff to the j-th linear coefficient.
        inline int                               addLinObjectives               ( SparseMatrix const& mx );              //!< \brief Assumes vector input, so ((mx.cols == 1) || (mx.rows == 1)).
        inline std::vector<Scalar>        const& getLinObjectives               ()        const { return _linObjs; }     //!< \brief Read-only getter for linear objective vector.
        inline SparseMatrix                      getLinObjectivesMatrix         ()        const;                         //!< \brief Returns a 1 column SparseMatrix temporary.
        // Quadratic Objecitve function (Q_o)
        inline int                               addQObjective                  ( int i, int j, Scalar coeff );          //!< \brief Adds input to quadratic objective matrix, duplicate entries get summed in a lazy fashion (when update() is called).
        inline int                               addQObjectives                 ( SparseMatrix const& mx );              //!< \brief Adds input matrix to quadratic objective matrix, duplicate entries get summed in a lazy fashion (when update() is called).
        inline int                               setQObjectives                 ( SparseMatrix const& mx );              //!< \brief Clears current objectives, and sets input as quadratic objective matrix.
        inline SparseEntries              const& getQuadraticObjectives         ()        const { return _quadObjList; } //!< \brief Read-only getter for quadratic objective matrix's unordered entry set.
        inline SparseMatrix                      getQuadraticObjectivesMatrix   ()        const;                         //!< \brief Returns all Qo entries assembled to a Sparsematrix, temporarily.
        // Hessian of Lagrangian
        //! \brief  estimateHessianOfLagrangian  Estimates Hessian of Lagrangian from Quadratic objective matrix.
        //! \return SparseMatrix                 Lower triangle only, assumed to be symmetric. (TODO: add quadratic constraints).
        inline SparseMatrix                      estimateHessianOfLagrangian    ()        const;

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //// Constraints //////////////////////////////////////////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////

        /// Formulation: lower_bound_i <= 1/2 * X' * Qi * X + A(i,:) * X <= upper_bound_i
        /// Note: There can only be one objective function, but there can be m constraints, indexed by i.
        /// TODO: assumed, not tested yet...
        /// Tested:   Q0     AND A with one row (A vector)
        /// Tested:              A matrix alone
        /// UNtested: Q0, Q1 AND A matrix

        //! \brief                      Append a constraint to the system, and optionally a linear line to the constraint matrix A. The line's row_id is the same, as the constr_id of the quadratic constraints.
        //! \param bound_type           The bound type of the i-th constraint (i-th constraint is assumed to mean: Qi and A(i,:) )
        //! \param lower_bound          The lower bound of the i-th constraint (i-th constraint is assumed to mean: Qi and A(i,:) )
        //! \param upper_bound          The upper bound of the i-th constraint (i-th constraint is assumed to mean: Qi and A(i,:) )
        //! \param coeffs               Pointer to a dense vector to be appended to linear constraint matrix A.
        //! \return                     EXIT_SUCCESS/EXIT_FAILURE, if bound!=SG_INFINITY was set to an unbounded type (i.e. SG_BOUND::FREE, lower bound of SG_BOUND::LESS_EQ, or upper bound of SG_BOUND::GREATER_EQ)
        inline int                               addConstraint                  ( BOUND                bound_type
                                                                                , Scalar               lower_bound // = -getINF()
                                                                                , Scalar               upper_bound // = +getINF()
                                                                                , SparseMatrix       * row_vector     = NULL
                                                                                , LINEARITY            c_lin = LINEARITY::LINEAR
                                                                                );
        inline int                               addConstraint                  ( Variable const& constr );

        //! \brief                      Add linear constraint using an std::vector
        inline int                               addLinConstraint               ( BOUND                bound_type
                                                                                , Scalar               lower_bound // = -getINF()
                                                                                , Scalar               upper_bound //  = +getINF()
                                                                                , std::vector<Scalar>  coeffs
                                                                                , LINEARITY            c_lin = LINEARITY::LINEAR
                                                                                );

        inline int                               addLinConstraints              ( SparseMatrix const& mx );                     //!< \brief Adds full linear constraint matrix. Please call addLinConstraint before.
        inline int                               setLinConstraints              ( SparseMatrix const& mx );                     //!< \brief Sets full linear constraint matrix, clears previous lin constraints. Please call addLinConstraint before.

        inline BOUND                      const  getConstraintBoundType         ( int i ) const { return _bkc[i]; }             //!< \brief Returns bound type of i-th constraint (free, lower, upper, range, fixed).
        inline Scalar                     const  getConstraintLowerBound        ( int i ) const { return _blc[i]; }             //!< \brief Returns lower bound of i-th constraint.
        inline Scalar                     const  getConstraintUpperBound        ( int i ) const { return _buc[i]; }             //!< \brief Returns upper bound of i-th constraint.
        inline LINEARITY                  const  getConstraintLinearity         ( int j ) const { return _lin_c[j]; }
        inline size_t                            getConstraintCount             ()        const { return _bkc.size(); }         //!< \brief Returns the number of linear constraint lines currently in the system.
        inline SparseMatrix                      estimateJacobianOfConstraints  ()        const;

        //! \brief addQConstraint       Append an entry to a quadratic constraint matrix. Duplicate entries are summed in a lazy fashion (when <i>update()</i>) is called.
        //! \brief                      The corresponding linear constraints and bounds are stored in A.row( constr_id ), and _b[k|l|u]c[constr_id].
        /// \param constr_id            Append an entry to Qi, where i == constr_id. Create entry by <i>addLinConstraint</i> before a quadratic entry can be added.
        /// \param i                    Row index of entry in Qi.
        /// \param j                    Column index of entry in Qi.
        /// \param coeff                Entry value to insert into Qi.
        /// \return                     EXIT_SUCCESS/EXIT_FAILURE, if constr_id >= getConstraintCount().
        inline int                               addQConstraint                 ( int constr_id, int i, int j, Scalar coeff );
        inline int                               addQConstraints                ( SparseMatrix const& mx );

        inline SparseEntries              const& getLinConstraints              ()        const { return _linConstrList; }      //!< \brief Returns all A entries, unordered.
        inline SparseMatrix                      getLinConstraintsMatrix        ()        const;                                //!< \brief Returns all A entries assembled to a Sparsematrix, temporarily.
        inline SparseEntries              const& getQuadraticConstraints        ( int i ) const { return _quadConstrList[i]; }  //!< \brief Returns all Qi entries, unordered.
        inline std::vector<SparseEntries> const& getQuadraticConstraints        ()        const { return _quadConstrList; }     //!< \brief Returns vector of Qi entries, vector is ordered by i, but the entries are unordered.
        inline SparseMatrix                      getQuadraticConstraintsMatrix  ( int i ) const;                                //!< \brief Returns all Qi entries assembled to a Sparsematrix, temporarily.

        //! \brief Prints inner state before optimization
        inline int                               printProblem                   ()        const;
        //! \brief Setter to store latest solution, to be called from optimize().
        inline void                              setSolution                    ( VectorX const& sol ) { _x = sol; }

        inline int                               write                          ( std::string const& path ) const; //!< \brief Serialize to disk.
        inline int                               read                           ( std::string const& proj_path ); //!< \brief Read from disk.

    protected:
        // Variables: _blx <= X <= _buc
        std::vector<BOUND>          _bkx;                 //!< \brief Variable bound keys
        std::vector<Scalar>         _blx;                 //!< \brief Variable numerical value of lower bounds
        std::vector<Scalar>         _bux;                 //!< \brief Variable numerical value of upper bounds
        std::vector<VAR_TYPE>       _type_x;              //!< \brief Variable type
        std::vector<LINEARITY>      _lin_x;               //!< \brief Variable linearity

        // Objective: 1/2 * X' * Q_o * X + q_o * X + c
        Scalar                      _cfix;                //!< \brief Constant bias to objective function                                   (c)
        std::vector<Scalar>         _linObjs;             //!< \brief Linear variable coeffs (objective function)                           (q_o)
        SparseEntries               _quadObjList;         //!< \brief Accumulator for QObjective values                                     (Q_o)

        // Constraints of the form: _blc <= 1/2 * X' * Qi * X + A(i,:) * X <= _buc //TODO: is it really A(i,:)?
        std::vector<BOUND>          _bkc;                 //!< \brief Constraints bound keys
        std::vector<Scalar>         _blc;                 //!< \brief Constraints numerical value of lower bounds
        std::vector<Scalar>         _buc;                 //!< \brief Constraints numerical value of upper bounds
        std::vector<LINEARITY>      _lin_c;               //!< \brief Constraints linearity
        SparseEntries               _linConstrList;       //!< \brief Stores values to construct the A linear constraint matrix, unordered. (A, A.cols = n = # of vars, A.rows = m = # of constraints)
        std::vector<SparseEntries>  _quadConstrList;      //!< \brief Accumulator for Quadratic constraint values                           [Q0, Q1, ... Qi ... Qm ]

        bool                        _updated;             //!< \brief This has to be flipped to true by calling update() for optimize() to run.
        VectorX                     _x;                   //!< \brief Optimize() stores the solution here when finishing.

        inline std::string getAuxName() const { return "aux.csv"; }
        inline std::string getqoName() const { return "qo.csv"; }
        inline std::string getQoName() const { return "Qo.csv"; }
        inline std::string getAName() const { return "A.csv"; }
        inline std::string getQiName( int i ) const { char id[255]; sprintf( id, "Q%d.csv", i ); return std::string( id ); }
        inline int         _parseAuxFile( std::string const& aux_file_path );
}; // ...class SGOpt

} // ... namespace qcqpcpp

#ifndef QCQPCPP_INC_SGOPTPROBLEM_HPP
#   define QCQPCPP_INC_SGOPTPROBLEM_HPP
#   include "qcqpcpp/impl/optProblem.hpp"
#endif // QCQPCPP_INC_SGOPTPROBLEM_HPP

#include "qcqpcpp/io/io.h"
#include "sys/stat.h"
namespace  qcqpcpp
{
    template <typename _Scalar> int
    OptProblem<_Scalar>::write( std::string const& path ) const
    {
        std::cout << "creating " << path << std::endl;
        mkdir( path.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH );
        std::vector<std::string> paths;

        // vars
        paths.push_back( getAuxName() );
        {
            std::ofstream aux_file( path + "/" + paths.back() );
            if ( !aux_file.is_open() )
            {
                std::cerr << "[" << __func__ << "]: " << "could not open " << path + "/" + paths.back() << std::endl;
                return EXIT_FAILURE;
            }

            aux_file << "# vars, constraints, objective bias (c)" << std::endl;
            aux_file << this->getVarCount() << "," << this->getConstraintCount() << "," << this->getObjectiveBias() << std::endl;
            aux_file << "# bound_type, var_type, lower_boud, upper_bound, linearity" << std::endl;
            for ( int j = 0; j != this->getVarCount(); ++j )
            {
                aux_file << this->getVarBoundType(j) << ","
                         << this->getVarType( j ) << ","
                         << this->getVarLowerBound( j ) << ","
                         << this->getVarUpperBound( j ) << ","
                         << this->getVarLinearity( j )
                         << std::endl;
            } // for each variable
            aux_file << "# bound_type, lower_boud, upper_bound, linearity" << std::endl;
            for ( int i = 0; i != this->getConstraintCount(); ++i )
            {
                aux_file << this->getConstraintBoundType( i ) << ","
                         << this->getConstraintLowerBound( i ) << ","
                         << this->getConstraintUpperBound( i ) << ","
                         << this->getConstraintLinearity( i )
                         << std::endl;
            } // for each variable

            std::cout << "wrote " << paths.back() << std::endl;
        } // vars file

        // qo
        paths.push_back(getqoName());
        io::writeSparseMatrix( this->getLinObjectivesMatrix(), path + "/" + paths.back(), 0 );
        std::cout << "wrote " << paths.back() << std::endl;

        // Qo
        paths.push_back(getQoName());
        io::writeSparseMatrix( this->getQuadraticObjectivesMatrix(), path + "/" + paths.back(), 0 );
        std::cout << "wrote " << paths.back() << std::endl;

        // A
        paths.push_back(getAName());
        io::writeSparseMatrix( this->getLinConstraintsMatrix(), path + "/" + paths.back(), 0 );
        std::cout << "wrote " << paths.back() << std::endl;

        // Qi
        for ( int i = 0; i != this->getConstraintCount(); ++i )
        {
            paths.push_back( getQiName(i) );
            io::writeSparseMatrix( this->getQuadraticConstraintsMatrix(i), path + "/" + paths.back(), 0 );
            std::cout << "wrote " << paths.back() << std::endl;
        }

        std::string proj_path = path + "/problem.proj";
        std::ofstream proj_file( proj_path );
        if ( !proj_file.is_open() )
        {
            std::cerr << "[" << __func__ << "]: " << "could not open " << proj_path << std::endl;
        }
        for ( size_t i = 0; i != paths.size(); ++i )
            proj_file << paths[i] << std::endl;

        proj_file.close();
        std::cout << "project written to " << proj_path << std::endl;

        return EXIT_FAILURE;
    } // ... Optproblem::write

    template <typename _Scalar> int
    OptProblem<_Scalar>::_parseAuxFile( std::string const& aux_file_path )
    {
        int err = EXIT_SUCCESS;

        // open file
        std::ifstream aux_file( aux_file_path );
        if ( !aux_file.is_open() )
        {
            std::cerr << "[" << __func__ << "]: " << "Coult not open aux_file " << aux_file_path << std::endl;
            return EXIT_FAILURE;
        }

        // parse lines
        int         vars        = -1, // "uninited"
                    constraints = -1; // "uninited"
        int         lid         = -1; // -1 == there's an extra line to read first, the header
        std::string line;             // tmp storage of read line
        while ( getline(aux_file, line) )
        {
            // skip comment
            if ( line[0] == '#') continue;

            // parse line
            std::istringstream iss( line );
            std::string        word;
            int                word_id = 0; // line parse state
            Variable           tmp;         // output storage variable
            while ( std::getline(iss, word, ',') )
            {
                if ( vars < 0 ) // header line
                {
                    vars = atoi( word.c_str() );        // parse number of variables
                    std::getline(iss, word, ',');
                    constraints = atoi( word.c_str() ); // parse number of constraints
                    std::getline(iss, word, ',');
                    this->_cfix = atof( word.c_str() ); // parse constant objective function bias

                    std::cout << "bias is now " << this->getObjectiveBias()
                              << " reading " << vars << " vars and " << constraints << " constraints"
                              << std::endl;

                    std::getline( iss, word );          // clear rest of line
                }
                else // variable or constraint line
                {
                    if ( lid < vars ) // variable line
                    {
                        switch ( word_id ) // bound_type, var_type, lower_boud, upper_bound, linearity
                        {
                            case 0: // Variable bound type
                                tmp._bkx = static_cast<BOUND>( atoi(word.c_str()) );
                                break;
                            case 1: // Variable type (cont/integer/binary)
                                tmp._type_x = static_cast<VAR_TYPE>( atoi(word.c_str()) );
                                break;
                            case 2: // Variable lower bound
                                tmp._blx = atof( word.c_str() );
                                break;
                            case 3: // Variable upper bound
                                tmp._bux = atof( word.c_str() );
                                break;
                            case 4: // Variable linearity
                                tmp._lin_x = static_cast<LINEARITY>( atoi(word.c_str()) );
                                break;
                            default:
                                std::cerr << "[" << __func__ << "]: " << "Unknown VAR word_id..." << word_id << std::endl;
                                break;
                        } //...switch word_id
                    } //...if variable
                    else // constraint line
                    {
                        switch( word_id ) // bound_type, lower_boud, upper_bound, linearity
                        {
                            case 0: // Constraint bound type
                                tmp._bkx = static_cast<BOUND>( atoi(word.c_str()) );
                                break;
                            case 1: // Constraint lower bound
                                tmp._blx = atof( word.c_str() );
                                break;
                            case 2: // Constraint upper bound
                                tmp._bux = atof( word.c_str() );
                                break;
                            case 3: // Constraint linearity
                                tmp._lin_x = static_cast<LINEARITY>( atoi(word.c_str()) );
                                break;
                            default:
                                std::cerr << "[" << __func__ << "]: " << "Unknown CNSTR word_id..." << word_id << ", word: " << word << std::endl;
                                break;
                        } //...switch word_id
                    } //...if constraint

                    // increment line parse status
                    ++word_id;
                } //... non-header line
            } //...while line has more words

            // add variable/constraint
            if ( lid >= 0 )
            {
                if ( lid < vars ) this->addVariable  ( tmp ); // variable
                else              this->addConstraint( tmp ); // constraint
            } // if non-header-line

            // increase read line count (did not count comments)
            ++lid;
        } //... while file has more lines

        // close file
        aux_file.close();

        return err;
    } //...OptProblem::_parseAuxFile()

    template <typename _Scalar> int
    OptProblem<_Scalar>::read( std::string const& proj_file_path )
    {
        int err = EXIT_SUCCESS;

        // Parse project path
        std::string proj_path = ".";
        int slash_pos = proj_file_path.rfind("/");
        if ( slash_pos != std::string::npos )
        {
            std::cout << "slahspos: " << slash_pos << ", size: " << proj_file_path.size() << std::endl;
            proj_path = proj_file_path.substr( 0, slash_pos );
        }

        // Open project file
        std::vector<std::string> paths;
        {
            // open
            std::ifstream proj_file( proj_file_path );
            if ( !proj_file.is_open() )
            {
                std::cerr << "[" << __func__ << "]: " << "could not open " << proj_file_path << std::endl;
                return EXIT_FAILURE;
            }

            // parse project file
            std::string line;
            while ( getline(proj_file, line) )
            {
                if ( line[0] == '#') continue;

                paths.push_back( line );
            }
            proj_file.close();
        } //... parse project file

        // Parse project files
        for ( size_t i = 0; i != paths.size(); ++i )
        {
            std::string fname = paths[i];
            std::cout << "reading " << proj_path + "/" + paths[i] << std::endl; fflush(stdout);
            // aux file or sparse matrix
            if ( fname.find("aux") != std::string::npos )
            {
                this->_parseAuxFile( proj_path + "/" + paths[i] );
            }
            else // sparse matrix
            {
                // read
                SparseMatrix mx = io::readSparseMatrix<_Scalar>( proj_path + "/" + paths[i], 0 );

                // save
                if ( fname.find("qo") != std::string::npos ) // lin objective
                {
                    std::cout << "read " << proj_path + "/" + paths[i] << " as lin obj mx" << std::endl;
                    err += this->addLinObjectives( mx );
                    std::cout << "problem varcount : " << this->getVarCount() << std::endl;
                } //...if qo
                else if ( fname.find("Qo") != std::string::npos ) // quadratic objective
                {
                    std::cout << "read " << proj_path + "/" + paths[i] << " as quadratic obj mx" << std::endl;
                    err += this->addQObjectives( mx );
                } //...if Qo
                else if ( fname.find("A") != std::string::npos ) // lin constraint
                {
                    std::cout << "read " << proj_path + "/" + paths[i] << " as linear constraints mx" << std::endl;
                    err += this->addLinConstraints( mx );
                } //...if A
                else if ( fname.find("Q") != std::string::npos ) // Qo will be parsed before, so that's fine
                {
                    std::cout << "read " << proj_path + "/" + paths[i] << " as quadratic constraints mx" << std::endl;
                    err += this->addQConstraints( mx );
                } //...ifelse sparse matrix type
            } //...sparse matrix
        } //...for project file

        // print summary
        this->printProblem();

        return err;
    } // ...OptProblem::read

} // ... namespace qcqpcpp

#endif // QCQPCPP_INC_SGOPTPROBLEM_H
