#include "qcqpcpp/bonminOptProblem.h"
#include "Eigen/Dense"

template <typename Scalar>
inline Eigen::Matrix<Scalar,3,1> pos( Eigen::Matrix<Scalar,6,1> const& line )
{
    return line.template head<3>();
}

template <typename Scalar>
inline Eigen::Matrix<Scalar,3,1> normal( Eigen::Matrix<Scalar,6,1> const& line )
{
    return line.template segment<3>(3).cross( Eigen::Matrix<Scalar,3,1>::UnitZ() );
}

template <typename Scalar>
inline Eigen::Matrix<Scalar,3,1> dir( Eigen::Matrix<Scalar,6,1> const& line )
{
    return line.template segment<3>( 3 );
}

template <typename Scalar>
inline Eigen::Matrix<Scalar,3,1> dir( Eigen::Matrix<Scalar,3,1> const& normal )
{
    return Eigen::Matrix<Scalar,3,1>::UnitZ().cross( normal );
}

template <typename Scalar>
inline Scalar getDistance( Eigen::Matrix<Scalar,6,1> const& line, Eigen::Matrix<Scalar,3,1> const& point )
{
    return (pos(line) - point).cross( dir(line) ).norm();
} //...getDistance()

int testQuadraticConstraintsBonmin(int argc, char** argv)
{
    typedef double                      Scalar;
    typedef Eigen::Matrix<Scalar,3,1>   PositionT;
    typedef Eigen::Matrix<Scalar,6,1>   PrimitiveT;
    typedef int                         LineKeyT; // todo: change to std::pair<int,int>
    typedef qcqpcpp::BonminOpt<Scalar>  OptProblemT;
    typedef std::map<int, Scalar>       IntScalarMap;
    int const Dims = 2;
    Scalar const bound = 1.e30;

    // prepare data
    std::vector<PositionT>          points;
    std::vector<PrimitiveT>         lines;
    std::map<int,std::vector<int> > lines_points; // first: lineId, second: [ pointId, ...]
    {
        points = { { .5, 0., 0. }
                 , { .5, .5, 0. }
                 , { .5, 1., 0. }
                 , { 0., .5, 0. }
                 , { .5, .5, 0. }
                 , { 1., .5, 0. }
                 };
        lines = { (PrimitiveT()<< .6, .5, 0., 1., 1., 0.).finished()
                  , (PrimitiveT()<< .6, .4, 0., 1., -1., 0.).finished()
                };

        lines_points[0] = {0,1,2};
        lines_points[1] = {3,4,5};
    }

    std::map< LineKeyT, std::vector<int> > prims_vars; // first: line id, second: var ids for [nx, ny, nz, d]
    IntScalarMap x0;         // first: var_id , second: start value

    OptProblemT problem;
    problem.printSolutionAtEndOfAlgorithm();
    problem.setDebug( true );

    // vars
    for ( int lid = 0; lid != lines.size(); ++lid )
    {
        // nx
        int var_id = problem.addVariable(qcqpcpp::OptProblem<Scalar>::BOUND::FREE, -bound, bound, qcqpcpp::OptProblem<Scalar>::VAR_TYPE::CONTINUOUS);
        prims_vars[lid].push_back( var_id );
        x0[var_id] = normal( lines[lid] )(0);

        // ny
        var_id = problem.addVariable(qcqpcpp::OptProblem<Scalar>::BOUND::FREE, -bound, bound, qcqpcpp::OptProblem<Scalar>::VAR_TYPE::CONTINUOUS);
        prims_vars[lid].push_back( var_id );
        x0[var_id] = normal( lines[lid] )(1);

        // nz
        // var_id = problem.addVariable(qcqpcpp::OptProblem<Scalar>::BOUND::FREE, -bound, bound, qcqpcpp::OptProblem<Scalar>::VAR_TYPE::CONTINUOUS);
        // vars[j].push_back( var_id );
        // x0[var_id] = normal( lines[j] )(2);

        // d
        var_id = problem.addVariable(qcqpcpp::OptProblem<Scalar>::BOUND::FREE, -bound, bound, qcqpcpp::OptProblem<Scalar>::VAR_TYPE::CONTINUOUS);
        prims_vars[lid].push_back( var_id );
        x0[var_id] = getDistance<Scalar>( lines[lid], PositionT::Zero() );
    }

    // constr
    for ( int lid = 0; lid != lines.size(); ++lid )
    {
        OptProblemT::SparseMatrix norm_constraint( problem.getVarCount(), problem.getVarCount() );

        // 1 * nx * nx + 1 * ny * ny + 1 *  nz * nz
        for ( int dim = 0; dim != Dims; ++dim )
            norm_constraint.insert( prims_vars[LineKeyT(lid)][dim]
                                  , prims_vars[LineKeyT(lid)][dim] ) = Scalar( 1. );

        // add constraint instance
        problem.addConstraint  ( OptProblemT::BOUND::EQUAL, /* = 1 */ Scalar(1.), /* <= 1 */ Scalar(1.), /* linear constraint coeffs: */ NULL );
        // add quadratic coefficients
        problem.addQConstraints( norm_constraint );
    }

    // objectives
    {
        Scalar coeff = Scalar( 0. );
        for ( int lid = 0; lid != lines.size(); ++lid )
        {
            for ( int pid_id = 0; pid_id != lines_points[lid].size(); ++pid_id )
            {
                const int pid = lines_points[lid][pid_id];

                // for each dimension: x,y,z
                for ( int dim = 0; dim != Dims; ++dim )
                {
                    // (p_x)^2 . (n_x)^2
                    coeff = points[pid]( dim );                                                                                               // (p_x)
                    coeff *= coeff;                                                                                                           // (p_x)^2
                    problem.addQObjective( prims_vars[LineKeyT(lid)][dim], prims_vars[LineKeyT(lid)][dim], coeff );                       // (p_x)^2 . (n_x)^2

                    // 2 . p_x . n_x . d
                    coeff = Scalar(2.) * points[pid]( dim );                                                                                // 2 . p_x
                    problem.addQObjective( /* n_x: */ prims_vars[LineKeyT(lid)][dim], /* d: */ prims_vars[LineKeyT(lid)][Dims], coeff );  // 2 . p_x . n_x . d
                }

                // d^2
                coeff = Scalar( 1. );
                problem.addQObjective( /* d: */ prims_vars[LineKeyT(lid)][Dims], /* d: */ prims_vars[LineKeyT(lid)][Dims], coeff );

                // 2 . px . py . nx . ny
                coeff = Scalar(2.) * points[pid](0) * points[pid](1);
                problem.addQObjective( prims_vars[LineKeyT(lid)][0], prims_vars[LineKeyT(lid)][1], coeff );

            } //...points
        } //...lines
    } //...objectives

    // starting point
    {
        OptProblemT::SparseMatrix mx( problem.getVarCount(), 1 );
        for ( IntScalarMap::const_iterator it = x0.begin(); it != x0.end(); ++it )
        {
            mx.insert( it->first, 0 ) = it->second;
        }
        std::cout << "[" << __func__ << "]: " << "starting point: " << mx << std::endl;
        problem.setStartingPoint( mx );
    }

    problem.update( true );
    problem.printProblem();
    problem.setAlgorithm( Bonmin::B_BB );
    std::vector<Scalar> x_out;
    problem.optimize( &x_out, qcqpcpp::OptProblem<Scalar>::MINIMIZE );

    for ( int i = 0; i != x_out.size(); ++i )
    {
        std::cout << "x[" << i << "]: " << x_out[i] << "\n";
    }

    std::vector< Eigen::Matrix<Scalar,6,1> > lines0, lines1;
    std::ofstream f0("x0.plot"), f1("x1.plot"), p("points.plot");
    for ( int i = 0; i < x_out.size(); i += Dims+1 )
    {
        Eigen::Matrix<Scalar,3,1> n0,n1;
        n0.setZero();
        n1.setZero();
        for ( int d = 0; d != Dims; ++d )
        {
            n0(d) = x0   .at( i+d );
            n1(d) = x_out.at( i+d );
        }

        lines0.push_back( (PrimitiveT() << Eigen::Matrix<Scalar,3,1>::Zero() - n0 * x0   .at( i+Dims )
                                                        , dir(n0) ).finished()
                        );
        std::cout << "n1: " << n1.transpose()
                  << ", x_out.at( i+Dims ): " << x_out.at( i+Dims )
                  << ", n1 * x_out.at( i+Dims ): " << (n1 * x_out.at( i+Dims )).transpose() << std::endl;
        lines1.push_back( (PrimitiveT() << Eigen::Matrix<Scalar,3,1>::Zero() - n1 * x_out.at( i+Dims )
                                                        , dir(n1) ).finished()
                        );

        f0 << lines0.back().head<3>().transpose() << "\n"
           << (lines0.back().head<3>() + lines0.back().segment<3>(3)).transpose() << "\n\n\n";
        f1 << lines1.back().head<3>().transpose() << "\n"
           << (lines1.back().head<3>() + lines1.back().segment<3>(3)).transpose() << "\n\n\n";
    }
    for ( int pid = 0; pid != points.size(); ++pid )
    {
        p << points[pid].transpose() << "\n";
    }
    f0.close();
    f1.close();
    p.close();

    system( "gnuplot -e \"plot 'x0.plot' w lines title 'initial', 'x1.plot' w lines title 'optimized', 'points.plot' w points\" -p" );


    return 1;
}
