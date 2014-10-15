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
    return line.template segment<3>(3).cross( Eigen::Matrix<Scalar,3,1>::UnitZ() ).normalized();
}

template <typename Scalar>
inline Eigen::Matrix<Scalar,3,1> dir( Eigen::Matrix<Scalar,6,1> const& line )
{
    return line.template segment<3>( 3 ).normalized();
}

template <typename Scalar>
inline Eigen::Matrix<Scalar,3,1> dir( Eigen::Matrix<Scalar,3,1> const& normal )
{
    return Eigen::Matrix<Scalar,3,1>::UnitZ().cross( normal ).normalized();
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
//        points = { { .1, 0., 0. }
//                 , { .1, .5, 0. }
//                 , { .1, 1., 0. }
//                 , { .25, .025, 0. }
//                 , { .5, .05, 0. }
//                 , { .75, .0, 0. }
//                 , { .9, 0., 0. }
//                 , { .92, .6, 0. }
//                 , { .81, .9, 0. }
//                 , { .27, .95, 0. }
//                 , { .45, .9, 0. }
//                 , { .7, .975, 0. }
//                 };
                points = { { .1 , .25, 0. }
                         , { .1 , .5, 0. }
                         , { .1 , .75, 0. }

                         , { .25, .1, 0. }
                         , { .5 , .1, 0. }
                         , { .75, .1, 0. }

                         , { .9 , .25, 0. }
                         , { .9 , .5, 0. }
                         , { .9 , .75, 0. }

                         , { .25, .9, 0. }
                         , { .5 , .9, 0. }
                         , { .75, .9, 0. }
                         };
        lines = { (PrimitiveT()<< .6, .5, 0., 1., 1., 0.).finished()
                , (PrimitiveT()<< .6, .4, 0., 1., -1., 0.).finished()
                , (PrimitiveT()<< .0, .0, 0., -1., 1., 0.).finished()
                , (PrimitiveT()<< .9, .5, 0., 1., 1., 0.).finished()
                };
//        lines = { (PrimitiveT()<< .2, .0, 0., .1, .9, 0.).finished()
//                , (PrimitiveT()<< .0, .1, 0., .9, .1, 0.).finished()
//                , (PrimitiveT()<< .8, .0, 0., .1, .9, 0.).finished()
//                , (PrimitiveT()<< .0, .9, 0., .9, .1, 0.).finished()
//                };
        for( int lid = 0; lid != lines.size(); ++lid )
            lines[lid].segment<3>(3).normalize();

        lines_points[0] = {0,1,2};
        lines_points[1] = {3,4,5};
        lines_points[2] = {6,7,8};
        lines_points[3] = {9,10,11};
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
        //x0[var_id] = normal( lines[lid] )(0);

        // ny
        var_id = problem.addVariable(qcqpcpp::OptProblem<Scalar>::BOUND::FREE, -bound, bound, qcqpcpp::OptProblem<Scalar>::VAR_TYPE::CONTINUOUS);
        prims_vars[lid].push_back( var_id );
        //x0[var_id] = normal( lines[lid] )(1);

        // nz
        // var_id = problem.addVariable(qcqpcpp::OptProblem<Scalar>::BOUND::FREE, -bound, bound, qcqpcpp::OptProblem<Scalar>::VAR_TYPE::CONTINUOUS);
        // vars[j].push_back( var_id );
        // x0[var_id] = normal( lines[j] )(2);

        // d
        var_id = problem.addVariable(qcqpcpp::OptProblem<Scalar>::BOUND::FREE, -bound, bound, qcqpcpp::OptProblem<Scalar>::VAR_TYPE::CONTINUOUS);
        prims_vars[lid].push_back( var_id );
        //x0[var_id] = getDistance<Scalar>( lines[lid], PositionT::Zero() );
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

    // perp constr
    std::vector< std::vector<int> > perps = { {0,1,0}, {0,2,1}, {1,3,1} };
    for ( int i = 0; i != perps.size(); ++i )
    {
        OptProblemT::SparseMatrix perp_constraint( problem.getVarCount(), problem.getVarCount() );
        const int lid = perps[i][0];
        const int lid1 = perps[i][1];
        const Scalar rhs = perps[i][2];

        // nx0 * nx1 + ny0 * ny1 = 0
        for ( int dim = 0; dim != Dims; ++dim )
            perp_constraint.insert( prims_vars[LineKeyT(lid1)][dim]
                                  , prims_vars[LineKeyT(lid)][dim] ) = Scalar( 1. );

        // add constraint instance
        problem.addConstraint  ( OptProblemT::BOUND::EQUAL, /* = 1 */ rhs, /* <= 1 */ rhs, /* linear constraint coeffs: */ NULL );
        // add quadratic coefficients
        problem.addQConstraints( perp_constraint );
    }

    problem.update( true );
    problem.printProblem();
    problem.setAlgorithm( Bonmin::B_BB );

    int err = 1;
    std::vector< Eigen::Matrix<Scalar,6,1> > lines0, lines1, *input = &lines;
    std::vector<Scalar> x_out;
    for ( int iter = 0; iter != 2; ++iter )
    {
        // starting point
        x0.clear();
        for ( int lid = 0; lid != (*input).size(); ++lid )
        {
            PositionT nrm = normal( (*input)[lid] );
            for ( int d = 0; d != Dims; ++d )
            {
                const int var_id = prims_vars[lid][d];
                std::cout << "inserting into " << var_id << std::endl;
                x0[var_id] = nrm(d);
            }

            std::cout << "adddded to " << prims_vars[lid][Dims] << std::endl;
            x0[ prims_vars[lid][Dims] ] = getDistance<Scalar>( lines[lid], PositionT::Zero() );
        }

        {
            OptProblemT::SparseMatrix mx( problem.getVarCount(), 1 );
            for ( IntScalarMap::const_iterator it = x0.begin(); it != x0.end(); ++it )
            {
                mx.insert( it->first, 0 ) = it->second;
            }
            std::cout << "[" << __func__ << "]: " << "starting point: " << mx << std::endl;
            problem.setStartingPoint( mx );
        }

        problem.optimize( &x_out, qcqpcpp::OptProblem<Scalar>::MINIMIZE );

        for ( int i = 0; i != x_out.size(); ++i )
        {
            std::cout << "x[" << i << "]: " << x_out[i] << "\n";
        }

        lines0.clear();
        lines = *input;
        lines1.clear();
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
            n0.normalize();
            n1.normalize();

            lines0.push_back( (PrimitiveT() << Eigen::Matrix<Scalar,3,1>::Zero() - n0 * x0.at( i+Dims )
                               , dir(n0) ).finished()
                              );

            lines1.push_back( (PrimitiveT() << Eigen::Matrix<Scalar,3,1>::Zero() - n1 * x_out.at( i+Dims )
                               , dir(n1) ).finished()
                              );

            // line
            f0 << lines0.back().head<3>().transpose() << "\n"
               << (lines0.back().head<3>() + dir(lines0.back())).transpose() << "\n\n\n";

            // normal
            PositionT p0 = lines0.back().head<3>() + dir(lines0.back()) * .5f;

            f0 << p0.transpose() << "\n"
               << (p0 + normal(lines0.back()) * 0.1f ).transpose() << "\n\n\n";

            f1 << lines1.back().head<3>().transpose() << "\n"
               << (lines1.back().head<3>() + lines1.back().segment<3>(3)).transpose() << "\n\n\n";
        }
        f0.close();
        f1.close();

        for ( int pid = 0; pid != points.size(); ++pid )
        {
            p << points[pid].transpose() << "\n";
        }

        p.close();

        system( "gnuplot -e \"plot 'x0.plot' w lines title 'initial', 'x1.plot' w lines title 'optimized', 'points.plot' w points\" -p" );

        // check
        {
            Scalar dot_deviations = 0.;
            for ( int i = 0; i != perps.size(); ++i )
            {
                // abs( desired dot - dot )
                dot_deviations += std::abs( Scalar(perps[i][2]) - normal( lines1[ perps[i][0] ] ).dot(normal(lines1[ perps[i][1] ])) );
            }
            std::cout << "deviation from prescribed dots: " << dot_deviations << std::endl;

            std::vector<Scalar> data_errors( lines.size() );
            int cnt = 0;
            for ( int lid = 0; lid != lines1.size(); ++lid )
            {
                for ( int pid_id = 0; pid_id != lines_points[lid].size(); ++pid_id )
                {
                    Scalar d = std::abs( getDistance( lines1[ lid ], points[ lines_points[lid][pid_id] ] ) );;
                    std::cout << "distance( lines1[" << lid << "], points[ " << lines_points[lid][pid_id] << "] ) = " << d << std::endl;
                    data_errors[lid] += d;
                    ++cnt;
                }
            }
            std::cout<<"data_errors:";
            Scalar sum = 0.;
            for(size_t vi=0;vi!=data_errors.size();++vi)
            {
                std::cout << data_errors[vi]/Scalar(lines_points[vi].size()) << " ";
                sum += data_errors[vi]/Scalar(lines_points[vi].size());
            }
            std::cout << "\n";
            std::cout << "sum of average errors per line: " << sum << std::endl;
            if ( sum < 1.e-3 )
                err = 0;
        }

        input = &lines1;
    } //for iter

    return err;
}
