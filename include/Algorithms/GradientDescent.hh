#pragma once

#include <FunctionBase/FunctionBaseSparse.hh>
#include "LineSearch.hh"
#include <iostream>

//== NAMESPACES ===============================================================

namespace AOPT {

    /* Performs a gradient descent on a given problem.
     * This can work with any Problem with a FunctionBase-style interface since
     * the gradient descent method is rather generic mathematically */
    class GradientDescent {
    public:
        typedef FunctionBaseSparse::Vec Vec; ///< Eigen::VectorXd


        /**
         * \param _problem a pointer to a specific Problem, which can be any type that
         *        has the same interface as FunctionBase's (i.e. with eval_f, eval_gradient, etc.)
         * \param _initial_x  the x starting point
         * \param _eps the stopping criterion below which we consider the method
         *             to be done
         * \param _max_iters a capping number of iterations in case you would end-up with a
         *             bad configuration where the successive attempts of finding the
         *             minimum kind of oscillate around the actual minimum without
         *             finding it
         *
         * \return the minimum found by the method. */
        template <class Problem>
        static Vec solve(Problem *_problem, const Vec& _initial_x, const double _eps = 1e-4, const int _max_iters = 1000000) {
            std::cout << "******** Gradient Descent ********" << std::endl;

            // squared epsilon for stopping criterion
            double e2 = _eps * _eps;

            // get starting point
            Vec x = _initial_x;

            // allocate gradient storage
            Vec g(_problem->n_unknowns());
            int iter(0);

            //------------------------------------------------------//
            //TODO: implement the gradient descent
            double alpha = 1.0;     // Initial step size
            double f_x = _problem->eval_f(x);  // Initial function value

            while (iter < _max_iters) {
                // Compute the gradient at the current point
                _problem->eval_gradient(x, g);

                // Check stopping criterion
                if (g.squaredNorm() < e2) {
                    std::cout << "Convergence achieved with gradient norm squared < epsilon^2: " << g.squaredNorm() << std::endl;
                    break;
                }

                // Perform backtracking line search to find a suitable step size
                double f_x = _problem->eval_f(x); // Replace this with your actual evaluation method
                Vec f_x_vec(1);                       // Create a Vec of size 1 to hold the function value
                f_x_vec(0) = f_x;                     // Assign the double to the Vec
                alpha = LineSearch::backtracking_line_search(_problem, x, g, f_x_vec, alpha);



                // Update x using the gradient and step size
                x = x - alpha * g;

                // Update function value and iteration count
                f_x = _problem->eval_f(x);
                iter++;

                // Output progress
                if (iter % 1000 == 0) {
                    std::cout << "Iteration " << iter << " | f(x) = " << f_x << " | ||g|| = " << g.norm() << std::endl;
                }
            }

            // Check if max iterations were reached
            if (iter == _max_iters) {
                std::cout << "Maximum iterations reached without convergence." << std::endl;
            }
            //------------------------------------------------------//

            return x;
        }
    };
}
