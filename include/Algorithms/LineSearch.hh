#pragma once

#include <FunctionBase/FunctionBaseSparse.hh>

//== NAMESPACES ===============================================================

namespace AOPT {

    //== CLASS DEFINITION =========================================================

    class LineSearch {
    public:
        typedef FunctionBaseSparse::Vec Vec;
        typedef FunctionBaseSparse::SMat SMat;

        /** Back-tracking line search method
         *
         * \param _problem a pointer to a specific Problem, which can be any type that
         *        has the same interface as FunctionBase's (i.e. with eval_f, eval_gradient, etc.)
         * \param _x starting point of the method. Should be of the same dimension as the Problem's
         * \param _g gradient at the starting point.
         * \param _dx delta x (search direction)
         * \param _t0 initial step of the method
         * \param _alpha and _tau variation constants, as stated by the method's definition
         * \return the final step t computed by the back-tracking line search */
        template <class Problem>
        static double backtracking_line_search(Problem *_problem,
                                               const Vec &_x,
                                               const Vec &_g,
                                               const Vec &_dx,
                                               const double _t0,
                                               const double _alpha = 0.5,
                                               const double _tau = 0.75) {
            
            // Initialize step size
            double t = _t0;

            // Evaluate the function at the starting point
            double f_x = _problem->eval_f(_x);

            // Compute the directional derivative (inner product of gradient and direction)
            double grad_dot_dx = _g.dot(_dx);

            // Backtracking loop
            Vec x_new;  // To store the updated x value
            while (true) {
                // Compute x_new = x + t * dx
                x_new = _x + t * _dx;

                // Evaluate the function at the new point
                double f_x_new = _problem->eval_f(x_new);

                // Check the Armijo condition
                if (f_x_new <= f_x + _alpha * t * grad_dot_dx) {
                    break; // Condition met, exit the loop
                }

                // If not, reduce t
                t *= _tau;
            }

            return t; // Return the final step size
        }

    private:
        
    };
    //=============================================================================
}
