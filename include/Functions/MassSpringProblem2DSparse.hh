#pragma once

#include <FunctionBase/FunctionBaseSparse.hh>
#include <FunctionBase/ParametricFunctionBase.hh>
#include <Functions/ConstrainedSpringElement2D.hh>

namespace AOPT {

    class MassSpringProblem2DSparse : public FunctionBaseSparse {
    public:
        using Vec = FunctionBaseSparse::Vec;
        using Mat = Eigen::MatrixXd;

        using Edge = std::pair<int, int>;
        using SMat = FunctionBaseSparse::SMat;
        using T = FunctionBaseSparse::T;

        MassSpringProblem2DSparse(ParametricFunctionBase& _spring, const int _n_unknowns) :
                FunctionBaseSparse(),
            n_(_n_unknowns),
            func_(_spring)
        {
            xe_.resize(func_.n_unknowns());
            ge_.resize(func_.n_unknowns());
            he_.resize(func_.n_unknowns(), func_.n_unknowns());
        }

        ~MassSpringProblem2DSparse() {}

        virtual int n_unknowns() override {
            return n_;
        }

        /** evaluates the spring element's energy, which is the sum of the energy
         * of all its springs.
         *
         * \param _x the problem's springs positions.
         *           It should contain the positions of all nodes of the system.
         *           i.e. (_x[2*i], _x[2*i+1]) is the position of the i-th node
         * \return the sum of the energy of all the springs */
        virtual double eval_f(const Vec &_x) override {
            double energy = 0;
            Vec coeff(2);

            //------------------------------------------------------//
            // TODO: (done!) assemble function values of all spring elements
            // use vector xe_ to store the local coordinates of two nodes of every spring
            // then pass it to func_.eval_f(...)
            for(size_t i = 0; i < springs_.size(); ++i) {
                xe_[0] = _x[2 * springs_[i].first];
                xe_[1] = _x[2 * springs_[i].first + 1];
                xe_[2] = _x[2 * springs_[i].second];
                xe_[3] = _x[2 * springs_[i].second + 1];

                coeff[0] = ks_[i];
                coeff[1] = ls_[i];

                energy += func_.eval_f(xe_, coeff);
            }
            //------------------------------------------------------//

            //------------------------------------------------------//
            // TODO: (done!) assemble function values of all the constrained spring elements
            // use cs_xe_ to store the coordinate of the attached node index
            // and cse_ to compute the attached node's functions
            for(size_t i = 0; i < attached_node_indices_.size(); ++i) {
                cs_xe_[0] = _x[2 * attached_node_indices_[i]];
                cs_xe_[1] = _x[2 * attached_node_indices_[i] + 1];
                
                Vec cs_coeff(1);
                cs_coeff[0] = weights_[i];

                energy += cse_.eval_f(cs_xe_, cs_coeff);

            }
            //------------------------------------------------------//
            
            return energy;
        }

        /** The problem's energy gradient is a composition of the individual gradient
         * of each of its springs.
         *
         * \param _x the problem's springs positions.
         *           It should contain the positions of all nodes of the system.
         *           i.e. (_x[2*i], _x[2*i+1]) is the position of the i-th node
         *
         * \param _g should contain the successive per-node gradients.
         *           i.e. (_g[2*i], _g[2*i+1]) is the sum of gradients of all the
         *           springs connected to the i-th node */
        virtual void eval_gradient(const Vec &_x, Vec &_g) override {
            _g.resize(n_unknowns());
            _g.setZero();
            Vec coeff(2);

            //------------------------------------------------------//
            // TODO: (done!) assemble local gradient vector to the global one
            // use ge_ to store the result of the local gradient
            for(size_t i=0; i<springs_.size(); ++i) {
                xe_[0] = _x[2 * springs_[i].first];
                xe_[1] = _x[2 * springs_[i].first + 1];
                xe_[2] = _x[2 * springs_[i].second];
                xe_[3] = _x[2 * springs_[i].second + 1];

                coeff[0] = ks_[i];
                coeff[1] = ls_[i];
                
                func_.eval_gradient(xe_, coeff, ge_);

                _g[2 * springs_[i].first] += ge_[0];
                _g[2 * springs_[i].first + 1] += ge_[1];
                _g[2 * springs_[i].second] += ge_[2];
                _g[2 * springs_[i].second + 1] += ge_[3];
            }
            //------------------------------------------------------//

            //------------------------------------------------------//
            // TODO: (done!) assemble local gradient vector of all the constrained spring elements to the global one
            // use cs_ge_ to store the gradient of the attached node index
            for(size_t i = 0; i < attached_node_indices_.size(); ++i) {
                cs_xe_[0] = _x[2 * attached_node_indices_[i]];
                cs_xe_[1] = _x[2 * attached_node_indices_[i] + 1];

                Vec cs_coeff(1);
                cs_coeff[0] = weights_[i];

                cse_.eval_gradient(cs_xe_, cs_coeff, cs_ge_);


                _g[2 * attached_node_indices_[i]] += cs_ge_[0];
                _g[2 * attached_node_indices_[i] + 1] += cs_ge_[1];
            }
            //------------------------------------------------------//
        }

        /** The problem's energy Hessian is a composition of the individual Hessian
         * of each of its springs.
         * This should be the same as in MassSpringProblem2DDense except for
         * how the resulting matrix is filled. Please refer to Eigen's official
         * documentation if the usage of sparse matrices is still not clear to you
         *
         * \param _x the problem's springs positions.
         *           It should contain the positions of all nodes of the system.
         *           i.e. (_x[2*i], _x[2*i+1]) is the position of the i-th node **/
        virtual void eval_hessian(const Vec &_x, SMat& _h) override {
            std::vector<T> triplets;
            triplets.reserve(16 * springs_.size());

            Vec coeff(2);
            //------------------------------------------------------//
            // TODO: (done!) assemble local hessian matrix to the global one
            // use he_ to store the local hessian matrix
            for(int i=0; i<ks_.size(); ++i) {
                int id0 = 2 * springs_[i].first;
                int id1 = 2 * springs_[i].first + 1;
                int id2 = 2 * springs_[i].second;
                int id3 = 2 * springs_[i].second + 1;

                xe_[0] = _x[id0];
                xe_[1] = _x[id1];
                xe_[2] = _x[id2];
                xe_[3] = _x[id3];

                coeff[0] = ks_[i];
                coeff[1] = ls_[i];

                func_.eval_hessian(xe_, coeff, he_);

                triplets.emplace_back(id0, id0, he_(0,0));
                triplets.emplace_back(id0, id1, he_(0,1));
                triplets.emplace_back(id0, id2, he_(0,2));
                triplets.emplace_back(id0, id3, he_(0,3));

                triplets.emplace_back(id1, id0, he_(1,0));
                triplets.emplace_back(id1, id1, he_(1,1));
                triplets.emplace_back(id1, id2, he_(1,2));
                triplets.emplace_back(id1, id3, he_(1,3));

                triplets.emplace_back(id2, id0, he_(2,0));
                triplets.emplace_back(id2, id1, he_(2,1));
                triplets.emplace_back(id2, id2, he_(2,2));
                triplets.emplace_back(id2, id3, he_(2,3));

                triplets.emplace_back(id3, id0, he_(3,0));
                triplets.emplace_back(id3, id1, he_(3,1));
                triplets.emplace_back(id3, id2, he_(3,2));
                triplets.emplace_back(id3, id3, he_(3,3));
            }
            //------------------------------------------------------//

            //------------------------------------------------------//
            // TODO: (done!) assemble local hessian matrix of all the constrained spring elements to the global one
            // use cs_he_ to store the hessian of the attached node index
            for (size_t i = 0; i < attached_node_indices_.size(); ++i) {
                int idx0 = 2 * attached_node_indices_[i];
                int idx1 = idx0 + 1;

                cs_xe_[0] = _x[idx0];
                cs_xe_[1] = _x[idx1];

                Vec cs_coeff(1);
                cs_coeff[0] = weights_[i];

                cse_.eval_hessian(cs_xe_, cs_coeff, cs_he_);


                triplets.emplace_back(idx0, idx0, cs_he_(0, 0));
                triplets.emplace_back(idx0, idx1, cs_he_(0, 1));
                triplets.emplace_back(idx1, idx0, cs_he_(1, 0));
                triplets.emplace_back(idx1, idx1, cs_he_(1, 1));
            }
            //------------------------------------------------------//

            _h.resize(n_unknowns(), n_unknowns());
            _h.setZero();
            _h.setFromTriplets(triplets.begin(), triplets.end());
        }

        void add_spring_element(const int _v_idx0, const int _v_idx1, const double _k = 1., const double _l = 1.) {
            if (2 * _v_idx0 > (int) n_ || _v_idx0 < 0 || 2 * _v_idx1 >= (int) n_ || _v_idx1 < 0)
                std::cout << "Warning: invalid spring element was added... " << _v_idx0 << " " << _v_idx1 << std::endl;
            else {
                springs_.emplace_back(_v_idx0, _v_idx1);
                ks_.push_back(_k);
                ls_.push_back(_l);
            }
        }

        void add_constrained_spring_element(const int _v_idx, const double _w = 1., const double _px = 0., const double _py = 0.) {
            if (2 * _v_idx > (int) n_ || _v_idx < 0)
                std::cout << "Warning: invalid node constraint element was added... " << _v_idx << std::endl;
            else {
                attached_node_indices_.push_back(_v_idx);
                weights_.push_back(_w);
                desired_points_.push_back(_px);
                desired_points_.push_back(_py);
            }
        }

    private:
        int n_;
        std::vector<Edge> springs_;
        ParametricFunctionBase& func_;
        std::vector<double> ks_;
        std::vector<double> ls_;
        Vec xe_;
        Vec ge_;
        Mat he_;
        std::vector<int> attached_node_indices_;
        ConstrainedSpringElement2D cse_;
        std::vector<double> weights_;
        std::vector<double> desired_points_;
        Vec cs_xe_;
        Vec cs_ge_;
        Mat cs_he_;
    };
}
