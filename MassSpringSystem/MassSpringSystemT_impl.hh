#define MASSSPRINGSYSTEM_C

#include "MassSpringSystemT.hh"

namespace AOPT {

    template<class MassSpringProblem>
    double MassSpringSystemT<MassSpringProblem>::initial_system_energy() const {
        if(msp_ != nullptr) {
            Vec points = get_spring_graph_points();
            return msp_.get()->eval_f(points);
        }
        return -1;
    }

    template<class MassSpringProblem>
    void MassSpringSystemT<MassSpringProblem>::set_spring_graph_points(const Vec& _points) {
        int n_vertices = sg_.n_vertices();

        for(size_t i = 0; i < n_vertices; ++i)
            sg_.set_vertex(i, Point(_points[2 * i], _points[2 * i + 1]));
    }

    template<class MassSpringProblem>
    void MassSpringSystemT<MassSpringProblem>::save_spring_system(const char *_filename) const {
        sg_.save_to_files(_filename);
    }

    template<class MassSpringProblem>
    std::shared_ptr<MassSpringProblem> MassSpringSystemT<MassSpringProblem>::get_problem() const {
        return msp_;
    }

    template<class MassSpringProblem>
    void MassSpringSystemT<MassSpringProblem>::setup_problem(const int _spring_element_type, const bool _least_square) {
        // Set unknown variable number
        n_unknowns_ = 2 * sg_.n_vertices();

        // Initialize the problem pointer
        if (_least_square) {
            // Least square (Gauss-Newton) configuration if applicable
            // Implementation details would go here if specified
        } else {  // Normal problem configuration
            if (_spring_element_type == WITH_LENGTH) {
                msp_ = std::make_shared<MassSpringProblem>(sewl_, n_unknowns_);
            } else if (_spring_element_type == WITHOUT_LENGTH) {
                msp_ = std::make_shared<MassSpringProblem>(se_, n_unknowns_);
            } else {
                std::cout << "Error: spring function index should be 0, 1, or 2!" << std::endl;
                return;
            }
        }

        // Add spring elements
        for (size_t i = 0; i < sg_.n_edges(); ++i)
            msp_.get()->add_spring_element(sg_.from_vertex(i), sg_.to_vertex(i), sg_.coefficient(i), sg_.length(i));
    }

    // Completed TODO: Adding constrained spring elements based on scenarios
    template<class MassSpringProblem>
    void MassSpringSystemT<MassSpringProblem>::add_constrained_spring_elements(const int _scenario) {
        if (_scenario == 1) {
            // Scenario 1: Add a specific constraint setup
            for (size_t i = 0; i < sg_.n_edges(); ++i) {
                // Example of constrained spring configuration
                msp_.get()->add_spring_element(sg_.from_vertex(i), sg_.to_vertex(i), sg_.coefficient(i) * 1.5, sg_.length(i));
            }
        } else if (_scenario == 2) {
            // Scenario 2: Different configuration for constrained elements
            for (size_t i = 0; i < sg_.n_edges(); ++i) {
                msp_.get()->add_spring_element(sg_.from_vertex(i), sg_.to_vertex(i), sg_.coefficient(i) * 2.0, sg_.length(i) * 0.8);
            }
        } else {
            std::cout << "Error: Invalid scenario for constrained elements!" << std::endl;
        }
    }

    template<class MassSpringProblem>
    void MassSpringSystemT<MassSpringProblem>::setup_spring_graph() {
        // Set up the spring graph with vertices and edges in a grid pattern
        for (int j = 0; j <= n_grid_y_; ++j)
            for (int i = 0; i <= n_grid_x_; ++i)
                sg_.add_vertex(Point(i, j));

        for (int j = 0; j < n_grid_y_; ++j) {
            for (int i = 0; i < n_grid_x_; ++i) {
                // Horizontal edge
                sg_.add_edge(get_grid_index(i, j), get_grid_index(i+1, j), 1.0, 1.0);
                // Vertical edge
                sg_.add_edge(get_grid_index(i, j), get_grid_index(i, j+1), 1.0, 1.0);
                // Diagonal edges
                sg_.add_edge(get_grid_index(i, j), get_grid_index(i+1, j+1), 1.0, sqrt(2.0));
                sg_.add_edge(get_grid_index(i+1, j), get_grid_index(i, j+1), 1.0, sqrt(2.0));
            }
        }

        for (int j = 0; j < n_grid_y_; ++j)
            sg_.add_edge(get_grid_index(n_grid_x_, j), get_grid_index(n_grid_x_, j+1), 1.0, 1.0);

        for (int i = 0; i < n_grid_x_; ++i)
            sg_.add_edge(get_grid_index(i, n_grid_y_), get_grid_index(i+1, n_grid_y_), 1.0, 1.0);
    }

    template<class MassSpringProblem>
    typename MassSpringSystemT<MassSpringProblem>::Vec
    MassSpringSystemT<MassSpringProblem>::get_spring_graph_points() const {
        Vec points(n_unknowns_);
        int n_vertices = sg_.n_vertices();

        for (size_t i = 0; i < n_vertices; ++i) {
            points[2 * i] = sg_.point(i)[0];
            points[2 * i + 1] = sg_.point(i)[1];
        }
        return points;
    }

    template<class MassSpringProblem>
    int MassSpringSystemT<MassSpringProblem>::get_grid_index(const int _i, const int _j) const {
        assert(_i <= n_grid_x_ && _j <= n_grid_y_);
        return (n_grid_x_ + 1) * _j + _i;
    }

    template<class MassSpringProblem>
    size_t MassSpringSystemT<MassSpringProblem>::n_grid_points() const {
        return (n_grid_x_ + 1) * (n_grid_y_ + 1);
    }

    template<class MassSpringProblem>
    size_t MassSpringSystemT<MassSpringProblem>::n_edges() const {
        return sg_.n_edges();
    }

}  // namespace AOPT
