#include <iostream>
#include <string>
#include <Utils/StopWatch.hh>

#include <Algorithms/GradientDescent.hh>
#include <Utils/OptimizationStatistic.hh>
#include <Utils/RandomNumberGenerator.hh>
#include <MassSpringSystemT.hh>
#include <Utils/DerivativeChecker.hh>

std::vector<AOPT::GradientDescent::Vec> get_start_points(int n_grid_x, int n_grid_y) {
    std::vector<AOPT::GradientDescent::Vec> start_pts;
    const int n_vertices = (n_grid_x + 1) * (n_grid_y + 1);

    // Generate the first start point: a regular grid with slight random perturbation
    AOPT::GradientDescent::Vec start1(2 * n_vertices);
    for (int i = 0; i < n_vertices; ++i) {
        start1[2 * i] = i / (n_grid_x + 1);
        start1[2 * i + 1] = i % (n_grid_x + 1);
    }
    AOPT::RandomNumberGenerator rng1(-0.2, 0.2);
    start1 += rng1.get_random_nd_vector(2 * n_vertices);
    start_pts.push_back(start1);

    // Generate the second start point: random positions within a larger range
    AOPT::RandomNumberGenerator rng2(-10.0, 10.0);
    start_pts.push_back(rng2.get_random_nd_vector(2 * n_vertices));

    return start_pts;
}

int main(int _argc, const char* _argv[]) {
    // Set fixed parameters for tests
    int n_grid_x = 8;
    int n_grid_y = 8;
    int max_iter = 10000;

    // Loop over both function types (0: without length, 1: with length)
    for (int func_index = 0; func_index <= 1; ++func_index) {
        // Set the filename prefix to indicate function type
        std::string filename_prefix = (func_index == 0) ? "without_length" : "with_length";
        
        // Loop over two constrained scenarios
        for (int scenario = 1; scenario <= 2; ++scenario) {
            std::string scenario_str = "_scenario_" + std::to_string(scenario);

            // Construct the mass-spring system for current function and scenario
            AOPT::MassSpringSystemT<AOPT::MassSpringProblem2DSparse> mss(n_grid_x, n_grid_y, func_index);
            mss.add_constrained_spring_elements(scenario);

            // Statistic for recording optimization process
            auto opt_stat = std::make_unique<AOPT::OptimizationStatistic>(mss.get_problem().get());

            // Generate two start points
            auto start_points = get_start_points(n_grid_x, n_grid_y);

            // Run gradient descent optimization for each start point
            for (size_t i = 0; i < start_points.size(); ++i) {
                mss.set_spring_graph_points(start_points[i]);
                auto initial_energy = mss.initial_system_energy();
                std::cout << "\nStarting Gradient Descent with function index " << func_index 
                          << ", scenario " << scenario << ", start point " << i + 1 << std::endl;
                std::cout << "Initial system energy: " << initial_energy << std::endl;

                // Save initial state for screenshots
                std::string initial_file = filename_prefix + scenario_str + "_start_" + std::to_string(i+1) + "_initial.csv";
                mss.save_spring_system(initial_file.c_str());

                // Run gradient descent
                opt_stat->start_recording();
                AOPT::GradientDescent::Vec optimized_points = AOPT::GradientDescent::solve(opt_stat.get(), start_points[i], 1e-4, max_iter);
                opt_stat->print_statistics();

                // Set optimized points and calculate final energy
                mss.set_spring_graph_points(optimized_points);
                auto final_energy = mss.initial_system_energy();
                std::cout << "Final system energy after optimization: " << final_energy << std::endl;

                // Save optimized state for screenshots
                std::string optimized_file = filename_prefix + scenario_str + "_start_" + std::to_string(i+1) + "_optimized.csv";
                mss.save_spring_system(optimized_file.c_str());
            }
        }
    }

    return 0;
}
