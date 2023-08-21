#pragma once

#include <iostream>
#include <vector>
#include "ortools/linear_solver/linear_solver.h"

namespace csp {

    /** @struct order_t
     *  @brief Basic order structure. Consists of the width of the pieces need to cut and their
     *  @brief number.
     */
    struct order_t {
        uint64_t num;    /*!< Number of the pieces */
        uint64_t width;   /*!< Width of the pieces */
    };


    class solver
    {
        public:            
            uint64_t blanks_width = 0;
            std::vector<order_t> orders_list; /*!< List of orders */

        private:
            uint64_t _iterations = 20; /*  simple constraint for computation time */
            bool _is_exact_cut = false;
            struct result {
                std::vector<double> solution_values;
                std::vector<double> marginal_values;
                double best_solution = 0;
                operations_research::MPSolver::ResultStatus result_status = operations_research::MPSolver::NOT_SOLVED;
            };

        public:
            solver(std::vector<order_t>& orders_list, uint64_t blanks_width);
            solver(std::vector<order_t>& orders_list, uint64_t blanks_width, uint64_t iterations);

            void set_iterations(uint64_t iterations);
            void set_exact_cut(bool is_exact_cut);
            void solve_large_model();

        private:
            void print_patterns(solver::result& result, std::vector<std::vector<uint64_t>>& patterns);
            std::vector<std::vector<uint64_t>> get_initial_patterns();
            result get_new_pattern(solver::result& last_result);
            result solve_master(std::vector<std::vector<uint64_t>>& patterns, bool is_integer);
            void init_solver(std::unique_ptr<operations_research::MPSolver>& solver, bool is_integer);
    };   
}
