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
            uint64_t blanks_width_in_units = 0;
            const double blanks_width_in_percent = 100.0;
            std::vector<order_t> orders_list; /*!< List of orders */

        private:
            struct result {
                std::vector<double> solution_values;
                std::vector<double> dual_values;
                operations_research::MPSolver::ResultStatus result_status = operations_research::MPSolver::NOT_SOLVED;
            };

        public:
            solver(std::vector<order_t>& orders_list, uint64_t blanks_width);

        public:
            std::vector<std::vector<uint64_t>> get_initial_patterns();
            result solve_master(std::vector<std::vector<uint64_t>>& patterns);

    };   
}
