#include "csp.h"
#include <ostream>
#include <cmath>
#include "ortools/base/logging.h"

namespace csp {
    solver::solver(std::vector<order_t>& orders_list, uint64_t blanks_width) : orders_list{ orders_list }, blanks_width_in_units{ blanks_width }
    {
        /* TODO: Check for empty input parameters */
    }

    /*
        Returns identity matrix of orders_list size.
        For example if orders_list size is 3 it will be:
        1 0 0
        0 1 0
        0 0 1
     */
    std::vector<std::vector<uint64_t>> solver::get_initial_patterns()
    {
        uint64_t orders_num = orders_list.size();
        std::vector<std::vector<uint64_t>> init_patterns( /* init vector of vectors as all zeroes */
            orders_num,
            std::vector<uint64_t>(orders_num));

        for (uint64_t i = 0; i < orders_list.size(); ++i) {
            for (uint64_t j = 0; j < orders_list.size(); ++j) {
                if (j == i) {
                    init_patterns[i][j] = 1;
                }
            }
        }

        return init_patterns;
    }

    solver::result solver::solve_master(std::vector<std::vector<uint64_t>>& patterns)
    {
        /* TODO: Check for empty patterns vector */

        uint64_t patterns_num = patterns[0].size(); /* Assuming the matrix is squared */
        result ret;

        std::unique_ptr<operations_research::MPSolver> solver(
            operations_research::MPSolver::CreateSolver("GLOP_LINEAR_PROGRAMMING"));

        if (!solver) {
            LOG(WARNING) << "GLOP_LINEAR_PROGRAMMING solver unavailable.";
            return ret;
        }

        std::vector<operations_research::MPVariable*> y;

        for (int i = 0; i < patterns_num; ++i) {
            y.push_back(solver->MakeIntVar(0.0, 1000.0, ""));
        }

        /* Create the objective - minimize total blanks used */
        operations_research::MPObjective* const objective = solver->MutableObjective();

        for (int i = 0; i < patterns_num; ++i) {
            objective->SetCoefficient(y[i], 1);
        }

        objective->SetMinimization();

        /* Add constraints that patterns (orders) must be met */
        const double infinity = solver->infinity();
        std::vector<operations_research::MPConstraint*> constraints;
        for (int i = 0; i < patterns_num; ++i) {
            constraints.push_back(solver->MakeRowConstraint(orders_list[i].num, infinity));
            for (int j = 0; j < patterns_num; ++j) {
                constraints.back()->SetCoefficient(y[j], patterns[i][j]);
            }
        }

        ret.result_status = solver->Solve();

        // Check that the problem has an optimal solution.
        if (ret.result_status != operations_research::MPSolver::OPTIMAL) {
            LOG(INFO) << "The problem does not have an optimal solution!";
            if (ret.result_status == operations_research::MPSolver::FEASIBLE) {
                LOG(INFO) << "A potentially suboptimal solution was found";
            }
            else {
                LOG(INFO) << "The solver could not solve the problem.";
                return ret;
            }
        }

        for (auto& x : y) {
            ret.solution_values.push_back((*x).solution_value());
        }

        for (auto& x : constraints) {
            ret.dual_values.push_back((*x).dual_value());
        }

        return ret;
    }
}

int main() {

    std::vector<csp::order_t> orders_list {
        { 20, 700 },
        { 40, 400 },
        { 32, 234 },
        { 256, 23 },
        { 12, 657 },
    };

    std::unique_ptr<csp::solver> solver = std::make_unique<csp::solver>(orders_list, 3000);

    std::vector<std::vector<uint64_t>> initial_patterns = solver->get_initial_patterns();

    solver->solve_master(initial_patterns);

    return EXIT_SUCCESS;
}