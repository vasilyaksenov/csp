#include "csp.h"
#include <ostream>
#include <cmath>
#include "ortools/base/logging.h"
#include <argumentum/argparse-h.h>
#include "rapidcsv/src/rapidcsv.h"

using namespace std;
using namespace operations_research;
using namespace argumentum;

namespace csp {
    solver::solver(std::vector<order_t>& orders_list, uint64_t blanks_width) : orders_list{ orders_list }, blanks_width{ blanks_width }
    {
        /* TODO: Check for empty input parameters
                 Check if any order of quantity 1 is greater than blanks */
    }

    void solver::print_patterns(solver::result& result, std::vector<std::vector<uint64_t>>& patterns)
    {
        size_t sol_val_num = result.solution_values.size();
        size_t patterns_len = patterns.size();
        size_t sol_val_j = 0;
        uint64_t sum = 0;

        orders_list;
        map<uint64_t, uint64_t> orders;

        for (auto& l : orders_list) {
            if (auto search = orders.find(l.width); search != orders.end()) {
                search->second += l.num; // add quantity if width already exist
            }
            else {
                orders.insert(make_pair(l.width, l.num)); // insert new width
            }            
        }

        cout << "Solution: " << endl;

        for (size_t j = 0; j < sol_val_num; ++j) {
            sol_val_j = static_cast<size_t>(result.solution_values[j]);
            if (sol_val_j > 0) {
                for (size_t z = 0; z < sol_val_j; ++z) {
                    sum = 0;
                    std::string delim = "";
                    std::cout << "[";
                    for (size_t i = 0; i < patterns_len; ++i) {
                        if (patterns[i][j] > 0) {
                            if (_is_exact_cut) {
                                auto order = orders.find(orders_list[i].width);
                                if (order->second >= patterns[i][j]) {
                                    sum += orders_list[i].width * patterns[i][j];
                                    std::cout << delim << orders_list[i].width << " * " << patterns[i][j];
                                    delim = ", ";
                                    order->second -= patterns[i][j];
                                }
                                else if (order->second != 0) {
                                    sum += orders_list[i].width * order->second;
                                    std::cout << delim << orders_list[i].width << " * " << order->second;
                                    delim = ", ";
                                    order->second = 0;
                                }
                            }
                            else {
                                sum += orders_list[i].width * patterns[i][j];
                                std::cout << delim << orders_list[i].width << " * " << patterns[i][j];
                                delim = ", ";
                            }
                        }                        
                    }
                    std::cout << "] : [ " << blanks_width - sum << " ]" << std::endl;
                }
            }
        }
    }

    void solver::set_iterations(optional<int64_t> iterations)
    {
        _iterations = iterations;
    }

    void solver::set_exact_cut(bool is_exact_cut)
    {
        _is_exact_cut = is_exact_cut;
    }

    void solver::solve_large_model()
    {
        /* TODO: check if orders_list is set! */

        size_t orders_num = orders_list.size();
        std::vector<std::vector<uint64_t>> patterns = get_initial_patterns();
        solver::result result;
        int z = 0;
        
        if (_iterations.has_value()) {
            uint64_t iterations = _iterations.value();
            solver::result result_new;
            for (int i = 0; i < iterations; ++i) {
                result = solve_master(patterns, false);
                result_new = get_new_pattern(result);

                /* add i-th cut of new pattern to i-thp pattern */
                for (uint64_t j = 0; j < orders_list.size(); ++j) {
                    patterns[j].push_back(static_cast<uint64_t>(std::round(result_new.solution_values[j])));
                }
            }
        }
        else {
            solver::result result_new;
            bool run = true;
            vector<double> solution_values_old(orders_num, 0.0);            

            while (run) {
                result = solve_master(patterns, false);
                result_new = get_new_pattern(result);

                /* add i-th cut of new pattern to i-thp pattern */
                for (uint64_t j = 0; j < orders_list.size(); ++j) {
                    patterns[j].push_back(static_cast<uint64_t>(std::round(result_new.solution_values[j])));
                }

                /* if new pattern is equal to old one then end iteration */
                if (result_new.solution_values == solution_values_old) {
                    run = false;
                }
                solution_values_old = result_new.solution_values;
            }
        }

        result = solve_master(patterns, true);

        print_patterns(result, patterns);
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
        size_t orders_num = orders_list.size();
        /* init vector of vectors as all zeroes */
        std::vector<std::vector<uint64_t>> init_patterns(
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

    /* Cutting stock sub-problem */
    solver::result solver::get_new_pattern(solver::result& last_result)
    {
        result ret;
        std::unique_ptr<operations_research::MPSolver> solver;
        init_solver(solver, true);

        size_t marg_val_num =  last_result.marginal_values.size();

        std::vector<operations_research::MPVariable*> new_pattern;

        for (uint64_t i = 0; i < marg_val_num; ++i) {
            new_pattern.push_back(solver->MakeIntVar(0.0, blanks_width, ""));
        }

        /* Create the objective - maximizes the sum of the values times the number of occurrence of that roll in a pattern */
        operations_research::MPObjective* const objective = solver->MutableObjective();

        for (int i = 0; i < marg_val_num; ++i) {
            objective->SetCoefficient(new_pattern[i], last_result.marginal_values[i]);
        }

        objective->SetMaximization();

        /* Ensuring that the pattern stays within the total width of the blank */
        const double infinity = solver->infinity();
        /* Order_length * Order_quantity <= Blank_width */
        operations_research::MPConstraint* const constraints = solver->MakeRowConstraint(-infinity, blanks_width);
        for (int i = 0; i < marg_val_num; ++i) {
            constraints->SetCoefficient(new_pattern[i], static_cast<double>(orders_list[i].width));
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

        for (auto& x : new_pattern) {
            if (x->integer()) {
                ret.solution_values.push_back((*x).solution_value());
            }
            else {
                ret.solution_values.push_back(ceil((*x).solution_value()));
            }
        }
        
        ret.best_solution = objective->Value();

        return ret;
    }

    /* cutting stock master problem */
    solver::result solver::solve_master(std::vector<std::vector<uint64_t>>& patterns, bool is_integer)
    {
        /* TODO: Check for empty patterns vector */

        size_t patterns_num = patterns.size();
        size_t patterns_len = patterns[0].size();
        result ret;

        std::unique_ptr<operations_research::MPSolver> solver;
        init_solver(solver, is_integer);

        std::vector<operations_research::MPVariable*> y;

        for (int i = 0; i < patterns_len; ++i) {
            y.push_back(solver->MakeIntVar(0.0, 100000.0, ""));
        }

        /* Create the objective - minimize total blanks used */
        operations_research::MPObjective* const objective = solver->MutableObjective();

        for (int i = 0; i < patterns_len; ++i) {
            objective->SetCoefficient(y[i], 1);
        }

        objective->SetMinimization();

        /* Add constraints that patterns must be met */
        const double infinity = solver->infinity();
        std::vector<operations_research::MPConstraint*> constraints;
        for (int i = 0; i < patterns_num; ++i) {

                constraints.push_back(solver->MakeRowConstraint(orders_list[i].num,
                                       infinity)); /* cut additional pieces if blank have space left */

                //constraints.push_back(solver->MakeRowConstraint(orders_list[i].num,
                //    orders_list[i].num)); /* cut exactly the number of pieces ordered */

            for (int j = 0; j < patterns_len; ++j) {
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
            ret.solution_values.push_back(static_cast<uint64_t>(round((*x).solution_value())));
        }

        for (auto& x : constraints) {
            if (is_integer) {
                ret.marginal_values.push_back(0);
            }
            else {
                ret.marginal_values.push_back((*x).dual_value());
            }
            
        }

        return ret;
    }

    void solver::init_solver(std::unique_ptr<operations_research::MPSolver>& solver, bool is_integer)
    {
        if (is_integer) {
            solver = std::unique_ptr<operations_research::MPSolver>(
                operations_research::MPSolver::CreateSolver("CBC_MIXED_INTEGER_PROGRAMMING"));
        }
        else {
            solver = std::unique_ptr<operations_research::MPSolver>(
                operations_research::MPSolver::CreateSolver("GLOP_LINEAR_PROGRAMMING"));
        }

        if (!solver) {
            LOG(WARNING) << "Solver unavailable.";
        }
    }
}

int main(int argc, char** argv) {
    string filename;
    int blank_width;
    optional<int64_t> iterations;
    bool exact_cut = false;

    /* Parse input parameters */
    auto parser = argument_parser{};
    auto params = parser.params();
    parser.config().program(argv[0]).description("Cutting Stock Promblem");
    params.add_parameter(blank_width, "--width", "-w")
        .help("Blank width")
        .nargs(1)
        .required(true);
    params.add_parameter(filename, "--file", "-f")
        .help("Input file name")
        .nargs(1)
        .required(true);
    params.add_parameter(iterations, "--iterations", "-i")
        .help("Simple constraint for computation time (default: until get best solution)")
        .nargs(1);
    params.add_parameter(exact_cut, "--exact", "-e")
        .nargs(0)
        .help("Cut exactly as demanded (default: use leftovers to cut random pieces from the order)");

    if (!parser.parse_args(argc, argv, 1))
        return 1;

    cout << "Blank width is " << blank_width << "." << endl;

    rapidcsv::Document doc(filename);

    vector<int64_t> quantity = doc.GetColumn<int64_t>("quantity");
    vector<int64_t> width = doc.GetColumn<int64_t>("width");

    cout << "Found " << quantity.size() << " order(s)." << endl;

    uint64_t bad_orders_num = 0;
    bool bad_order_f = false;
    std::vector<csp::order_t> orders_list;

    for (uint64_t i = 0; i < quantity.size(); ++i) {
        if (width[i] > blank_width) {
            cout << "Order number " << i << " width " << width[i]
                << " doesn't' fit into blank: " << blank_width << endl;
            bad_order_f = true;
        }
        if (width[i] <= 0) {
            cout << "Order number " << i << " have invalid width: " << width[i] << endl;
            bad_order_f = true;
        }
        if (quantity[i] <= 0) {
            cout << "Order number " << i << " have invalid quantity: " << quantity[i] << endl;
            bad_order_f = true;
        }

        if (bad_order_f) {
            ++bad_orders_num;
            bad_order_f = false;
            continue;
        }

        orders_list.push_back(csp::order_t{
            .num = static_cast<uint64_t>(quantity[i]),
                .width = static_cast<uint64_t>(width[i]) });
    }

    if (bad_orders_num > 0) cout << "Found " << bad_orders_num << " order(s) with invalid parameters." << endl;
    if (orders_list.size() == 0) {
        cout << "No valid orders!" << endl;
        return EXIT_SUCCESS;
    }
    if (orders_list.size() < quantity.size()) {
        cout << endl << "Will use " << orders_list.size() << " valid order(s)." << endl;
    }

    if (iterations.has_value()) {
        cout << "Will generate patterns over " << iterations.value() << " iterations." << endl;
    }
    else {
        cout << "Will generate patterns until get best solution." << endl;
    }

    cout << "Computating..." << endl << endl;

    /* Get solution */
    std::unique_ptr<csp::solver> solver = std::make_unique<csp::solver>(orders_list, blank_width);
    solver->set_iterations(iterations);
    solver->set_exact_cut(exact_cut);
    solver->solve_large_model();

    return EXIT_SUCCESS;
}