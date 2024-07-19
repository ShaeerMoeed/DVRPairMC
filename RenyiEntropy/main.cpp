#include <iostream>
#include "MonteCarloIntegrator.cpp"
#include <chrono>
#include <filesystem>
#include <thread>

class InputParamters{
public:
    int P;
    double g;
    double T;
    std::string startConfig;

    InputParamters(const int &num_beads, const double &coupling_strength, const double &simulation_temperature,
                   const std::string &start_config){
        P = num_beads;
        g = coupling_strength;
        T = simulation_temperature;
        startConfig = start_config;
    }
};

class MCSimulator{
public:

    std::string rootFolder = "/Users/shaeermoeed/Github/DVRPairMC/";
    std::string densityFolder = rootFolder + "ProbabilityDensities/";
    int num_digits = 3;
    int N_A = 2;
    int N_B = 2;
    // TODO: Use max_states = 10 for prod
    int max_states = 4;
    int num_grid_pts = 2 * max_states + 1;
    int mc_steps = 2000000;
    int sampler_seed = 123456789;
    int sampler_seed_prime = 456789012;
    int swap_seed = 890123456;

    std::string dt_folder;

    MCSimulator(){

        // Maybe worth passing in the datetime string using the python wrapper
        auto time = std::chrono::floor<std::chrono::seconds>(std::chrono::system_clock::now());
        auto today = std::chrono::floor<std::chrono::days>(time);
        std::string datetime = std::format("{0:%d_%m_%Y}_{1:%H_%M_%S}",
                                           std::chrono::year_month_day{today},
                                           std::chrono::hh_mm_ss{time-today});
        // TODO: Remove Test from path for prod
        dt_folder = rootFolder + "Results/Entropy/Test/" + datetime;

        if (not std::filesystem::is_directory(dt_folder)){
            std::filesystem::create_directory(dt_folder);
        }
    }

    void simulateMC(const InputParamters &sim_params) const{

        double g = sim_params.g;
        double T=  sim_params.T;
        int P = sim_params.P;
        std::string startConfig = sim_params.startConfig;

        std::string parameter_folder_name = "g_" + std::to_string(g).substr(0, num_digits) +
                                            "_T_" + std::to_string(T).substr(0, num_digits) +
                                            "_P_"+ std::to_string(P) + "_NA_" + std::to_string(N_A) +
                                            "_NB_" + std::to_string(N_B) + "_l_" + std::to_string(max_states) +
                                            "_Start_" + startConfig;

        std::string parameterFolder = dt_folder + "/" + parameter_folder_name;

        if (not std::filesystem::is_directory(parameterFolder)){
            std::filesystem::create_directory(parameterFolder);
        }

        std::string outFolder = parameterFolder + "/";

        if (not std::filesystem::is_directory(outFolder)){
            std::filesystem::create_directory(outFolder);
        }

        //TODO: Set diagnostics to false for prod
        MonteCarloIntegrator mcSimulator = MonteCarloIntegrator(N_A, N_B,P,
        max_states, num_grid_pts, T, g, 0,
        densityFolder,num_digits,mc_steps,startConfig,
        111, sampler_seed, sampler_seed_prime, swap_seed, false,outFolder);

        mcSimulator.integrateSwappedUnswapped();

        int swappedMoves = mcSimulator.CountSwappedMoves;
        int unswappedMoves = mcSimulator.CountUnswappedMoves;

        std::cout << "P " << P << "\n";
        std::cout << "N " << swappedMoves << "\n";
        std::cout << "D " << unswappedMoves << "\n";

    }

    void iterateParamters(const std::vector<InputParamters> &simulation_params) const{

        for (int i = 0; i < simulation_params.size(); i++){
            simulateMC(simulation_params.at(i));
        }
    }

};



int main() {

    // TODO: Decide whether some post-processing (averages, errors etc.) should be computed by cpp (NO)
    // TODO: Add utilities for reading input file for parameters
    // TODO: Add utilities to main.cpp for creating the parent directory inside a results folder
    // TODO: Write post-processing utilities in python

    std::chrono::time_point<std::chrono::system_clock> start = std::chrono::system_clock::now();

    double sim_T = 0.1;
    std::vector<double> list_g {0.4, 0.5, 0.6, 1.0, 1.1, 1.2, 1.3};
    std::vector<int> list_p {14, 16, 18};

    int g_index = 0;
    std::vector<InputParamters> simulation_params_1;
    for (int & i : list_p){
        simulation_params_1.emplace_back(i, list_g.at(g_index), sim_T, "Zero");
        simulation_params_1.emplace_back(i, list_g.at(g_index), sim_T, "Pi");
    }
    g_index += 1;

    std::vector<InputParamters> simulation_params_2;
    for (int & i : list_p){
        simulation_params_2.emplace_back(i, list_g.at(g_index), sim_T, "Zero");
        simulation_params_2.emplace_back(i, list_g.at(g_index), sim_T, "Pi");
    }
    g_index += 1;

    std::vector<InputParamters> simulation_params_3;
    for (int & i : list_p){
        simulation_params_3.emplace_back(i, list_g.at(g_index), sim_T, "Zero");
        simulation_params_3.emplace_back(i, list_g.at(g_index), sim_T, "Pi");
    }
    g_index += 1;

    std::vector<InputParamters> simulation_params_4;
    for (int & i : list_p){
        simulation_params_4.emplace_back(i, list_g.at(g_index), sim_T, "Zero");
        simulation_params_4.emplace_back(i, list_g.at(g_index), sim_T, "Pi");
    }
    g_index += 1;

    std::vector<InputParamters> simulation_params_5;
    for (int & i : list_p){
        simulation_params_5.emplace_back(i, list_g.at(g_index), sim_T, "Zero");
        simulation_params_5.emplace_back(i, list_g.at(g_index), sim_T, "Pi");
    }
    g_index += 1;

    std::vector<InputParamters> simulation_params_6;
    for (int & i : list_p){
        simulation_params_6.emplace_back(i, list_g.at(g_index), sim_T, "Zero");
        simulation_params_6.emplace_back(i, list_g.at(g_index), sim_T, "Pi");
    }
    g_index += 1;

    std::vector<InputParamters> simulation_params_7;
    for (int & i : list_p){
        simulation_params_7.emplace_back(i, list_g.at(g_index), sim_T, "Zero");
        simulation_params_7.emplace_back(i, list_g.at(g_index), sim_T, "Pi");
    }
    g_index += 1;
    /*
    std::vector<InputParamters> simulation_params_8;
    for (int & i : list_p){
        simulation_params_8.emplace_back(i, list_g.at(g_index), sim_T, "Zero");
        simulation_params_8.emplace_back(i, list_g.at(g_index), sim_T, "Pi");
    }
    g_index += 1;

    std::vector<InputParamters> simulation_params_9;
    for (int & i : list_p){
        simulation_params_9.emplace_back(i, list_g.at(g_index), sim_T, "Zero");
        simulation_params_9.emplace_back(i, list_g.at(g_index), sim_T, "Pi");
    }
    g_index += 1;

    std::vector<InputParamters> simulation_params_10;
    for (int & i : list_p){
        simulation_params_10.emplace_back(i, list_g.at(g_index), sim_T, "Zero");
        simulation_params_10.emplace_back(i, list_g.at(g_index), sim_T, "Pi");
    }
    g_index += 1;

    std::vector<InputParamters> simulation_params_11;
    for (int & i : list_p){
        simulation_params_11.emplace_back(i, list_g.at(g_index), sim_T, "Zero");
        simulation_params_11.emplace_back(i, list_g.at(g_index), sim_T, "Pi");
    }
    */

    MCSimulator mcSimulations_1 = MCSimulator();
    MCSimulator mcSimulations_2 = MCSimulator();
    MCSimulator mcSimulations_3 = MCSimulator();
    MCSimulator mcSimulations_4 = MCSimulator();
    MCSimulator mcSimulations_5 = MCSimulator();
    MCSimulator mcSimulations_6 = MCSimulator();
    MCSimulator mcSimulations_7 = MCSimulator();
    /*
    MCSimulator mcSimulations_8 = MCSimulator();
    MCSimulator mcSimulations_9 = MCSimulator();
    MCSimulator mcSimulations_10 = MCSimulator();
    MCSimulator mcSimulations_11 = MCSimulator();
    */

    std::thread thread_1(&MCSimulator::iterateParamters, &mcSimulations_1, simulation_params_1);
    std::thread thread_2(&MCSimulator::iterateParamters, &mcSimulations_2, simulation_params_2);
    std::thread thread_3(&MCSimulator::iterateParamters, &mcSimulations_3, simulation_params_3);
    std::thread thread_4(&MCSimulator::iterateParamters, &mcSimulations_4, simulation_params_4);
    std::thread thread_5(&MCSimulator::iterateParamters, &mcSimulations_5, simulation_params_5);
    std::thread thread_6(&MCSimulator::iterateParamters, &mcSimulations_6, simulation_params_6);
    std::thread thread_7(&MCSimulator::iterateParamters, &mcSimulations_7, simulation_params_7);

    /*
    std::thread thread_8(&MCSimulator::iterateParamters, &mcSimulations_8, simulation_params_8);
    std::thread thread_9(&MCSimulator::iterateParamters, &mcSimulations_9, simulation_params_9);
    std::thread thread_10(&MCSimulator::iterateParamters, &mcSimulations_10, simulation_params_10);
    std::thread thread_11(&MCSimulator::iterateParamters, &mcSimulations_11, simulation_params_11);
    */

    thread_1.join();
    thread_2.join();
    thread_3.join();
    thread_4.join();
    thread_5.join();
    thread_6.join();
    thread_7.join();

    /*
    thread_8.join();
    thread_9.join();
    thread_10.join();
    thread_11.join();
    */

    std::chrono::time_point<std::chrono::system_clock> end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;

    std::cout << "elapsed time: " << elapsed_seconds.count() << "s\n";

    /*
    std::vector<InputParamters> simulation_params;

    simulation_params.emplace_back(4, 1.0, 0.1, "Zero");
    simulation_params.emplace_back(4, 1.0, 0.1, "Pi");

    MCSimulator mcSimulations = MCSimulator();

    mcSimulations.iterateParamters(simulation_params);
     */


    return 0;
}
