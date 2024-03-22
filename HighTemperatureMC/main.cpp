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

    InputParamters(const int &num_beads, const double &coupling_strength, const double &simulation_temperature){
        P = num_beads;
        g = coupling_strength;
        T = simulation_temperature;
    }
};

class MCSimulator{
public:

    std::string rootFolder = "/Users/shaeermoeed/Github/DVRPairMC/";
    std::string densityFolder = rootFolder + "ProbabilityDensities/";
    int num_digits = 3;
    int block_num = 1;
    int num_blocks = 1;
    int N = 100;
    int max_states = 10;
    int num_grid_pts = 2 * max_states + 1;
    int mc_steps = 100000;
    int sampler_seed = 123456789;

    std::string dt_folder;

    MCSimulator(){

        // Maybe worth passing in the datetime string using the python wrapper
        auto time = std::chrono::floor<std::chrono::seconds>(std::chrono::system_clock::now());
        auto today = std::chrono::floor<std::chrono::days>(time);
        std::string datetime = std::format("{0:%d_%m_%Y}_{1:%H_%M_%S}",
                                           std::chrono::year_month_day{today},
                                           std::chrono::hh_mm_ss{time-today});
        dt_folder = rootFolder + + "Results/" + datetime;

        if (not std::filesystem::is_directory(dt_folder)){
            std::filesystem::create_directory(dt_folder);
        }
    }

    void simulateMC(const InputParamters &sim_params){

        double g = sim_params.g;
        double T=  sim_params.T;
        int P = sim_params.P;

        std::string parameter_folder_name = "g_" + std::to_string(g).substr(0, num_digits) +
                                            "_T_" + std::to_string(T).substr(0, num_digits) +
                                            "_P_"+ std::to_string(P) + "_N_" + std::to_string(N) +
                                            "_l_" + std::to_string(max_states);

        std::string parameterFolder = dt_folder + "/" + parameter_folder_name;

        if (not std::filesystem::is_directory(parameterFolder)){
            std::filesystem::create_directory(parameterFolder);
        }

        std::string blockFolder = parameterFolder + "/" + "Block_" + std::to_string(block_num);

        if (not std::filesystem::is_directory(blockFolder)){
            std::filesystem::create_directory(blockFolder);
        }

        std::string outFolder = blockFolder + "/";

        if (not std::filesystem::is_directory(outFolder)){
            std::filesystem::create_directory(outFolder);
        }

        MonteCarloIntegrator mcSimulator = MonteCarloIntegrator(N, P, max_states, num_grid_pts,
                                                                T, g,
                                                                densityFolder,
                                                                num_digits, mc_steps,
                                                                block_num, num_blocks, false,
                                                                111, sampler_seed,
                                                                false, outFolder);

        mcSimulator.integrate();

    }

    void iterateParamters(const std::vector<InputParamters> &simulation_params){

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
    /*
    std::chrono::time_point<std::chrono::system_clock> start = std::chrono::system_clock::now();
    std::vector<InputParamters> simulation_params;
    MCSimulator mcSimulations = MCSimulator();
    simulation_params.push_back(InputParamters(1, 0.2, 1.0));
    simulation_params.push_back(InputParamters(2, 0.2, 1.0));
    simulation_params.push_back(InputParamters(6, 0.2, 1.0));

    simulation_params.push_back(InputParamters(1, 1.0, 1.0));
    simulation_params.push_back(InputParamters(2, 1.0, 1.0));
    simulation_params.push_back(InputParamters(6, 1.0, 1.0));

    simulation_params.push_back(InputParamters(1, 2.0, 1.0));
    simulation_params.push_back(InputParamters(2, 2.0, 1.0));
    simulation_params.push_back(InputParamters(6, 2.0, 1.0));

    mcSimulations.iterateParamters(simulation_params);

    std::chrono::time_point<std::chrono::system_clock> end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;

    std::cout << "elapsed time: " << elapsed_seconds.count() << "s\n";
    */

    std::chrono::time_point<std::chrono::system_clock> start = std::chrono::system_clock::now();

    std::vector<InputParamters> simulation_params_1;
    std::vector<InputParamters> simulation_params_2;
    std::vector<InputParamters> simulation_params_3;
    std::vector<InputParamters> simulation_params_4;
    std::vector<InputParamters> simulation_params_5;
    std::vector<InputParamters> simulation_params_6;
    std::vector<InputParamters> simulation_params_7;
    std::vector<InputParamters> simulation_params_8;
    std::vector<InputParamters> simulation_params_9;
    std::vector<InputParamters> simulation_params_10;

    MCSimulator mcSimulations_1 = MCSimulator();
    MCSimulator mcSimulations_2 = MCSimulator();
    MCSimulator mcSimulations_3 = MCSimulator();
    MCSimulator mcSimulations_4 = MCSimulator();
    MCSimulator mcSimulations_5 = MCSimulator();
    MCSimulator mcSimulations_6 = MCSimulator();
    MCSimulator mcSimulations_7 = MCSimulator();
    MCSimulator mcSimulations_8 = MCSimulator();
    MCSimulator mcSimulations_9 = MCSimulator();
    MCSimulator mcSimulations_10 = MCSimulator();

    double sim_T = 0.1;
    int p_1 = 10;
    int p_2 = 20;
    int p_3 = 40;

    simulation_params_1.push_back(InputParamters(p_1, 0.2, sim_T));
    simulation_params_1.push_back(InputParamters(p_2, 0.2, sim_T));
    simulation_params_1.push_back(InputParamters(p_3, 0.2, sim_T));

    simulation_params_2.push_back(InputParamters(p_1, 0.4, sim_T));
    simulation_params_2.push_back(InputParamters(p_2, 0.4, sim_T));
    simulation_params_2.push_back(InputParamters(p_3, 0.4, sim_T));

    simulation_params_3.push_back(InputParamters(p_1, 0.6, sim_T));
    simulation_params_3.push_back(InputParamters(p_2, 0.6, sim_T));
    simulation_params_3.push_back(InputParamters(p_3, 0.6, sim_T));

    simulation_params_4.push_back(InputParamters(p_1, 0.8, sim_T));
    simulation_params_4.push_back(InputParamters(p_2, 0.8, sim_T));
    simulation_params_4.push_back(InputParamters(p_3, 0.8, sim_T));

    simulation_params_5.push_back(InputParamters(p_1, 1.0, sim_T));
    simulation_params_5.push_back(InputParamters(p_2, 1.0, sim_T));
    simulation_params_5.push_back(InputParamters(p_3, 1.0, sim_T));

    simulation_params_6.push_back(InputParamters(p_1, 1.2, sim_T));
    simulation_params_6.push_back(InputParamters(p_2, 1.2, sim_T));
    simulation_params_6.push_back(InputParamters(p_3, 1.2, sim_T));

    simulation_params_7.push_back(InputParamters(p_1, 1.4, sim_T));
    simulation_params_7.push_back(InputParamters(p_2, 1.4, sim_T));
    simulation_params_7.push_back(InputParamters(p_3, 1.4, sim_T));

    simulation_params_8.push_back(InputParamters(p_1, 1.6, sim_T));
    simulation_params_8.push_back(InputParamters(p_2, 1.6, sim_T));
    simulation_params_8.push_back(InputParamters(p_3, 1.6, sim_T));

    simulation_params_9.push_back(InputParamters(p_1, 1.8, sim_T));
    simulation_params_9.push_back(InputParamters(p_2, 1.8, sim_T));
    simulation_params_9.push_back(InputParamters(p_3, 1.8, sim_T));

    simulation_params_10.push_back(InputParamters(p_1, 2.0, sim_T));
    simulation_params_10.push_back(InputParamters(p_2, 2.0, sim_T));
    simulation_params_10.push_back(InputParamters(p_3, 2.0, sim_T));

    std::thread thread_1(&MCSimulator::iterateParamters, &mcSimulations_1, simulation_params_1);
    std::thread thread_2(&MCSimulator::iterateParamters, &mcSimulations_2, simulation_params_2);
    std::thread thread_3(&MCSimulator::iterateParamters, &mcSimulations_3, simulation_params_3);
    std::thread thread_4(&MCSimulator::iterateParamters, &mcSimulations_4, simulation_params_4);
    std::thread thread_5(&MCSimulator::iterateParamters, &mcSimulations_5, simulation_params_5);
    std::thread thread_6(&MCSimulator::iterateParamters, &mcSimulations_6, simulation_params_6);
    std::thread thread_7(&MCSimulator::iterateParamters, &mcSimulations_7, simulation_params_7);
    std::thread thread_8(&MCSimulator::iterateParamters, &mcSimulations_8, simulation_params_8);
    std::thread thread_9(&MCSimulator::iterateParamters, &mcSimulations_9, simulation_params_9);
    std::thread thread_10(&MCSimulator::iterateParamters, &mcSimulations_10, simulation_params_10);

    thread_1.join();
    thread_2.join();
    thread_3.join();
    thread_4.join();
    thread_5.join();
    thread_6.join();
    thread_7.join();
    thread_8.join();
    thread_9.join();
    thread_10.join();

    std::chrono::time_point<std::chrono::system_clock> end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;

    std::cout << "elapsed time: " << elapsed_seconds.count() << "s\n";


    /*
    std::vector<InputParamters> simulation_params;

    for (int i = 0; i < 10; i++){
        simulation_params.push_back(InputParamters(1, 0.2 * (i+1), 1.0));
        simulation_params.push_back(InputParamters(2, 0.2 * (i+1), 1.0));
        simulation_params.push_back(InputParamters(6, 0.2 * (i+1), 1.0));
    }

    MCSimulator mcSimulations = MCSimulator();

    mcSimulations.iterateParamters(simulation_params);
    */

    return 0;
}
