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
    int N;
    bool trackSgn;

    InputParamters(const int &num_beads, const double &coupling_strength, const double &simulation_temperature,
                   const int &num_rotors, const bool &track_sgn){
        P = num_beads;
        g = coupling_strength;
        T = simulation_temperature;
        N = num_rotors;
        trackSgn = track_sgn;
    }
};

class MCSimulator{
public:

    std::string rootFolder;
    std::string densityFolder;
    int num_digits = 3;
    int block_num = 1;
    int num_blocks = 1;
    //int N = 150;
    //int max_states = 14;
    int max_states = 10;
    int num_grid_pts = 2 * max_states + 1;
    int mc_steps = 50000;
    //int mc_steps = 10000;
    int sampler_seed = 901234567;
    int startPos = 1;

    std::string dt_folder;

    MCSimulator(const std::string &root_folder, const std::string &datetime_folder){

        rootFolder = root_folder;
        dt_folder = datetime_folder;
        densityFolder = rootFolder + "ProbabilityDensities/";

        // TODO: Remove Test from path for prod
        /*
        dt_folder = rootFolder + "Results/PIGS/Test/" + datetime;
        */

        if (not std::filesystem::is_directory(dt_folder)){
            throw std::runtime_error("dt folder must exist\n");
        }
    }

    void simulateMC(const InputParamters &sim_params){

        double g = sim_params.g;
        double T=  sim_params.T;
        int P = sim_params.P;
        int N = sim_params.N;

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

        // TODO: Uncomment for prod runs (commented for testing entropy)
        MonteCarloIntegrator mcSimulator = MonteCarloIntegrator(N, P, max_states, num_grid_pts,
                                                                T, g,
                                                                densityFolder,
                                                                num_digits, mc_steps,
                                                                block_num, num_blocks, startPos,
                                                                111, sampler_seed,
                                                                false, outFolder);

        // Turning diagnostics on for testing entropy
        // TODO: Remove
        /*
        MonteCarloIntegrator mcSimulator = MonteCarloIntegrator(N, P, max_states, num_grid_pts,
                                                                T, g,
                                                                densityFolder,
                                                                num_digits, mc_steps,
                                                                block_num, num_blocks, startPos,
                                                                111, sampler_seed,
                                                                true, outFolder);
        */
        if (sim_params.trackSgn){
            mcSimulator.integrateWithSgn();
        }
        else {
            mcSimulator.integrateWithoutSgn();
        }

    }

    void iterateParamters(const std::vector<InputParamters> &simulation_params){

        for (int i = 0; i < simulation_params.size(); i++){
            simulateMC(simulation_params.at(i));
        }
    }

};



int main() {

    // Maybe worth passing in the datetime string using the python wrapper
    auto time = std::chrono::floor<std::chrono::seconds>(std::chrono::system_clock::now());
    auto today = std::chrono::floor<std::chrono::days>(time);
    std::string datetime = std::format("{0:%d_%m_%Y}_{1:%H_%M_%S}",
                                       std::chrono::year_month_day{today},
                                       std::chrono::hh_mm_ss{time-today});

    std::string root_folder = "/Users/shaeermoeed/Github/DVRPairMC/";
    std::string dt_folder = root_folder + "Results/PIGS/" + datetime;

    // TODO: Remove Test from path for prod
    /*
    dt_folder = rootFolder + "Results/PIGS/Test/" + datetime;
    */

    if (not std::filesystem::is_directory(dt_folder)){
        std::filesystem::create_directory(dt_folder);
    }

    // TODO: Decide whether some post-processing (averages, errors etc.) should be computed by cpp (NO)
    // TODO: Add utilities for reading input file for parameters
    // TODO: Add utilities to main.cpp for creating the parent directory inside a results folder
    // TODO: Write post-processing utilities in python

    // For Multi-threaded runs (uncomment for production runs and change N and mc_steps back in class above)
    // TODO: Uncomment
    // TODO: It is possible that multiple datetime folders are created because multiple instances of MCSimulator are
    // instantiated. Should make dt_folder in main to prevent that from happening.
    // ------------------------------------------------------------------------------------------------

    /*
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

    MCSimulator mcSimulations_1 = MCSimulator(root_folder, dt_folder);
    MCSimulator mcSimulations_2 = MCSimulator(root_folder, dt_folder);
    MCSimulator mcSimulations_3 = MCSimulator(root_folder, dt_folder);
    MCSimulator mcSimulations_4 = MCSimulator(root_folder, dt_folder);
    MCSimulator mcSimulations_5 = MCSimulator(root_folder, dt_folder);
    MCSimulator mcSimulations_6 = MCSimulator(root_folder, dt_folder);
    MCSimulator mcSimulations_7 = MCSimulator(root_folder, dt_folder);
    MCSimulator mcSimulations_8 = MCSimulator(root_folder, dt_folder);
    MCSimulator mcSimulations_9 = MCSimulator(root_folder, dt_folder);
    MCSimulator mcSimulations_10 = MCSimulator(root_folder, dt_folder);

    double sim_T = 0.1;
    int p_1 = 50;
    int p_2 = 60;
    int p_3 = 70;

    simulation_params_1.push_back(InputParamters(p_1, 0.1, sim_T, 150));
    simulation_params_1.push_back(InputParamters(p_2, 0.1, sim_T, 150));
    simulation_params_1.push_back(InputParamters(p_3, 0.1, sim_T, 150));

    simulation_params_2.push_back(InputParamters(p_1, 0.2, sim_T, 150));
    simulation_params_2.push_back(InputParamters(p_2, 0.2, sim_T, 150));
    simulation_params_2.push_back(InputParamters(p_3, 0.2, sim_T, 150));

    simulation_params_3.push_back(InputParamters(p_1, 0.3, sim_T, 150));
    simulation_params_3.push_back(InputParamters(p_2, 0.3, sim_T, 150));
    simulation_params_3.push_back(InputParamters(p_3, 0.3, sim_T, 150));

    simulation_params_4.push_back(InputParamters(p_1, 0.4, sim_T, 150));
    simulation_params_4.push_back(InputParamters(p_2, 0.4, sim_T, 150));
    simulation_params_4.push_back(InputParamters(p_3, 0.4, sim_T, 150));

    simulation_params_5.push_back(InputParamters(p_1, 0.5, sim_T, 150));
    simulation_params_5.push_back(InputParamters(p_2, 0.5, sim_T, 150));
    simulation_params_5.push_back(InputParamters(p_3, 0.5, sim_T, 150));

    simulation_params_6.push_back(InputParamters(p_1, 0.6, sim_T, 150));
    simulation_params_6.push_back(InputParamters(p_2, 0.6, sim_T, 150));
    simulation_params_6.push_back(InputParamters(p_3, 0.6, sim_T, 150));

    simulation_params_7.push_back(InputParamters(p_1, 0.7, sim_T, 150));
    simulation_params_7.push_back(InputParamters(p_2, 0.7, sim_T, 150));
    simulation_params_7.push_back(InputParamters(p_3, 0.7, sim_T, 150));

    simulation_params_8.push_back(InputParamters(p_1, 0.8, sim_T, 150));
    simulation_params_8.push_back(InputParamters(p_2, 0.8, sim_T, 150));
    simulation_params_8.push_back(InputParamters(p_3, 0.8, sim_T, 150));

    simulation_params_9.push_back(InputParamters(p_1, 0.9, sim_T, 150));
    simulation_params_9.push_back(InputParamters(p_2, 0.9, sim_T, 150));
    simulation_params_9.push_back(InputParamters(p_3, 0.9, sim_T, 150));

    simulation_params_10.push_back(InputParamters(p_1, 1.0, sim_T, 150));
    simulation_params_10.push_back(InputParamters(p_2, 1.0, sim_T, 150));
    simulation_params_10.push_back(InputParamters(p_3, 1.0, sim_T, 150));

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
    */
    /*
    std::vector<InputParamters> simulation_params_1;
    MCSimulator mcSimulations_1 = MCSimulator(root_folder, dt_folder);
    std::vector<InputParamters> simulation_params_2;
    MCSimulator mcSimulations_2 = MCSimulator(root_folder, dt_folder);
    double sim_T = 0.1;
    int p_1 = 60;

    simulation_params_1.push_back(InputParamters(p_1, 1.0, sim_T,
                                                 150, false));
    std::thread thread_1(&MCSimulator::iterateParamters, &mcSimulations_1, simulation_params_1);

    simulation_params_2.push_back(InputParamters(p_1, 2.0, sim_T,
                                                 150, false));
    std::thread thread_2(&MCSimulator::iterateParamters, &mcSimulations_2, simulation_params_2);
    thread_1.join();
    thread_2.join();
    */
    /*
     * Account for sign while computing binder ratio and correlation for g=0.5 and g=0.6
    double sim_T = 0.1;
    int p_1 = 80;
    int p_2 = 90;
    int p_3 = 100;

    std::vector<InputParamters> simulation_params_1;
    std::vector<InputParamters> simulation_params_2;
    std::vector<InputParamters> simulation_params_3;
    std::vector<InputParamters> simulation_params_4;
    std::vector<InputParamters> simulation_params_5;
    std::vector<InputParamters> simulation_params_6;

    MCSimulator mcSimulations_1 = MCSimulator(root_folder, dt_folder);
    MCSimulator mcSimulations_2 = MCSimulator(root_folder, dt_folder);
    MCSimulator mcSimulations_3 = MCSimulator(root_folder, dt_folder);
    MCSimulator mcSimulations_4 = MCSimulator(root_folder, dt_folder);
    MCSimulator mcSimulations_5 = MCSimulator(root_folder, dt_folder);
    MCSimulator mcSimulations_6 = MCSimulator(root_folder, dt_folder);

    simulation_params_1.push_back(InputParamters(p_1, 0.5, sim_T, 150));
    simulation_params_2.push_back(InputParamters(p_2, 0.5, sim_T, 150));
    simulation_params_3.push_back(InputParamters(p_3, 0.5, sim_T, 150));
    simulation_params_4.push_back(InputParamters(p_1, 0.6, sim_T, 150));
    simulation_params_5.push_back(InputParamters(p_2, 0.6, sim_T, 150));
    simulation_params_6.push_back(InputParamters(p_3, 0.6, sim_T, 150));

    std::thread thread_1(&MCSimulator::iterateParamters, &mcSimulations_1, simulation_params_1);
    std::thread thread_2(&MCSimulator::iterateParamters, &mcSimulations_2, simulation_params_2);
    std::thread thread_3(&MCSimulator::iterateParamters, &mcSimulations_3, simulation_params_3);
    std::thread thread_4(&MCSimulator::iterateParamters, &mcSimulations_4, simulation_params_4);
    std::thread thread_5(&MCSimulator::iterateParamters, &mcSimulations_5, simulation_params_5);
    std::thread thread_6(&MCSimulator::iterateParamters, &mcSimulations_6, simulation_params_6);

    thread_1.join();
    thread_2.join();
    thread_3.join();
    thread_4.join();
    thread_5.join();
    thread_6.join();
     */

    // For pair primitive energy comparison
    /*
    double sim_T = 0.1;
    int p_1 = 10;
    int p_2 = 14;
    int p_3 = 18;
    int p_4 = 22;
    int p_5 = 26;
    int p_6 = 30;
    int p_7 = 34;
    int p_8 = 38;

    std::vector<InputParamters> simulation_params_1;
    std::vector<InputParamters> simulation_params_2;
    std::vector<InputParamters> simulation_params_3;
    std::vector<InputParamters> simulation_params_4;
    std::vector<InputParamters> simulation_params_5;
    std::vector<InputParamters> simulation_params_6;
    std::vector<InputParamters> simulation_params_7;
    std::vector<InputParamters> simulation_params_8;

    MCSimulator mcSimulations_1 = MCSimulator(root_folder, dt_folder);
    MCSimulator mcSimulations_2 = MCSimulator(root_folder, dt_folder);
    MCSimulator mcSimulations_3 = MCSimulator(root_folder, dt_folder);
    MCSimulator mcSimulations_4 = MCSimulator(root_folder, dt_folder);
    MCSimulator mcSimulations_5 = MCSimulator(root_folder, dt_folder);
    MCSimulator mcSimulations_6 = MCSimulator(root_folder, dt_folder);
    MCSimulator mcSimulations_7 = MCSimulator(root_folder, dt_folder);
    MCSimulator mcSimulations_8 = MCSimulator(root_folder, dt_folder);

    simulation_params_1.push_back(InputParamters(p_1, 1.0, sim_T, 100));
    simulation_params_2.push_back(InputParamters(p_2, 1.0, sim_T, 100));
    simulation_params_3.push_back(InputParamters(p_3, 1.0, sim_T, 100));
    simulation_params_4.push_back(InputParamters(p_4, 1.0, sim_T, 100));
    simulation_params_5.push_back(InputParamters(p_5, 1.0, sim_T, 100));
    simulation_params_6.push_back(InputParamters(p_6, 1.0, sim_T, 100));
    simulation_params_7.push_back(InputParamters(p_7, 1.0, sim_T, 100));
    simulation_params_8.push_back(InputParamters(p_8, 1.0, sim_T, 100));

    std::thread thread_1(&MCSimulator::iterateParamters, &mcSimulations_1, simulation_params_1);
    std::thread thread_2(&MCSimulator::iterateParamters, &mcSimulations_2, simulation_params_2);
    std::thread thread_3(&MCSimulator::iterateParamters, &mcSimulations_3, simulation_params_3);
    std::thread thread_4(&MCSimulator::iterateParamters, &mcSimulations_4, simulation_params_4);
    std::thread thread_5(&MCSimulator::iterateParamters, &mcSimulations_5, simulation_params_5);
    std::thread thread_6(&MCSimulator::iterateParamters, &mcSimulations_6, simulation_params_6);
    std::thread thread_7(&MCSimulator::iterateParamters, &mcSimulations_7, simulation_params_7);
    std::thread thread_8(&MCSimulator::iterateParamters, &mcSimulations_8, simulation_params_8);

    thread_1.join();
    thread_2.join();
    thread_3.join();
    thread_4.join();
    thread_5.join();
    thread_6.join();
    thread_7.join();
    thread_8.join();
    */


    // For chemical potential
    double sim_T = 0.1;
    int p_1 = 20;
    int p_2 = 40;
    int p_3 = 60;
    double g = 0.5;
    int N_init = 2;

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
    std::vector<InputParamters> simulation_params_11;
    std::vector<InputParamters> simulation_params_12;
    std::vector<InputParamters> simulation_params_13;
    std::vector<InputParamters> simulation_params_14;
    std::vector<InputParamters> simulation_params_15;
    std::vector<InputParamters> simulation_params_16;
    std::vector<InputParamters> simulation_params_17;
    std::vector<InputParamters> simulation_params_18;
    std::vector<InputParamters> simulation_params_19;
    std::vector<InputParamters> simulation_params_20;

    MCSimulator mcSimulations_1 = MCSimulator(root_folder, dt_folder);
    MCSimulator mcSimulations_2 = MCSimulator(root_folder, dt_folder);
    MCSimulator mcSimulations_3 = MCSimulator(root_folder, dt_folder);
    MCSimulator mcSimulations_4 = MCSimulator(root_folder, dt_folder);
    MCSimulator mcSimulations_5 = MCSimulator(root_folder, dt_folder);
    MCSimulator mcSimulations_6 = MCSimulator(root_folder, dt_folder);
    MCSimulator mcSimulations_7 = MCSimulator(root_folder, dt_folder);
    MCSimulator mcSimulations_8 = MCSimulator(root_folder, dt_folder);
    MCSimulator mcSimulations_9 = MCSimulator(root_folder, dt_folder);
    MCSimulator mcSimulations_10 = MCSimulator(root_folder, dt_folder);
    MCSimulator mcSimulations_11 = MCSimulator(root_folder, dt_folder);
    MCSimulator mcSimulations_12 = MCSimulator(root_folder, dt_folder);
    MCSimulator mcSimulations_13 = MCSimulator(root_folder, dt_folder);
    MCSimulator mcSimulations_14 = MCSimulator(root_folder, dt_folder);
    MCSimulator mcSimulations_15 = MCSimulator(root_folder, dt_folder);
    MCSimulator mcSimulations_16 = MCSimulator(root_folder, dt_folder);
    MCSimulator mcSimulations_17 = MCSimulator(root_folder, dt_folder);
    MCSimulator mcSimulations_18 = MCSimulator(root_folder, dt_folder);
    MCSimulator mcSimulations_19 = MCSimulator(root_folder, dt_folder);
    MCSimulator mcSimulations_20 = MCSimulator(root_folder, dt_folder);

    simulation_params_1.push_back(InputParamters(p_1, g, sim_T, N_init+1, false));
    simulation_params_1.push_back(InputParamters(p_2, g, sim_T, N_init+1, false));
    simulation_params_1.push_back(InputParamters(p_3, g, sim_T, N_init+1, false));

    simulation_params_2.push_back(InputParamters(p_1, g, sim_T, N_init+2, false));
    simulation_params_2.push_back(InputParamters(p_2, g, sim_T, N_init+2, false));
    simulation_params_2.push_back(InputParamters(p_3, g, sim_T, N_init+2, false));

    simulation_params_3.push_back(InputParamters(p_1, g, sim_T, N_init+3, false));
    simulation_params_3.push_back(InputParamters(p_2, g, sim_T, N_init+3, false));
    simulation_params_3.push_back(InputParamters(p_3, g, sim_T, N_init+3, false));

    simulation_params_4.push_back(InputParamters(p_1, g, sim_T, N_init+4, false));
    simulation_params_4.push_back(InputParamters(p_2, g, sim_T, N_init+4, false));
    simulation_params_4.push_back(InputParamters(p_3, g, sim_T, N_init+4, false));

    simulation_params_5.push_back(InputParamters(p_1, g, sim_T, N_init+5, false));
    simulation_params_5.push_back(InputParamters(p_2, g, sim_T, N_init+5, false));
    simulation_params_5.push_back(InputParamters(p_3, g, sim_T, N_init+5, false));

    simulation_params_6.push_back(InputParamters(p_1, g, sim_T, N_init+6, false));
    simulation_params_6.push_back(InputParamters(p_2, g, sim_T, N_init+6, false));
    simulation_params_6.push_back(InputParamters(p_3, g, sim_T, N_init+6, false));

    simulation_params_7.push_back(InputParamters(p_1, g, sim_T, N_init+7, false));
    simulation_params_7.push_back(InputParamters(p_2, g, sim_T, N_init+7, false));
    simulation_params_7.push_back(InputParamters(p_3, g, sim_T, N_init+7, false));

    simulation_params_8.push_back(InputParamters(p_1, g, sim_T, N_init+8, false));
    simulation_params_8.push_back(InputParamters(p_2, g, sim_T, N_init+8, false));
    simulation_params_8.push_back(InputParamters(p_3, g, sim_T, N_init+8, false));

    simulation_params_9.push_back(InputParamters(p_1, g, sim_T, N_init+9, false));
    simulation_params_9.push_back(InputParamters(p_2, g, sim_T, N_init+9, false));
    simulation_params_9.push_back(InputParamters(p_3, g, sim_T, N_init+9, false));

    simulation_params_10.push_back(InputParamters(p_1, g, sim_T, N_init+10, false));
    simulation_params_10.push_back(InputParamters(p_2, g, sim_T, N_init+10, false));
    simulation_params_10.push_back(InputParamters(p_3, g, sim_T, N_init+10, false));

    simulation_params_11.push_back(InputParamters(p_1, g, sim_T, N_init+11, false));
    simulation_params_11.push_back(InputParamters(p_2, g, sim_T, N_init+11, false));
    simulation_params_11.push_back(InputParamters(p_3, g, sim_T, N_init+11, false));

    simulation_params_12.push_back(InputParamters(p_1, g, sim_T, N_init+12, false));
    simulation_params_12.push_back(InputParamters(p_2, g, sim_T, N_init+12, false));
    simulation_params_12.push_back(InputParamters(p_3, g, sim_T, N_init+12, false));

    simulation_params_13.push_back(InputParamters(p_1, g, sim_T, N_init+13, false));
    simulation_params_13.push_back(InputParamters(p_2, g, sim_T, N_init+13, false));
    simulation_params_13.push_back(InputParamters(p_3, g, sim_T, N_init+13, false));

    simulation_params_14.push_back(InputParamters(p_1, g, sim_T, N_init+14, false));
    simulation_params_14.push_back(InputParamters(p_2, g, sim_T, N_init+14, false));
    simulation_params_14.push_back(InputParamters(p_3, g, sim_T, N_init+14, false));

    simulation_params_15.push_back(InputParamters(p_1, g, sim_T, N_init+15, false));
    simulation_params_15.push_back(InputParamters(p_2, g, sim_T, N_init+15, false));
    simulation_params_15.push_back(InputParamters(p_3, g, sim_T, N_init+15, false));

    simulation_params_16.push_back(InputParamters(p_1, g, sim_T, N_init+16, false));
    simulation_params_16.push_back(InputParamters(p_2, g, sim_T, N_init+16, false));
    simulation_params_16.push_back(InputParamters(p_3, g, sim_T, N_init+16, false));

    simulation_params_17.push_back(InputParamters(p_1, g, sim_T, N_init+17, false));
    simulation_params_17.push_back(InputParamters(p_2, g, sim_T, N_init+17, false));
    simulation_params_17.push_back(InputParamters(p_3, g, sim_T, N_init+17, false));

    simulation_params_18.push_back(InputParamters(p_1, g, sim_T, N_init+18, false));
    simulation_params_18.push_back(InputParamters(p_2, g, sim_T, N_init+18, false));
    simulation_params_18.push_back(InputParamters(p_3, g, sim_T, N_init+18, false));

    simulation_params_19.push_back(InputParamters(p_1, g, sim_T, N_init+19, false));
    simulation_params_19.push_back(InputParamters(p_2, g, sim_T, N_init+19, false));
    simulation_params_19.push_back(InputParamters(p_3, g, sim_T, N_init+19, false));

    simulation_params_20.push_back(InputParamters(p_1, g, sim_T, N_init+20, false));
    simulation_params_20.push_back(InputParamters(p_2, g, sim_T, N_init+20, false));
    simulation_params_20.push_back(InputParamters(p_3, g, sim_T, N_init+20, false));

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
    std::thread thread_11(&MCSimulator::iterateParamters, &mcSimulations_11, simulation_params_11);
    std::thread thread_12(&MCSimulator::iterateParamters, &mcSimulations_12, simulation_params_12);
    std::thread thread_13(&MCSimulator::iterateParamters, &mcSimulations_13, simulation_params_13);
    std::thread thread_14(&MCSimulator::iterateParamters, &mcSimulations_14, simulation_params_14);
    std::thread thread_15(&MCSimulator::iterateParamters, &mcSimulations_15, simulation_params_15);
    std::thread thread_16(&MCSimulator::iterateParamters, &mcSimulations_16, simulation_params_16);
    std::thread thread_17(&MCSimulator::iterateParamters, &mcSimulations_17, simulation_params_17);
    std::thread thread_18(&MCSimulator::iterateParamters, &mcSimulations_18, simulation_params_18);
    /*
    std::thread thread_19(&MCSimulator::iterateParamters, &mcSimulations_19, simulation_params_19);
    std::thread thread_20(&MCSimulator::iterateParamters, &mcSimulations_20, simulation_params_20);
    */

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
    thread_11.join();
    thread_12.join();
    thread_13.join();
    thread_14.join();
    thread_15.join();
    thread_16.join();
    thread_17.join();
    thread_18.join();
    /*
    thread_19.join();
    thread_20.join();
    */

    /*
    std::vector<InputParamters> simulation_params_1;
    MCSimulator mcSimulations_1 = MCSimulator();

    double sim_T = 1.0;
    int p_1 = 6;
    double g = 1.0;
    int N = 3;

    simulation_params_1.push_back(InputParamters(p_1, g, sim_T, N));
    std::thread thread_1(&MCSimulator::iterateParamters, &mcSimulations_1, simulation_params_1);
    thread_1.join();
    */


    return 0;
}
