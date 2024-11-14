#include<vector>
#include<random>
#include "ProbabilityDensities.cpp"
#include "Estimators.cpp"

class MonteCarloIntegrator{
private:
    int numRotors;
    int numBeads;
    int numGridPts;

    ProbabilityDensities probabilityTables;
    Estimators observables;

    int simulationSteps;

    std::mt19937_64 initialConfigGenerator;
    std::mt19937_64 sampleGenerator;

    int startConfigurationSetting;
    int initialConfigSeed;
    int samplerSeed;
    std::vector<std::vector<int>> initialConfigurations;

    bool diagnostics;
    std::vector<std::vector<std::vector<int>>> phiConfigurations;

public:
    MonteCarloIntegrator(const int &num_rotors, const int &num_beads, const int &max_states, const int &num_grid_pts,
                         const double &simulation_temperature, const double &coupling_strength,
                         const std::string &density_directory, const int &density_file_digits,
                         const int &sim_steps, const int &block_num, const int &num_blocks, const int &start_positions,
                         const int &random_start_seed, const int &sampler_seed, const bool &track_diagnostics,
                         const std::string &out_dir):
                         probabilityTables(coupling_strength, simulation_temperature, max_states, num_beads,
                                           num_rotors, density_directory, density_file_digits),
                         observables(sim_steps, num_rotors, num_beads, coupling_strength, max_states,
                                     simulation_temperature, block_num, num_blocks, out_dir){

        numRotors = num_rotors;
        numBeads = num_beads;
        startConfigurationSetting = start_positions;
        initialConfigSeed = random_start_seed;
        samplerSeed = sampler_seed;
        numGridPts = num_grid_pts;
        diagnostics = track_diagnostics;

        std::uniform_int_distribution<> initialConfigDist(0, numGridPts-1);

        simulationSteps = sim_steps;
        sampleGenerator.seed(samplerSeed);
        double seedThrowaway = 0.0;
        // TODO: Uncomment for prod
        for (int i = 0; i < 2000; i++){
            seedThrowaway += initialConfigDist(sampleGenerator);
        }


        // TODO: Comment for prod
        //testSampler(1000, out_dir);

        if (startConfigurationSetting == 0){
            initialConfigGenerator.seed(initialConfigSeed);
            seedThrowaway = 0.0;
            for (int i = 0; i < 2000; i++){
                seedThrowaway += initialConfigDist(initialConfigGenerator);
            }
            seedThrowaway /= 2000.0;
            std::cout << "Average of seed throw away values: " << seedThrowaway << "\n";

            std::vector<int> beadInitialConfigs;
            for (int i = 0; i < numBeads+1; i++){
                beadInitialConfigs.clear();
                for (int j = 0; j < numRotors; j++){
                    beadInitialConfigs.push_back(initialConfigDist(initialConfigGenerator));
                }
                initialConfigurations.push_back(beadInitialConfigs);
            }
        }
        else if(startConfigurationSetting == 1){
            std::vector<int> beadInitialConfigs;
            for (int i = 0; i < numBeads+1; i++){
                beadInitialConfigs.clear();
                for (int j = 0; j < numRotors; j++){
                    beadInitialConfigs.push_back(0);
                }
                initialConfigurations.push_back(beadInitialConfigs);
            }
        }
        else if(startConfigurationSetting == -1){
            int pi_config_val = int(floor(double(numGridPts)/double(2)));
            std::vector<int> beadInitialConfigs;
            for (int i = 0; i < numBeads+1; i++){
                beadInitialConfigs.clear();
                for (int j = 0; j < numRotors; j++){
                    beadInitialConfigs.push_back(pi_config_val);
                }
                initialConfigurations.push_back(beadInitialConfigs);
            }
        }
        if (diagnostics){
            phiConfigurations.push_back(initialConfigurations);
        }

        //observables.updateOrientationalCorrelation(initialConfigurations);
        //observables.updateBinderRatio(initialConfigurations);
        //observables.updatePotentialEnergy(initialConfigurations);
        int configSign = probabilityTables.extractConfigSign(initialConfigurations);
        observables.updateAllProps(initialConfigurations, configSign);

    }

    void testSampler(const int &num_samples, const std::string &out_dir){

        std::vector<int> samples;
        for (int i = 0; i < num_samples; i++){
            std::vector<double> stepProbabilities;
            stepProbabilities.reserve(numGridPts);
            for (int j = 0; j < numGridPts; j++){
                stepProbabilities.push_back(1.0);
            }
            std::discrete_distribution<int> distribution(stepProbabilities.begin(), stepProbabilities.end());
            samples.push_back(distribution(sampleGenerator));
        }

        std::string filePath = out_dir + "/" + "Test Sampler Output PIGS.csv";
        std::string header = "Step, Sample\n";

        std::ofstream ofs(filePath);
        ofs << header;
        for (int i = 0; i < num_samples; i++){
            ofs << std::to_string(i+1);
            ofs << "," << std::to_string(samples.at(i));
            ofs << "\n";
        }
    }

    std::vector<std::vector<int>> gibbsSample(const std::vector<std::vector<int>> &step_configs){

        std::vector<std::vector<int>> nextConfigs = step_configs;
        std::vector<double> stepProbabilities;

        for (int j = 0; j < numRotors; j++){
            stepProbabilities.clear();
            for (int k = 0; k < numGridPts; k++){
                nextConfigs.at(0).at(j) = k;
                double prob_k = probabilityTables.pairProbabilityFirstBead(nextConfigs, j);
                stepProbabilities.push_back(prob_k);
            }
            std::discrete_distribution<int> distribution(stepProbabilities.begin(),stepProbabilities.end());
            nextConfigs.at(0).at(j) = distribution(sampleGenerator);
        }

        for (int i = 1; i < numBeads; i++){
            for (int j = 0; j < numRotors; j++){
                stepProbabilities.clear();
                for (int k = 0; k < numGridPts; k++){
                    nextConfigs.at(i).at(j) = k;
                    double prob_k = probabilityTables.pairProbabilityMiddleBead(nextConfigs, i,
                                                                          j);
                    stepProbabilities.push_back(prob_k);
                }
                std::discrete_distribution<int> distribution(stepProbabilities.begin(),stepProbabilities.end());
                nextConfigs.at(i).at(j) = distribution(sampleGenerator);
            }
        }

        for (int j = 0; j < numRotors; j++){
            stepProbabilities.clear();
            for (int k = 0; k < numGridPts; k++){
                nextConfigs.at(numBeads).at(j) = k;
                double prob_k = probabilityTables.pairProbabilityLastBead(nextConfigs, j);
                stepProbabilities.push_back(prob_k);
            }
            std::discrete_distribution<int> distribution(stepProbabilities.begin(),stepProbabilities.end());
            nextConfigs.at(numBeads).at(j) = distribution(sampleGenerator);
        }

        return nextConfigs;
    }

    void integrateWithSgn(){

        std::vector<std::vector<int>> currentStepConfigs = initialConfigurations;
        for (int i = 1; i < simulationSteps; i++){

            std::vector<std::vector<int>> prevStepConfigs = currentStepConfigs;

            currentStepConfigs = gibbsSample(prevStepConfigs);

            if (diagnostics){
                phiConfigurations.push_back(currentStepConfigs);
            }

            //observables.updateOrientationalCorrelation(currentStepConfigs);
            //observables.updateBinderRatio(currentStepConfigs);
            //observables.updatePotentialEnergy(currentStepConfigs);

            // Should have two different code paths (one which keeps track of the sign and one that does not)
            int configSign = probabilityTables.extractConfigSign(currentStepConfigs);
            observables.updateAllProps(currentStepConfigs, configSign);
        }

        observables.outputStepData();

        if (diagnostics){
            observables.outputDiagnostics(phiConfigurations);
        }

    }

    void integrateWithoutSgn(){

        int configSign=1;
        std::vector<std::vector<int>> currentStepConfigs = initialConfigurations;
        for (int i = 1; i < simulationSteps; i++){

            std::vector<std::vector<int>> prevStepConfigs = currentStepConfigs;

            currentStepConfigs = gibbsSample(prevStepConfigs);

            if (diagnostics){
                phiConfigurations.push_back(currentStepConfigs);
            }

            //observables.updateOrientationalCorrelation(currentStepConfigs);
            //observables.updateBinderRatio(currentStepConfigs);
            //observables.updatePotentialEnergy(currentStepConfigs);

            // Should have two different code paths (one which keeps track of the sign and one that does not)
            //int configSign = probabilityTables.extractConfigSign(currentStepConfigs);

            observables.updateAllProps(currentStepConfigs, configSign);
        }

        observables.outputStepData();

        if (diagnostics){
            observables.outputDiagnostics(phiConfigurations);
        }

    }
};
