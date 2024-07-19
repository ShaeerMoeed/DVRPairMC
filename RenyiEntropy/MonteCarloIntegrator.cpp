#include<vector>
#include<random>
#include "ProbabilityDensities.cpp"
#include "OutputWriter.cpp"

class MonteCarloIntegrator{
private:
    int numRotorsTotal;
    int numRotorsSectorA;
    int numRotorsSectorB;
    int numBeads;
    int numGridPts;
    int numberSwappedRotorsInit;
    int numberSwappedRotors;
    bool Swapped;

    ProbabilityDensities probabilityTables;
    OutputWriter OutWriter;

    int simulationSteps;

    std::mt19937_64 initialConfigGenerator;
    std::mt19937_64 initialConfigGeneratorPrime;
    std::mt19937_64 sampleGenerator;
    std::mt19937_64 sampleGeneratorPrime;
    std::mt19937_64 swapGenerator;

    std::string startConfigsSetting;
    int initialConfigSeed;
    int samplerSeed;
    int samplerSeedPrime;
    int swapSeed;
    std::vector<std::vector<int>> initialConfigurationsA;
    std::vector<std::vector<int>> initialConfigurationsB;
    std::vector<std::vector<int>> initialConfigurationsAPrime;
    std::vector<std::vector<int>> initialConfigurationsBPrime;
    std::vector<std::vector<std::vector<int>>> initialConfigurations;

    bool diagnostics;
    std::vector<std::vector<std::vector<std::vector<int>>>> phiConfigurations;

public:
    int CountSwappedMoves;
    int CountUnswappedMoves;
    std::vector<int> CountSwappedMovesVector;
    std::vector<int> CountUnswappedMovesVector;

    MonteCarloIntegrator(const int &num_rotors_A, const int &num_rotors_B, const int &num_beads, const int &max_states,
                     const int &num_grid_pts, const double &simulation_temperature, const double &coupling_strength,
                     const int num_swapped_rotors, const std::string &density_directory, const int &density_file_digits,
                     const int &sim_steps, const std::string &start_config, const int &random_start_seed,
                     const int &sampler_seed, const int &sampler_seed_prime, const int &swap_seed,
                     const bool &track_diagnostics, const std::string &out_dir):
                     probabilityTables(coupling_strength, simulation_temperature, max_states, num_beads,
                     density_directory, density_file_digits),
                     OutWriter(sim_steps, num_rotors_A, num_rotors_B, num_beads,
                     coupling_strength, max_states, simulation_temperature, out_dir){

        numRotorsSectorA = num_rotors_A;
        numRotorsSectorB = num_rotors_B;
        numRotorsTotal = numRotorsSectorA + numRotorsSectorB;

        numberSwappedRotors = num_swapped_rotors;
        numberSwappedRotorsInit = num_swapped_rotors;
        // TODO: Set Swapped to false for prod
        Swapped = false;

        if (numberSwappedRotorsInit != 0){
            if (numRotorsSectorA != numberSwappedRotorsInit){
                throw std::runtime_error("Number of initial swapped rotors (" +
                std::to_string(numberSwappedRotorsInit) +
                ") is not equal to number of rotors in sector A (" + std::to_string(numRotorsSectorA) + ") \n");
            }
        }

        CountSwappedMoves = 0;
        CountUnswappedMoves = 0;
        if (Swapped){
            CountSwappedMoves += 1;
        }
        else{
            CountUnswappedMoves += 1;
        }
        CountSwappedMovesVector.push_back(CountSwappedMoves);
        CountUnswappedMovesVector.push_back(CountUnswappedMoves);

        numBeads = num_beads;
        startConfigsSetting = start_config;
        initialConfigSeed = random_start_seed;
        samplerSeed = sampler_seed;
        samplerSeedPrime = sampler_seed_prime;
        swapSeed = swap_seed;
        numGridPts = num_grid_pts;
        diagnostics = track_diagnostics;

        std::uniform_int_distribution<> initialConfigDist(0, numGridPts-1);
        std::uniform_int_distribution<> swapDist(0, 1);

        simulationSteps = sim_steps;
        sampleGenerator.seed(samplerSeed);
        sampleGeneratorPrime.seed(samplerSeedPrime);
        swapGenerator.seed(swapSeed);
        double seedThrowaway = 0.0;

        // TODO: Uncomment for prod
        for (int i = 0; i < 2000; i++){
            seedThrowaway += initialConfigDist(sampleGenerator);
            seedThrowaway += initialConfigDist(sampleGeneratorPrime);
            seedThrowaway += swapDist(swapGenerator);
        }

        // TODO: Comment for prod
        //testSampler(1000, out_dir);

        std::vector<int> beadInitialConfigsA;
        std::vector<int> beadInitialConfigsB;
        std::vector<int> beadInitialConfigsAPrime;
        std::vector<int> beadInitialConfigsBPrime;
        if (startConfigsSetting == "random"){
            initialConfigGenerator.seed(initialConfigSeed);
            initialConfigGeneratorPrime.seed(initialConfigSeed);
            seedThrowaway = 0.0;
            for (int i = 0; i < 2000; i++){
                seedThrowaway += initialConfigDist(initialConfigGenerator);
                seedThrowaway += initialConfigDist(initialConfigGeneratorPrime);
            }
            seedThrowaway /= 2000.0;
            std::cout << "Average of seed throw away values: " << seedThrowaway << "\n";

            for (int i = 0; i < numBeads+1; i++){
                beadInitialConfigsA.clear();
                beadInitialConfigsB.clear();
                beadInitialConfigsAPrime.clear();
                beadInitialConfigsBPrime.clear();
                for (int j = 0; j < numRotorsSectorA; j++){
                    beadInitialConfigsA.push_back(initialConfigDist(initialConfigGenerator));
                    beadInitialConfigsAPrime.push_back(initialConfigDist(initialConfigGeneratorPrime));
                }
                for (int j = 0; j < numRotorsSectorB; j++){
                    beadInitialConfigsB.push_back(initialConfigDist(initialConfigGenerator));
                    beadInitialConfigsBPrime.push_back(initialConfigDist(initialConfigGeneratorPrime));
                }
                initialConfigurationsA.push_back(beadInitialConfigsA);
                initialConfigurationsB.push_back(beadInitialConfigsB);
                initialConfigurationsAPrime.push_back(beadInitialConfigsAPrime);
                initialConfigurationsBPrime.push_back(beadInitialConfigsBPrime);
            }
        }
        else if (startConfigsSetting == "Zero"){
            for (int i = 0; i < numBeads+1; i++){
                beadInitialConfigsA.clear();
                beadInitialConfigsAPrime.clear();
                beadInitialConfigsB.clear();
                beadInitialConfigsBPrime.clear();
                for (int j = 0; j < numRotorsSectorA; j++){
                    beadInitialConfigsA.push_back(0);
                    beadInitialConfigsAPrime.push_back(0);
                }
                for (int j = 0; j < numRotorsSectorB; j++){
                    beadInitialConfigsB.push_back(0);
                    beadInitialConfigsBPrime.push_back(0);
                }
                initialConfigurationsA.push_back(beadInitialConfigsA);
                initialConfigurationsB.push_back(beadInitialConfigsB);
                initialConfigurationsAPrime.push_back(beadInitialConfigsAPrime);
                initialConfigurationsBPrime.push_back(beadInitialConfigsBPrime);
            }
        }
        else if (startConfigsSetting == "Pi"){
            int pi_config_val = int(floor(double(numGridPts)/double(2)));
            for (int i = 0; i < numBeads+1; i++){
                beadInitialConfigsA.clear();
                beadInitialConfigsAPrime.clear();
                beadInitialConfigsB.clear();
                beadInitialConfigsBPrime.clear();
                for (int j = 0; j < numRotorsSectorA; j++){
                    beadInitialConfigsA.push_back(pi_config_val);
                    beadInitialConfigsAPrime.push_back(pi_config_val);
                }
                for (int j = 0; j < numRotorsSectorB; j++){
                    beadInitialConfigsB.push_back(pi_config_val);
                    beadInitialConfigsBPrime.push_back(pi_config_val);
                }
                initialConfigurationsA.push_back(beadInitialConfigsA);
                initialConfigurationsB.push_back(beadInitialConfigsB);
                initialConfigurationsAPrime.push_back(beadInitialConfigsAPrime);
                initialConfigurationsBPrime.push_back(beadInitialConfigsBPrime);
            }
        }

        initialConfigurations.push_back(initialConfigurationsA);
        initialConfigurations.push_back(initialConfigurationsB);
        initialConfigurations.push_back(initialConfigurationsAPrime);
        initialConfigurations.push_back(initialConfigurationsBPrime);

        if (diagnostics){
            phiConfigurations.push_back(initialConfigurations);
        }
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

        std::string filePath = out_dir + "/" + "Test Sampler Output Entropy.csv";
        std::string header = "Step, Sample\n";

        std::ofstream ofs(filePath);
        ofs << header;
        for (int i = 0; i < num_samples; i++){
            ofs << std::to_string(i+1);
            ofs << "," << std::to_string(samples.at(i));
            ofs << "\n";
        }
    }

    std::vector<int> generateSampleFirstBeadFirstRotor(const int &index_bead_plus, const int &index_rotor_plus,
                                        const int &index_bead_plus_rotor_plus, const int &index_bead_plus_prime,
                                        const int &index_rotor_plus_prime, const int &index_bead_plus_rotor_plus_prime){

        std::vector<double> stepProbabilities;
        std::vector<double> stepProbabilitiesPrime;

        std::vector<int> samples;

        for (int k = 0; k < numGridPts; k++){
            double prob_k = probabilityTables.pairProbBead1Rotor1(k, index_bead_plus, index_rotor_plus,
                                                                                        index_bead_plus_rotor_plus);
            stepProbabilities.push_back(prob_k);

            double prob_k_prime = probabilityTables.pairProbBead1Rotor1(k,
                                                        index_bead_plus_prime,
                                                        index_rotor_plus_prime,
                                                        index_bead_plus_rotor_plus_prime);
            stepProbabilitiesPrime.push_back(prob_k_prime);
        }

        std::discrete_distribution<int> distribution(stepProbabilities.begin(), stepProbabilities.end());
        samples.push_back(distribution(sampleGenerator));

        std::discrete_distribution<int> distribution_prime(stepProbabilitiesPrime.begin(),stepProbabilitiesPrime.end());
        samples.push_back(distribution_prime(sampleGeneratorPrime));

        return samples;
    }

    std::vector<int> generateSampleMidBeadFirstRotor(const int &index_bead_plus, const int &index_bead_minus,
                                                     const int &index_rotor_plus, const int &index_bead_plus_rotor_plus,
                                                     const int &index_bead_minus_rotor_plus, const int &index_bead_plus_prime,
                                                     const int &index_bead_minus_prime, const int &index_rotor_plus_prime,
                                                     const int &index_bead_plus_rotor_plus_prime,
                                                     const int &index_bead_minus_rotor_plus_prime){

        std::vector<double> stepProbabilities;
        std::vector<double> stepProbabilitiesPrime;

        std::vector<int> samples;

        for (int k = 0; k < numGridPts; k++){
            double prob_k = probabilityTables.pairProbBeadMidRotor1(k, index_bead_plus,
                                                                    index_bead_minus,
                                                                    index_rotor_plus, index_bead_plus_rotor_plus,
                                                                    index_bead_minus_rotor_plus);
            stepProbabilities.push_back(prob_k);

            double prob_k_prime = probabilityTables.pairProbBeadMidRotor1(k,
                                                      index_bead_plus_prime,
                                                      index_bead_minus_prime,
                                                      index_rotor_plus_prime,
                                                      index_bead_plus_rotor_plus_prime,
                                                      index_bead_minus_rotor_plus_prime);
            stepProbabilitiesPrime.push_back(prob_k_prime);
        }

        std::discrete_distribution<int> distribution(stepProbabilities.begin(), stepProbabilities.end());
        samples.push_back(distribution(sampleGenerator));

        std::discrete_distribution<int> distribution_prime(stepProbabilitiesPrime.begin(),stepProbabilitiesPrime.end());
        samples.push_back(distribution_prime(sampleGeneratorPrime));

        return samples;
    }

    std::vector<int> generateSampleLastBeadFirstRotor(const int &index_bead_minus, const int &index_rotor_plus,
                                      const int &index_bead_minus_rotor_plus, const int &index_bead_minus_prime,
                                      const int &index_rotor_plus_prime, const int &index_bead_minus_rotor_plus_prime){

        std::vector<double> stepProbabilities;
        std::vector<double> stepProbabilitiesPrime;

        std::vector<int> samples;

        for (int k = 0; k < numGridPts; k++){
            double prob_k = probabilityTables.pairProbBeadLastRotor1(k, index_bead_minus,
                                                                    index_rotor_plus, index_bead_minus_rotor_plus);
            stepProbabilities.push_back(prob_k);

            double prob_k_prime = probabilityTables.pairProbBeadLastRotor1(k,
                                                              index_bead_minus_prime,
                                                              index_rotor_plus_prime,
                                                              index_bead_minus_rotor_plus_prime);
            stepProbabilitiesPrime.push_back(prob_k_prime);
        }

        std::discrete_distribution<int> distribution(stepProbabilities.begin(), stepProbabilities.end());
        samples.push_back(distribution(sampleGenerator));

        std::discrete_distribution<int> distribution_prime(stepProbabilitiesPrime.begin(),stepProbabilitiesPrime.end());
        samples.push_back(distribution_prime(sampleGeneratorPrime));

        return samples;
    }

    std::vector<int> generateSampleFirstBeadMidRotor(const int &index_bead_plus, const int &index_rotor_plus,
                            const int &index_rotor_minus, const int &index_bead_plus_rotor_plus,
                            const int &index_bead_plus_rotor_minus, const int &index_bead_plus_prime,
                            const int &index_rotor_plus_prime, const int &index_rotor_minus_prime,
                            const int &index_bead_plus_rotor_plus_prime, const int &index_bead_plus_rotor_minus_prime){

        std::vector<double> stepProbabilities;
        std::vector<double> stepProbabilitiesPrime;

        std::vector<int> samples;

        for (int k = 0; k < numGridPts; k++){
            double prob_k = probabilityTables.pairProbBead1RotorMid(k, index_bead_plus, index_rotor_plus,
                                                      index_bead_plus_rotor_plus, index_rotor_minus,
                                                      index_bead_plus_rotor_minus);
            stepProbabilities.push_back(prob_k);

            double prob_k_prime = probabilityTables.pairProbBead1RotorMid(k, index_bead_plus_prime,
              index_rotor_plus_prime,index_bead_plus_rotor_plus_prime,
              index_rotor_minus_prime,index_bead_plus_rotor_minus_prime);
            stepProbabilitiesPrime.push_back(prob_k_prime);
        }

        std::discrete_distribution<int> distribution(stepProbabilities.begin(), stepProbabilities.end());
        samples.push_back(distribution(sampleGenerator));

        std::discrete_distribution<int> distribution_prime(stepProbabilitiesPrime.begin(),
                                                           stepProbabilitiesPrime.end());
        samples.push_back(distribution_prime(sampleGeneratorPrime));

        return samples;
    }

    std::vector<int> generateSampleMidBeadMidRotor(const int &index_rotor_plus, const int &index_rotor_minus,
                           const int &index_bead_plus, const int &index_bead_plus_rotor_plus,
                           const int &index_bead_plus_rotor_minus, const int &index_bead_minus,
                           const int &index_bead_minus_rotor_plus, const int &index_bead_minus_rotor_minus,
                           const int &index_rotor_plus_prime, const int &index_rotor_minus_prime,
                           const int &index_bead_plus_prime, const int &index_bead_plus_rotor_plus_prime,
                           const int &index_bead_plus_rotor_minus_prime, const int &index_bead_minus_prime,
                           const int &index_bead_minus_rotor_plus_prime, const int &index_bead_minus_rotor_minus_prime){

        std::vector<double> stepProbabilities;
        std::vector<double> stepProbabilitiesPrime;

        std::vector<int> samples;

        for (int k = 0; k < numGridPts; k++){
            double prob_k = probabilityTables.pairProbBeadMidRotorMid(k, index_rotor_plus, index_rotor_minus,
                                      index_bead_plus, index_bead_plus_rotor_plus, index_bead_plus_rotor_minus, index_bead_minus,
                                      index_bead_minus_rotor_plus, index_bead_minus_rotor_minus);
            stepProbabilities.push_back(prob_k);

            double prob_k_prime = probabilityTables.pairProbBeadMidRotorMid(k, index_rotor_plus_prime,
                    index_rotor_minus_prime,index_bead_plus_prime,
                    index_bead_plus_rotor_plus_prime,
                    index_bead_plus_rotor_minus_prime,index_bead_minus_prime,
                    index_bead_minus_rotor_plus_prime,
                    index_bead_minus_rotor_minus_prime);
            stepProbabilitiesPrime.push_back(prob_k_prime);
        }

        std::discrete_distribution<int> distribution(stepProbabilities.begin(), stepProbabilities.end());
        samples.push_back(distribution(sampleGenerator));

        std::discrete_distribution<int> distribution_prime(stepProbabilitiesPrime.begin(),
                                                           stepProbabilitiesPrime.end());
        samples.push_back(distribution_prime(sampleGeneratorPrime));

        return samples;
    }

    std::vector<int> generateSampleMidBeadMidRotorSwappedA(const int &index_rotor_plus, const int &index_rotor_plus_p,
       const int &index_rotor_minus, const int &index_bead_plus, const int &index_bead_plus_rotor_plus,
       const int &index_bead_plus_rotor_minus, const int &index_bead_minus, const int &index_bead_minus_rotor_plus,
       const int &index_bead_minus_rotor_minus, const int &index_rotor_plus_prime, const int &index_rotor_plus_p_prime,
       const int &index_rotor_minus_prime, const int &index_bead_plus_prime, const int &index_bead_plus_rotor_plus_prime,
       const int &index_bead_plus_rotor_minus_prime, const int &index_bead_minus_prime, const int &index_bead_minus_rotor_plus_prime,
       const int &index_bead_minus_rotor_minus_prime){

        std::vector<double> stepProbabilities;
        std::vector<double> stepProbabilitiesPrime;

        std::vector<int> samples;

        for (int k = 0; k < numGridPts; k++){
            double prob_k = probabilityTables.pairProbBeadMidRotorMidSwappedA(k, index_rotor_plus, index_rotor_plus_p,
                 index_rotor_minus,index_bead_plus, index_bead_plus_rotor_plus, index_bead_plus_rotor_minus, index_bead_minus,
                 index_bead_minus_rotor_plus, index_bead_minus_rotor_minus);
            stepProbabilities.push_back(prob_k);

            double prob_k_prime = probabilityTables.pairProbBeadMidRotorMidSwappedA(k, index_rotor_plus_prime,
                           index_rotor_plus_p_prime, index_rotor_minus_prime,
                           index_bead_plus_prime,index_bead_plus_rotor_plus_prime,
                           index_bead_plus_rotor_minus_prime,index_bead_minus_prime,
                           index_bead_minus_rotor_plus_prime,
                           index_bead_minus_rotor_minus_prime);
            stepProbabilitiesPrime.push_back(prob_k_prime);
        }

        std::discrete_distribution<int> distribution(stepProbabilities.begin(), stepProbabilities.end());
        samples.push_back(distribution(sampleGenerator));

        std::discrete_distribution<int> distribution_prime(stepProbabilitiesPrime.begin(),
                                                           stepProbabilitiesPrime.end());
        samples.push_back(distribution_prime(sampleGeneratorPrime));

        return samples;
    }

    std::vector<int> generateSampleMidBeadMidRotorSwappedB(const int &index_rotor_plus, const int &index_rotor_minus,
       const int &index_rotor_minus_p, const int &index_bead_plus, const int &index_bead_plus_rotor_plus,
       const int &index_bead_plus_rotor_minus, const int &index_bead_minus, const int &index_bead_minus_rotor_plus,
       const int &index_bead_minus_rotor_minus, const int &index_rotor_plus_prime, const int &index_rotor_minus_prime,
       const int &index_rotor_minus_p_prime, const int &index_bead_plus_prime, const int &index_bead_plus_rotor_plus_prime,
       const int &index_bead_plus_rotor_minus_prime, const int &index_bead_minus_prime, const int &index_bead_minus_rotor_plus_prime,
       const int &index_bead_minus_rotor_minus_prime){

        std::vector<double> stepProbabilities;
        std::vector<double> stepProbabilitiesPrime;

        std::vector<int> samples;

        for (int k = 0; k < numGridPts; k++){
            double prob_k = probabilityTables.pairProbBeadMidRotorMidSwappedB(k, index_rotor_plus,
                                                                              index_rotor_minus, index_rotor_minus_p,
                                                                              index_bead_plus, index_bead_plus_rotor_plus,
                                                                              index_bead_plus_rotor_minus, index_bead_minus,
                                                                              index_bead_minus_rotor_plus,
                                                                              index_bead_minus_rotor_minus);
            stepProbabilities.push_back(prob_k);

            double prob_k_prime = probabilityTables.pairProbBeadMidRotorMidSwappedB(k,
                index_rotor_plus_prime,index_rotor_minus_prime,
                index_rotor_minus_p_prime,index_bead_plus_prime,
                index_bead_plus_rotor_plus_prime,
                index_bead_plus_rotor_minus_prime,
                index_bead_minus_prime,
                index_bead_minus_rotor_plus_prime,
                index_bead_minus_rotor_minus_prime);
            stepProbabilitiesPrime.push_back(prob_k_prime);
        }

        std::discrete_distribution<int> distribution(stepProbabilities.begin(), stepProbabilities.end());
        samples.push_back(distribution(sampleGenerator));

        std::discrete_distribution<int> distribution_prime(stepProbabilitiesPrime.begin(),
                                                           stepProbabilitiesPrime.end());
        samples.push_back(distribution_prime(sampleGeneratorPrime));

        return samples;
    }

    std::vector<int> generateSampleLastBeadMidRotor(const int &index_rotor_plus, const int &index_rotor_minus,
                        const int &index_bead_minus, const int &index_bead_minus_rotor_plus,
                        const int &index_bead_minus_rotor_minus,
                        const int &index_rotor_plus_prime, const int &index_rotor_minus_prime,
                        const int &index_bead_minus_prime,
                        const int &index_bead_minus_rotor_plus_prime,
                        const int &index_bead_minus_rotor_minus_prime){

        std::vector<double> stepProbabilities;
        std::vector<double> stepProbabilitiesPrime;

        std::vector<int> samples;

        for (int k = 0; k < numGridPts; k++){
            double prob_k = probabilityTables.pairProbBeadLastRotorMid(k, index_bead_minus, index_rotor_plus,
                                                                       index_bead_minus_rotor_plus, index_rotor_minus,
                                                                       index_bead_minus_rotor_minus);
            stepProbabilities.push_back(prob_k);

            double prob_k_prime = probabilityTables.pairProbBeadLastRotorMid(k,
                                                     index_bead_minus_prime,
                                                     index_rotor_plus_prime,
                                                     index_bead_minus_rotor_plus_prime,
                                                     index_rotor_minus_prime,
                                                     index_bead_minus_rotor_minus_prime);
            stepProbabilitiesPrime.push_back(prob_k_prime);
        }

        std::discrete_distribution<int> distribution(stepProbabilities.begin(), stepProbabilities.end());
        samples.push_back(distribution(sampleGenerator));

        std::discrete_distribution<int> distribution_prime(stepProbabilitiesPrime.begin(),
                                                           stepProbabilitiesPrime.end());
        samples.push_back(distribution_prime(sampleGeneratorPrime));

        return samples;
    }

    std::vector<int> generateSampleFirstBeadLastRotor(const int &index_bead_plus, const int &index_rotor_minus,
                                     const int &index_bead_plus_rotor_minus, const int &index_bead_plus_prime,
                                     const int &index_rotor_minus_prime, const int &index_bead_plus_rotor_minus_prime){

        std::vector<double> stepProbabilities;
        std::vector<double> stepProbabilitiesPrime;

        std::vector<int> samples;

        for (int k = 0; k < numGridPts; k++){
            double prob_k = probabilityTables.pairProbBead1RotorN(k, index_bead_plus, index_rotor_minus,
                                                                    index_bead_plus_rotor_minus);
            stepProbabilities.push_back(prob_k);

            double prob_k_prime = probabilityTables.pairProbBead1RotorN(k,
                                        index_bead_plus_prime,index_rotor_minus_prime,
                                        index_bead_plus_rotor_minus_prime);
            stepProbabilitiesPrime.push_back(prob_k_prime);
        }

        std::discrete_distribution<int> distribution(stepProbabilities.begin(), stepProbabilities.end());
        samples.push_back(distribution(sampleGenerator));

        std::discrete_distribution<int> distribution_prime(stepProbabilitiesPrime.begin(),stepProbabilitiesPrime.end());
        samples.push_back(distribution_prime(sampleGeneratorPrime));

        return samples;
    }

    std::vector<int> generateSampleMidBeadLastRotor(const int &index_bead_plus, const int &index_rotor_minus,
                                          const int &index_bead_plus_rotor_minus, const int &index_bead_minus,
                                          const int &index_bead_minus_rotor_minus, const int &index_bead_plus_prime,
                                          const int &index_rotor_minus_prime, const int &index_bead_plus_rotor_minus_prime,
                                          const int &index_bead_minus_prime,const int &index_bead_minus_rotor_minus_prime){

        std::vector<double> stepProbabilities;
        std::vector<double> stepProbabilitiesPrime;

        std::vector<int> samples;

        for (int k = 0; k < numGridPts; k++){
            double prob_k = probabilityTables.pairProbBeadMidRotorN(k, index_bead_plus, index_bead_minus,
                                                                    index_rotor_minus, index_bead_plus_rotor_minus,
                                                                    index_bead_minus_rotor_minus);
            stepProbabilities.push_back(prob_k);

            double prob_k_prime = probabilityTables.pairProbBeadMidRotorN(k,
                                                      index_bead_plus_prime,
                                                      index_bead_minus_prime,
                                                      index_rotor_minus_prime,
                                                      index_bead_plus_rotor_minus_prime,
                                                      index_bead_minus_rotor_minus_prime);
            stepProbabilitiesPrime.push_back(prob_k_prime);
        }

        std::discrete_distribution<int> distribution(stepProbabilities.begin(), stepProbabilities.end());
        samples.push_back(distribution(sampleGenerator));

        std::discrete_distribution<int> distribution_prime(stepProbabilitiesPrime.begin(),stepProbabilitiesPrime.end());
        samples.push_back(distribution_prime(sampleGeneratorPrime));

        return samples;
    }

    std::vector<int> generateSampleLastBeadLastRotor(const int &index_rotor_minus, const int &index_bead_minus,
                                    const int &index_bead_minus_rotor_minus, const int &index_rotor_minus_prime,
                                    const int &index_bead_minus_prime, const int &index_bead_minus_rotor_minus_prime){

        std::vector<double> stepProbabilities;
        std::vector<double> stepProbabilitiesPrime;

        std::vector<int> samples;

        for (int k = 0; k < numGridPts; k++){
            double prob_k = probabilityTables.pairProbBeadLastRotorN(k, index_bead_minus, index_rotor_minus,
                                                                    index_bead_minus_rotor_minus);
            stepProbabilities.push_back(prob_k);

            double prob_k_prime = probabilityTables.pairProbBeadLastRotorN(k,
                                                       index_bead_minus_prime,
                                                       index_rotor_minus_prime,
                                                       index_bead_minus_rotor_minus_prime);
            stepProbabilitiesPrime.push_back(prob_k_prime);
        }

        std::discrete_distribution<int> distribution(stepProbabilities.begin(), stepProbabilities.end());
        samples.push_back(distribution(sampleGenerator));

        std::discrete_distribution<int> distribution_prime(stepProbabilitiesPrime.begin(),stepProbabilitiesPrime.end());
        samples.push_back(distribution_prime(sampleGeneratorPrime));

        return samples;
    }

    void sampleExchangeSwappedUnswapped(const std::vector<std::vector<int>> &configsA,
                                        const std::vector<std::vector<int>> &configsB,
                                        const std::vector<std::vector<int>> &configsAPrime,
                                        const std::vector<std::vector<int>> &configsBPrime){

        int midBead = numBeads/2;

        int index_mid_bead_first_rotor = configsA.at(midBead).at(0);
        int index_mid_bead_second_rotor = configsA.at(midBead).at(1);
        int index_mid_bead_third_rotor = configsB.at(midBead).at(0);

        int index_mid_bead_plus_first_rotor = configsA.at(midBead+1).at(0);
        int index_mid_bead_plus_second_rotor = configsA.at(midBead+1).at(1);
        int index_mid_bead_plus_third_rotor = configsB.at(midBead+1).at(0);

        int index_mid_bead_first_rotor_prime = configsAPrime.at(midBead).at(0);
        int index_mid_bead_second_rotor_prime = configsAPrime.at(midBead).at(1);
        int index_mid_bead_third_rotor_prime = configsBPrime.at(midBead).at(0);

        int index_mid_bead_plus_first_rotor_prime = configsAPrime.at(midBead+1).at(0);
        int index_mid_bead_plus_second_rotor_prime = configsAPrime.at(midBead+1).at(1);
        int index_mid_bead_plus_third_rotor_prime = configsBPrime.at(midBead+1).at(0);

        std::vector<double> swapProbs = probabilityTables.swapUnswapProbs(index_mid_bead_first_rotor,
          index_mid_bead_second_rotor, index_mid_bead_third_rotor,
          index_mid_bead_plus_first_rotor,
          index_mid_bead_plus_second_rotor,
          index_mid_bead_plus_third_rotor,
          index_mid_bead_first_rotor_prime,
          index_mid_bead_second_rotor_prime,
          index_mid_bead_third_rotor_prime,
          index_mid_bead_plus_first_rotor_prime,
          index_mid_bead_plus_second_rotor_prime,
          index_mid_bead_plus_third_rotor_prime);

        std::discrete_distribution<int> distribution(swapProbs.begin(), swapProbs.end());
        if (distribution(swapGenerator) == 0){
            Swapped = true;
        }
        else {
            Swapped = false;
        }

    }

    void updateFirstBead(std::vector<std::vector<int>> &nextConfigsA, std::vector<std::vector<int>> &nextConfigsB,
                         std::vector<std::vector<int>> &nextConfigsAPrime, std::vector<std::vector<int>> &nextConfigsBPrime){

        // p = 0, n = 0
        int index_bead_plus_A = nextConfigsA.at(1).at(0);
        int index_rotor_plus_A = nextConfigsA.at(0).at(1);
        int index_bead_plus_rotor_plus_A = nextConfigsA.at(1).at(1);

        int index_bead_plus_APrime = nextConfigsAPrime.at(1).at(0);
        int index_rotor_plus_APrime = nextConfigsAPrime.at(0).at(1);
        int index_bead_plus_rotor_plus_APrime = nextConfigsAPrime.at(1).at(1);

        std::vector<int> generated_samples = generateSampleFirstBeadFirstRotor(index_bead_plus_A,
                           index_rotor_plus_A,index_bead_plus_rotor_plus_A,
                           index_bead_plus_APrime, index_rotor_plus_APrime,
                           index_bead_plus_rotor_plus_APrime);

        nextConfigsA.at(0).at(0) = generated_samples.at(0);
        nextConfigsAPrime.at(0).at(0) = generated_samples.at(1);

        // p = 0, 0 < n < N_A - 1
        for (int j = 1; j < numRotorsSectorA-1; j++) {
            index_bead_plus_A = nextConfigsA.at(1).at(j);
            index_rotor_plus_A = nextConfigsA.at(0).at(j + 1);
            int index_rotor_minus_A = nextConfigsA.at(0).at(j - 1);
            index_bead_plus_rotor_plus_A = nextConfigsA.at(1).at(j + 1);
            int index_bead_plus_rotor_minus_A = nextConfigsA.at(1).at(j - 1);

            index_bead_plus_APrime = nextConfigsAPrime.at(1).at(j);
            index_rotor_plus_APrime = nextConfigsAPrime.at(0).at(j + 1);
            int index_rotor_minus_APrime = nextConfigsAPrime.at(0).at(j - 1);
            index_bead_plus_rotor_plus_APrime = nextConfigsAPrime.at(1).at(j + 1);
            int index_bead_plus_rotor_minus_APrime = nextConfigsAPrime.at(1).at(j - 1);

            generated_samples.clear();
            generated_samples = generateSampleFirstBeadMidRotor(index_bead_plus_A,
                index_rotor_plus_A,index_rotor_minus_A,
                index_bead_plus_rotor_plus_A,index_bead_plus_rotor_minus_A,
                index_bead_plus_APrime,index_rotor_plus_APrime,
                index_rotor_minus_APrime,index_bead_plus_rotor_plus_APrime,
                index_bead_plus_rotor_minus_APrime);

            nextConfigsA.at(0).at(j) = generated_samples.at(0);
            nextConfigsAPrime.at(0).at(j) = generated_samples.at(1);
        }

        // p = 0, n = N_A - 1
        index_bead_plus_A = nextConfigsA.at(1).at(numRotorsSectorA - 1);
        index_rotor_plus_A = nextConfigsB.at(0).at(0);
        int index_rotor_minus_A = nextConfigsA.at(0).at(numRotorsSectorA - 2);
        index_bead_plus_rotor_plus_A = nextConfigsB.at(1).at(0);
        int index_bead_plus_rotor_minus_A = nextConfigsA.at(1).at(numRotorsSectorA - 2);

        index_bead_plus_APrime = nextConfigsAPrime.at(1).at(numRotorsSectorA - 1);
        index_rotor_plus_APrime = nextConfigsBPrime.at(0).at(0);
        int index_rotor_minus_APrime = nextConfigsAPrime.at(0).at(numRotorsSectorA - 2);
        index_bead_plus_rotor_plus_APrime = nextConfigsBPrime.at(1).at(0);
        int index_bead_plus_rotor_minus_APrime = nextConfigsAPrime.at(1).at(numRotorsSectorA - 2);

        generated_samples.clear();
        generated_samples = generateSampleFirstBeadMidRotor(index_bead_plus_A,index_rotor_plus_A,
                        index_rotor_minus_A, index_bead_plus_rotor_plus_A,
                        index_bead_plus_rotor_minus_A, index_bead_plus_APrime,
                        index_rotor_plus_APrime, index_rotor_minus_APrime,
                        index_bead_plus_rotor_plus_APrime,
                        index_bead_plus_rotor_minus_APrime);

        nextConfigsA.at(0).at(numRotorsSectorA - 1) = generated_samples.at(0);
        nextConfigsAPrime.at(0).at(numRotorsSectorA - 1) = generated_samples.at(1);

        // p = 0, n = N_A (first rotor in sector B)
        int index_bead_plus_B = nextConfigsB.at(1).at(0);
        int index_rotor_plus_B = nextConfigsB.at(0).at(1);
        int index_rotor_minus_B = nextConfigsA.at(0).at(numRotorsSectorA - 1);
        int index_bead_plus_rotor_plus_B = nextConfigsB.at(1).at(1);
        int index_bead_plus_rotor_minus_B = nextConfigsA.at(1).at(numRotorsSectorA - 1);

        int index_bead_plus_BPrime = nextConfigsBPrime.at(1).at(0);
        int index_rotor_plus_BPrime = nextConfigsBPrime.at(0).at(1);
        int index_rotor_minus_BPrime = nextConfigsAPrime.at(0).at(numRotorsSectorA - 1);
        int index_bead_plus_rotor_plus_BPrime = nextConfigsBPrime.at(1).at(1);
        int index_bead_plus_rotor_minus_BPrime = nextConfigsAPrime.at(1).at(numRotorsSectorA - 1);

        generated_samples.clear();
        generated_samples = generateSampleFirstBeadMidRotor(index_bead_plus_B,index_rotor_plus_B,
                            index_rotor_minus_B, index_bead_plus_rotor_plus_B,
                            index_bead_plus_rotor_minus_B, index_bead_plus_BPrime,
                            index_rotor_plus_BPrime, index_rotor_minus_BPrime,
                            index_bead_plus_rotor_plus_BPrime,
                            index_bead_plus_rotor_minus_BPrime);

        nextConfigsB.at(0).at(0) = generated_samples.at(0);
        nextConfigsBPrime.at(0).at(0) = generated_samples.at(1);

        // p = 0, N_A < n < N_B - 1 = 2 * N_A - 1
        for (int j = 1; j < numRotorsSectorB-1; j++) {
            index_bead_plus_B = nextConfigsB.at(1).at(j);
            index_rotor_plus_B = nextConfigsB.at(0).at(j + 1);
            index_rotor_minus_B = nextConfigsB.at(0).at(j - 1);
            index_bead_plus_rotor_plus_B = nextConfigsB.at(1).at(j + 1);
            index_bead_plus_rotor_minus_B = nextConfigsB.at(1).at(j - 1);

            index_bead_plus_BPrime = nextConfigsBPrime.at(1).at(j);
            index_rotor_plus_BPrime = nextConfigsBPrime.at(0).at(j + 1);
            index_rotor_minus_BPrime = nextConfigsBPrime.at(0).at(j - 1);
            index_bead_plus_rotor_plus_BPrime = nextConfigsBPrime.at(1).at(j + 1);
            index_bead_plus_rotor_minus_BPrime = nextConfigsBPrime.at(1).at(j - 1);

            generated_samples.clear();
            generated_samples = generateSampleFirstBeadMidRotor(index_bead_plus_B, index_rotor_plus_B,
                        index_rotor_minus_B, index_bead_plus_rotor_plus_B,
                        index_bead_plus_rotor_minus_B, index_bead_plus_BPrime,
                        index_rotor_plus_BPrime, index_rotor_minus_BPrime,
                        index_bead_plus_rotor_plus_BPrime,
                        index_bead_plus_rotor_minus_BPrime);

            nextConfigsB.at(0).at(j) = generated_samples.at(0);
            nextConfigsBPrime.at(0).at(j) = generated_samples.at(1);
        }

        // p = 0, n = N_B - 1 = 2 * N_A - 1
        index_bead_plus_B = nextConfigsB.at(1).at(numRotorsSectorB - 1);
        index_rotor_minus_B = nextConfigsB.at(0).at(numRotorsSectorB - 2);
        index_bead_plus_rotor_minus_B = nextConfigsB.at(1).at(numRotorsSectorB - 2);

        index_bead_plus_BPrime = nextConfigsBPrime.at(1).at(numRotorsSectorB - 1);
        index_rotor_minus_BPrime = nextConfigsBPrime.at(0).at(numRotorsSectorB - 2);
        index_bead_plus_rotor_minus_BPrime = nextConfigsBPrime.at(1).at(numRotorsSectorB - 2);

        generated_samples.clear();
        generated_samples = generateSampleFirstBeadLastRotor(index_bead_plus_B,index_rotor_minus_B,
                 index_bead_plus_rotor_minus_B, index_bead_plus_BPrime,
                 index_rotor_minus_BPrime,index_bead_plus_rotor_minus_BPrime);

        nextConfigsB.at(0).at(numRotorsSectorB - 1) = generated_samples.at(0);
        nextConfigsBPrime.at(0).at(numRotorsSectorB - 1) = generated_samples.at(1);

    }

    void updateFinalBead(std::vector<std::vector<int>> &nextConfigsA, std::vector<std::vector<int>> &nextConfigsB,
                     std::vector<std::vector<int>> &nextConfigsAPrime, std::vector<std::vector<int>> &nextConfigsBPrime){

        // p = P, n = 0
        int index_bead_minus_A = nextConfigsA.at(numBeads-1).at(0);
        int index_rotor_plus_A = nextConfigsA.at(numBeads).at(1);
        int index_bead_minus_rotor_plus_A = nextConfigsA.at(numBeads-1).at(1);

        int index_bead_minus_APrime = nextConfigsAPrime.at(numBeads-1).at(0);
        int index_rotor_plus_APrime = nextConfigsAPrime.at(numBeads).at(1);
        int index_bead_minus_rotor_plus_APrime = nextConfigsAPrime.at(numBeads-1).at(1);

        std::vector<int> generated_samples = generateSampleLastBeadFirstRotor(index_bead_minus_A,
          index_rotor_plus_A,index_bead_minus_rotor_plus_A,
          index_bead_minus_APrime, index_rotor_plus_APrime,
          index_bead_minus_rotor_plus_APrime);

        nextConfigsA.at(numBeads).at(0) = generated_samples.at(0);
        nextConfigsAPrime.at(numBeads).at(0) = generated_samples.at(1);

        // p = 0, 0 < n < N_A - 1
        for (int j = 1; j < numRotorsSectorA-1; j++) {
            index_bead_minus_A = nextConfigsA.at(numBeads-1).at(j);
            index_rotor_plus_A = nextConfigsA.at(numBeads).at(j + 1);
            int index_rotor_minus_A = nextConfigsA.at(numBeads).at(j - 1);
            index_bead_minus_rotor_plus_A = nextConfigsA.at(numBeads-1).at(j + 1);
            int index_bead_minus_rotor_minus_A = nextConfigsA.at(numBeads-1).at(j - 1);

            index_bead_minus_APrime = nextConfigsAPrime.at(numBeads-1).at(j);
            index_rotor_plus_APrime = nextConfigsAPrime.at(numBeads).at(j + 1);
            int index_rotor_minus_APrime = nextConfigsAPrime.at(numBeads).at(j - 1);
            index_bead_minus_rotor_plus_APrime = nextConfigsAPrime.at(numBeads-1).at(j + 1);
            int index_bead_minus_rotor_minus_APrime = nextConfigsAPrime.at(numBeads-1).at(j - 1);

            generated_samples.clear();
            generated_samples = generateSampleLastBeadMidRotor(index_rotor_plus_A,
                index_rotor_minus_A,index_bead_minus_A,
                index_bead_minus_rotor_plus_A,index_bead_minus_rotor_minus_A,
                index_rotor_plus_APrime,index_rotor_minus_APrime,
                index_bead_minus_APrime,index_bead_minus_rotor_plus_APrime,
                index_bead_minus_rotor_minus_APrime);

            nextConfigsA.at(numBeads).at(j) = generated_samples.at(0);
            nextConfigsAPrime.at(numBeads).at(j) = generated_samples.at(1);
        }

        // p = 0, n = N_A - 1
        index_bead_minus_A = nextConfigsA.at(numBeads-1).at(numRotorsSectorA-1);
        index_rotor_plus_A = nextConfigsB.at(numBeads).at(0);
        int index_rotor_minus_A = nextConfigsA.at(numBeads).at(numRotorsSectorA-2);
        index_bead_minus_rotor_plus_A = nextConfigsB.at(numBeads-1).at(0);
        int index_bead_minus_rotor_minus_A = nextConfigsA.at(numBeads-1).at(numRotorsSectorA-2);

        index_bead_minus_APrime = nextConfigsAPrime.at(numBeads-1).at(numRotorsSectorA-1);
        index_rotor_plus_APrime = nextConfigsBPrime.at(numBeads).at(0);
        int index_rotor_minus_APrime = nextConfigsAPrime.at(numBeads).at(numRotorsSectorA-2);
        index_bead_minus_rotor_plus_APrime = nextConfigsBPrime.at(numBeads-1).at(0);
        int index_bead_minus_rotor_minus_APrime = nextConfigsAPrime.at(numBeads-1).at(numRotorsSectorA-2);

        generated_samples.clear();
        generated_samples = generateSampleLastBeadMidRotor(index_rotor_plus_A,
           index_rotor_minus_A,index_bead_minus_A,
           index_bead_minus_rotor_plus_A,index_bead_minus_rotor_minus_A,
           index_rotor_plus_APrime,index_rotor_minus_APrime,
           index_bead_minus_APrime,index_bead_minus_rotor_plus_APrime,
           index_bead_minus_rotor_minus_APrime);

        nextConfigsA.at(numBeads).at(numRotorsSectorA-1) = generated_samples.at(0);
        nextConfigsAPrime.at(numBeads).at(numRotorsSectorA-1) = generated_samples.at(1);

        // p = 0, n = N_A (first rotor in sector B)
        int index_bead_minus_B = nextConfigsB.at(numBeads-1).at(0);
        int index_rotor_plus_B = nextConfigsB.at(numBeads).at(1);
        int index_rotor_minus_B = nextConfigsA.at(numBeads).at(numRotorsSectorA-1);
        int index_bead_minus_rotor_plus_B = nextConfigsB.at(numBeads-1).at(1);
        int index_bead_minus_rotor_minus_B = nextConfigsA.at(numBeads-1).at(numRotorsSectorA-1);

        int index_bead_minus_BPrime = nextConfigsBPrime.at(numBeads-1).at(0);
        int index_rotor_plus_BPrime = nextConfigsBPrime.at(numBeads).at(1);
        int index_rotor_minus_BPrime = nextConfigsAPrime.at(numBeads).at(numRotorsSectorA-1);
        int index_bead_minus_rotor_plus_BPrime = nextConfigsBPrime.at(numBeads-1).at(1);
        int index_bead_minus_rotor_minus_BPrime = nextConfigsAPrime.at(numBeads-1).at(numRotorsSectorA-1);

        generated_samples.clear();
        generated_samples = generateSampleLastBeadMidRotor(index_rotor_plus_B,
           index_rotor_minus_B,index_bead_minus_B,
           index_bead_minus_rotor_plus_B,index_bead_minus_rotor_minus_B,
           index_rotor_plus_BPrime,index_rotor_minus_BPrime,
           index_bead_minus_BPrime,index_bead_minus_rotor_plus_BPrime,
           index_bead_minus_rotor_minus_BPrime);

        nextConfigsB.at(numBeads).at(0) = generated_samples.at(0);
        nextConfigsBPrime.at(numBeads).at(0) = generated_samples.at(1);

        // p = 0, N_A < n < N_B - 1 = 2 * N_A - 1
        for (int j = 1; j < numRotorsSectorB-1; j++) {
            index_bead_minus_B = nextConfigsB.at(numBeads-1).at(j);
            index_rotor_plus_B = nextConfigsB.at(numBeads).at(j + 1);
            index_rotor_minus_B = nextConfigsB.at(numBeads).at(j - 1);
            index_bead_minus_rotor_plus_B = nextConfigsB.at(numBeads-1).at(j + 1);
            index_bead_minus_rotor_minus_B = nextConfigsB.at(numBeads-1).at(j - 1);

            index_bead_minus_BPrime = nextConfigsBPrime.at(numBeads-1).at(j);
            index_rotor_plus_BPrime = nextConfigsBPrime.at(numBeads).at(j + 1);
            index_rotor_minus_BPrime = nextConfigsBPrime.at(numBeads).at(j - 1);
            index_bead_minus_rotor_plus_BPrime = nextConfigsBPrime.at(numBeads-1).at(j + 1);
            index_bead_minus_rotor_minus_BPrime = nextConfigsBPrime.at(numBeads-1).at(j - 1);

            generated_samples.clear();
            generated_samples = generateSampleLastBeadMidRotor(index_rotor_plus_B,
               index_rotor_minus_B,index_bead_minus_B,
               index_bead_minus_rotor_plus_B,index_bead_minus_rotor_minus_B,
               index_rotor_plus_BPrime,index_rotor_minus_BPrime,
               index_bead_minus_BPrime,index_bead_minus_rotor_plus_BPrime,
               index_bead_minus_rotor_minus_BPrime);

            nextConfigsB.at(numBeads).at(j) = generated_samples.at(0);
            nextConfigsBPrime.at(numBeads).at(j) = generated_samples.at(1);
        }

        // p = 0, n = N_B - 1 = 2 * N_A - 1
        index_bead_minus_B = nextConfigsB.at(numBeads-1).at(numRotorsSectorB-1);
        index_rotor_minus_B = nextConfigsB.at(numBeads).at(numRotorsSectorB-2);
        index_bead_minus_rotor_minus_B = nextConfigsB.at(numBeads-1).at(numRotorsSectorB-2);

        index_bead_minus_BPrime = nextConfigsBPrime.at(numBeads-1).at(numRotorsSectorB-1);
        index_rotor_minus_BPrime = nextConfigsBPrime.at(numBeads).at(numRotorsSectorB-2);
        index_bead_minus_rotor_minus_BPrime = nextConfigsBPrime.at(numBeads-1).at(numRotorsSectorB-2);

        generated_samples.clear();
        generated_samples = generateSampleLastBeadLastRotor(index_rotor_minus_B, index_bead_minus_B,
        index_bead_minus_rotor_minus_B, index_rotor_minus_BPrime,
        index_bead_minus_BPrime,index_bead_minus_rotor_minus_BPrime);

        nextConfigsB.at(numBeads).at(numRotorsSectorB-1) = generated_samples.at(0);
        nextConfigsBPrime.at(numBeads).at(numRotorsSectorB-1) = generated_samples.at(1);

    }

    void updateMidBeadUnswapped(std::vector<std::vector<int>> &nextConfigsA, std::vector<std::vector<int>> &nextConfigsB,
                                std::vector<std::vector<int>> &nextConfigsAPrime, std::vector<std::vector<int>> &nextConfigsBPrime,
                                const int &beadNum){

        // n = 0
        int index_bead_plus_A = nextConfigsA.at(beadNum+1).at(0);
        int index_bead_minus_A = nextConfigsA.at(beadNum-1).at(0);
        int index_rotor_plus_A = nextConfigsA.at(beadNum).at(1);
        int index_bead_plus_rotor_plus_A = nextConfigsA.at(beadNum+1).at(1);
        int index_bead_minus_rotor_plus_A = nextConfigsA.at(beadNum-1).at(1);
        int index_bead_plus_APrime = nextConfigsAPrime.at(beadNum+1).at(0);
        int index_bead_minus_APrime = nextConfigsAPrime.at(beadNum-1).at(0);
        int index_rotor_plus_APrime = nextConfigsAPrime.at(beadNum).at(1);
        int index_bead_plus_rotor_plus_APrime = nextConfigsAPrime.at(beadNum+1).at(1);
        int index_bead_minus_rotor_plus_APrime = nextConfigsAPrime.at(beadNum-1).at(1);

        std::vector<int> generated_samples = generateSampleMidBeadFirstRotor(index_bead_plus_A,
         index_bead_minus_A,index_rotor_plus_A,
         index_bead_plus_rotor_plus_A,index_bead_minus_rotor_plus_A,
         index_bead_plus_APrime, index_bead_minus_APrime,
         index_rotor_plus_APrime, index_bead_plus_rotor_plus_APrime,
         index_bead_minus_rotor_plus_APrime);

        nextConfigsA.at(beadNum).at(0) = generated_samples.at(0);
        nextConfigsAPrime.at(beadNum).at(0) = generated_samples.at(1);

        // 1 <= n <= N_A-2
        for (int j = 1; j < numRotorsSectorA-1; j++) {
            index_bead_plus_A = nextConfigsA.at(beadNum+1).at(j);
            index_bead_minus_A = nextConfigsA.at(beadNum-1).at(j);
            index_rotor_plus_A = nextConfigsA.at(beadNum).at(j+1);
            int index_rotor_minus_A = nextConfigsA.at(beadNum).at(j-1);
            index_bead_plus_rotor_plus_A = nextConfigsA.at(beadNum+1).at(j+1);
            int index_bead_plus_rotor_minus_A = nextConfigsA.at(beadNum+1).at(j-1);
            index_bead_minus_rotor_plus_A = nextConfigsA.at(beadNum-1).at(j+1);
            int index_bead_minus_rotor_minus_A = nextConfigsA.at(beadNum-1).at(j-1);
            index_bead_plus_APrime = nextConfigsAPrime.at(beadNum+1).at(j);
            index_bead_minus_APrime = nextConfigsAPrime.at(beadNum-1).at(j);
            index_rotor_plus_APrime = nextConfigsAPrime.at(beadNum).at(j+1);
            int index_rotor_minus_APrime = nextConfigsAPrime.at(beadNum).at(j-1);
            index_bead_plus_rotor_plus_APrime = nextConfigsAPrime.at(beadNum+1).at(j+1);
            int index_bead_plus_rotor_minus_APrime = nextConfigsAPrime.at(beadNum+1).at(j-1);
            index_bead_minus_rotor_plus_APrime = nextConfigsAPrime.at(beadNum-1).at(j+1);
            int index_bead_minus_rotor_minus_APrime = nextConfigsAPrime.at(beadNum-1).at(j-1);

            generated_samples.clear();
            generated_samples = generateSampleMidBeadMidRotor(index_rotor_plus_A, index_rotor_minus_A,
            index_bead_plus_A,index_bead_plus_rotor_plus_A, index_bead_plus_rotor_minus_A,
            index_bead_minus_A,index_bead_minus_rotor_plus_A, index_bead_minus_rotor_minus_A,
            index_rotor_plus_APrime,index_rotor_minus_APrime, index_bead_plus_APrime,
            index_bead_plus_rotor_plus_APrime,index_bead_plus_rotor_minus_APrime,
            index_bead_minus_APrime,index_bead_minus_rotor_plus_APrime,
            index_bead_minus_rotor_minus_APrime);

            nextConfigsA.at(beadNum).at(j) = generated_samples.at(0);
            nextConfigsAPrime.at(beadNum).at(j) = generated_samples.at(1);
        }

        // n = N_A - 1
        index_bead_plus_A = nextConfigsA.at(beadNum+1).at(numRotorsSectorA-1);
        index_bead_minus_A = nextConfigsA.at(beadNum-1).at(numRotorsSectorA-1);
        index_rotor_plus_A = nextConfigsB.at(beadNum).at(0);
        int index_rotor_minus_A = nextConfigsA.at(beadNum).at(numRotorsSectorA-2);
        index_bead_plus_rotor_plus_A = nextConfigsB.at(beadNum+1).at(0);
        int index_bead_plus_rotor_minus_A = nextConfigsA.at(beadNum+1).at(numRotorsSectorA-2);
        index_bead_minus_rotor_plus_A = nextConfigsB.at(beadNum-1).at(0);
        int index_bead_minus_rotor_minus_A = nextConfigsA.at(beadNum-1).at(numRotorsSectorA-2);
        index_bead_plus_APrime = nextConfigsAPrime.at(beadNum+1).at(numRotorsSectorA-1);
        index_bead_minus_APrime = nextConfigsAPrime.at(beadNum-1).at(numRotorsSectorA-1);
        index_rotor_plus_APrime = nextConfigsBPrime.at(beadNum).at(0);
        int index_rotor_minus_APrime = nextConfigsAPrime.at(beadNum).at(numRotorsSectorA-2);
        index_bead_plus_rotor_plus_APrime = nextConfigsBPrime.at(beadNum+1).at(0);
        int index_bead_plus_rotor_minus_APrime = nextConfigsAPrime.at(beadNum+1).at(numRotorsSectorA-2);
        index_bead_minus_rotor_plus_APrime = nextConfigsBPrime.at(beadNum-1).at(0);
        int index_bead_minus_rotor_minus_APrime = nextConfigsAPrime.at(beadNum-1).at(numRotorsSectorA-2);

        generated_samples.clear();
        generated_samples = generateSampleMidBeadMidRotor(index_rotor_plus_A, index_rotor_minus_A,
                          index_bead_plus_A,index_bead_plus_rotor_plus_A,
                          index_bead_plus_rotor_minus_A,index_bead_minus_A,
                          index_bead_minus_rotor_plus_A,index_bead_minus_rotor_minus_A,
                          index_rotor_plus_APrime,index_rotor_minus_APrime,
                          index_bead_plus_APrime,index_bead_plus_rotor_plus_APrime,
                          index_bead_plus_rotor_minus_APrime,index_bead_minus_APrime,
                          index_bead_minus_rotor_plus_APrime,index_bead_minus_rotor_minus_APrime);

        nextConfigsA.at(beadNum).at(numRotorsSectorA - 1) = generated_samples.at(0);
        nextConfigsAPrime.at(beadNum).at(numRotorsSectorA - 1) = generated_samples.at(1);

        // n = N_A (first rotor in sector B)
        int index_bead_plus_B = nextConfigsB.at(beadNum+1).at(0);
        int index_bead_minus_B = nextConfigsB.at(beadNum-1).at(0);
        int index_rotor_plus_B = nextConfigsB.at(beadNum).at(1);
        int index_rotor_minus_B = nextConfigsA.at(beadNum).at(numRotorsSectorA-1);
        int index_bead_plus_rotor_plus_B = nextConfigsB.at(beadNum+1).at(1);
        int index_bead_plus_rotor_minus_B = nextConfigsA.at(beadNum+1).at(numRotorsSectorA-1);
        int index_bead_minus_rotor_plus_B = nextConfigsB.at(beadNum-1).at(1);
        int index_bead_minus_rotor_minus_B = nextConfigsA.at(beadNum-1).at(numRotorsSectorA-1);
        int index_bead_plus_BPrime = nextConfigsBPrime.at(beadNum+1).at(0);
        int index_bead_minus_BPrime = nextConfigsBPrime.at(beadNum-1).at(0);
        int index_rotor_plus_BPrime = nextConfigsBPrime.at(beadNum).at(1);
        int index_rotor_minus_BPrime = nextConfigsAPrime.at(beadNum).at(numRotorsSectorA-1);
        int index_bead_plus_rotor_plus_BPrime = nextConfigsBPrime.at(beadNum+1).at(1);
        int index_bead_plus_rotor_minus_BPrime = nextConfigsAPrime.at(beadNum+1).at(numRotorsSectorA-1);
        int index_bead_minus_rotor_plus_BPrime = nextConfigsBPrime.at(beadNum-1).at(1);
        int index_bead_minus_rotor_minus_BPrime = nextConfigsAPrime.at(beadNum-1).at(numRotorsSectorA-1);

        generated_samples.clear();
        generated_samples = generateSampleMidBeadMidRotor(index_rotor_plus_B, index_rotor_minus_B,
          index_bead_plus_B,index_bead_plus_rotor_plus_B,
          index_bead_plus_rotor_minus_B,index_bead_minus_B,
          index_bead_minus_rotor_plus_B,index_bead_minus_rotor_minus_B,
          index_rotor_plus_BPrime,index_rotor_minus_BPrime,
          index_bead_plus_BPrime,index_bead_plus_rotor_plus_BPrime,
          index_bead_plus_rotor_minus_BPrime,index_bead_minus_BPrime,
          index_bead_minus_rotor_plus_BPrime,index_bead_minus_rotor_minus_BPrime);

        nextConfigsB.at(beadNum).at(0) = generated_samples.at(0);
        nextConfigsBPrime.at(beadNum).at(0) = generated_samples.at(1);

        // N_A < n < N_B - 1 = 2 * N_A - 1
        for (int j = 1; j < numRotorsSectorB-1; j++) {
            index_bead_plus_B = nextConfigsB.at(beadNum+1).at(j);
            index_bead_minus_B = nextConfigsB.at(beadNum-1).at(j);
            index_rotor_plus_B = nextConfigsB.at(beadNum).at(j+1);
            index_rotor_minus_B = nextConfigsB.at(beadNum).at(j-1);
            index_bead_plus_rotor_plus_B = nextConfigsB.at(beadNum+1).at(j+1);
            index_bead_plus_rotor_minus_B = nextConfigsB.at(beadNum+1).at(j-1);
            index_bead_minus_rotor_plus_B = nextConfigsB.at(beadNum-1).at(j+1);
            index_bead_minus_rotor_minus_B = nextConfigsB.at(beadNum-1).at(j-1);
            index_bead_plus_BPrime = nextConfigsBPrime.at(beadNum+1).at(j);
            index_bead_minus_BPrime = nextConfigsBPrime.at(beadNum-1).at(j);
            index_rotor_plus_BPrime = nextConfigsBPrime.at(beadNum).at(j+1);
            index_rotor_minus_BPrime = nextConfigsBPrime.at(beadNum).at(j-1);
            index_bead_plus_rotor_plus_BPrime = nextConfigsBPrime.at(beadNum+1).at(j+1);
            index_bead_plus_rotor_minus_BPrime = nextConfigsBPrime.at(beadNum+1).at(j-1);
            index_bead_minus_rotor_plus_BPrime = nextConfigsBPrime.at(beadNum-1).at(j+1);
            index_bead_minus_rotor_minus_BPrime = nextConfigsBPrime.at(beadNum-1).at(j-1);

            generated_samples.clear();
            generated_samples = generateSampleMidBeadMidRotor(index_rotor_plus_B, index_rotor_minus_B,
              index_bead_plus_B,index_bead_plus_rotor_plus_B, index_bead_plus_rotor_minus_B,
              index_bead_minus_B,index_bead_minus_rotor_plus_B, index_bead_minus_rotor_minus_B,
              index_rotor_plus_BPrime,index_rotor_minus_BPrime, index_bead_plus_BPrime,
              index_bead_plus_rotor_plus_BPrime,index_bead_plus_rotor_minus_BPrime,
              index_bead_minus_BPrime,index_bead_minus_rotor_plus_BPrime,
              index_bead_minus_rotor_minus_BPrime);

            nextConfigsB.at(beadNum).at(j) = generated_samples.at(0);
            nextConfigsBPrime.at(beadNum).at(j) = generated_samples.at(1);
        }

        // n = N_B - 1 = 2*N_A - 1
        index_bead_plus_B = nextConfigsB.at(beadNum+1).at(numRotorsSectorB - 1);
        index_bead_minus_B = nextConfigsB.at(beadNum-1).at(numRotorsSectorB - 1);
        index_rotor_minus_B = nextConfigsB.at(beadNum).at(numRotorsSectorB - 2);
        index_bead_plus_rotor_minus_B = nextConfigsB.at(beadNum+1).at(numRotorsSectorB - 2);
        index_bead_minus_rotor_minus_B = nextConfigsB.at(beadNum-1).at(numRotorsSectorB - 2);
        index_bead_plus_BPrime = nextConfigsBPrime.at(beadNum+1).at(numRotorsSectorB - 1);
        index_bead_minus_BPrime = nextConfigsBPrime.at(beadNum-1).at(numRotorsSectorB - 1);
        index_rotor_minus_BPrime = nextConfigsBPrime.at(beadNum).at(numRotorsSectorB - 2);
        index_bead_plus_rotor_minus_BPrime = nextConfigsBPrime.at(beadNum+1).at(numRotorsSectorB - 2);
        index_bead_minus_rotor_minus_BPrime = nextConfigsBPrime.at(beadNum-1).at(numRotorsSectorB - 2);

        generated_samples.clear();
        generated_samples = generateSampleMidBeadLastRotor(index_bead_plus_B, index_rotor_minus_B,
                           index_bead_plus_rotor_minus_B, index_bead_minus_B,
                           index_bead_minus_rotor_minus_B,index_bead_plus_BPrime,
                           index_rotor_minus_BPrime,index_bead_plus_rotor_minus_BPrime,
                           index_bead_minus_BPrime,index_bead_minus_rotor_minus_BPrime);

        nextConfigsB.at(beadNum).at(numRotorsSectorB - 1) = generated_samples.at(0);
        nextConfigsBPrime.at(beadNum).at(numRotorsSectorB - 1) = generated_samples.at(1);

    }

    void updateMidBeadsSwapped(std::vector<std::vector<int>> &nextConfigsA, std::vector<std::vector<int>> &nextConfigsB,
                                std::vector<std::vector<int>> &nextConfigsAPrime, std::vector<std::vector<int>> &nextConfigsBPrime){

        // p = P/2, n = 0
        int beadNum = numBeads/2;
        int index_bead_plus_A = nextConfigsAPrime.at(beadNum+1).at(0);
        int index_bead_plus_rotor_plus_A = nextConfigsAPrime.at(beadNum+1).at(1);

        int index_rotor_plus_A = nextConfigsA.at(beadNum).at(1);

        int index_bead_minus_A = nextConfigsA.at(beadNum-1).at(0);
        int index_bead_minus_rotor_plus_A = nextConfigsA.at(beadNum-1).at(1);

        int index_bead_plus_APrime = nextConfigsA.at(beadNum+1).at(0);
        int index_bead_plus_rotor_plus_APrime = nextConfigsA.at(beadNum+1).at(1);

        int index_rotor_plus_APrime = nextConfigsAPrime.at(beadNum).at(1);

        int index_bead_minus_APrime = nextConfigsAPrime.at(beadNum-1).at(0);
        int index_bead_minus_rotor_plus_APrime = nextConfigsAPrime.at(beadNum-1).at(1);

        std::vector<int> generated_samples = generateSampleMidBeadFirstRotor(index_bead_plus_A,
         index_bead_minus_A,index_rotor_plus_A,
         index_bead_plus_rotor_plus_A,index_bead_minus_rotor_plus_A,
         index_bead_plus_APrime, index_bead_minus_APrime,
         index_rotor_plus_APrime, index_bead_plus_rotor_plus_APrime,
         index_bead_minus_rotor_plus_APrime);

        nextConfigsA.at(beadNum).at(0) = generated_samples.at(0);
        nextConfigsAPrime.at(beadNum).at(0) = generated_samples.at(1);

        // 1 <= n <= N_A-2
        for (int j = 1; j < numRotorsSectorA-1; j++) {
            index_bead_plus_A = nextConfigsAPrime.at(beadNum+1).at(j);
            int index_bead_plus_rotor_minus_A = nextConfigsAPrime.at(beadNum+1).at(j-1);
            index_bead_plus_rotor_plus_A = nextConfigsAPrime.at(beadNum+1).at(j+1);

            index_rotor_plus_A = nextConfigsA.at(beadNum).at(j+1);
            int index_rotor_minus_A = nextConfigsA.at(beadNum).at(j-1);

            index_bead_minus_A = nextConfigsA.at(beadNum-1).at(j);
            index_bead_minus_rotor_plus_A = nextConfigsA.at(beadNum-1).at(j+1);
            int index_bead_minus_rotor_minus_A = nextConfigsA.at(beadNum-1).at(j-1);

            index_bead_plus_APrime = nextConfigsA.at(beadNum+1).at(j);
            int index_bead_plus_rotor_minus_APrime = nextConfigsA.at(beadNum+1).at(j-1);
            index_bead_plus_rotor_plus_APrime = nextConfigsA.at(beadNum+1).at(j+1);

            index_rotor_plus_APrime = nextConfigsAPrime.at(beadNum).at(j+1);
            int index_rotor_minus_APrime = nextConfigsAPrime.at(beadNum).at(j-1);

            index_bead_minus_APrime = nextConfigsAPrime.at(beadNum-1).at(j);
            index_bead_minus_rotor_plus_APrime = nextConfigsAPrime.at(beadNum-1).at(j+1);
            int index_bead_minus_rotor_minus_APrime = nextConfigsAPrime.at(beadNum-1).at(j-1);

            generated_samples.clear();
            generated_samples = generateSampleMidBeadMidRotor(index_rotor_plus_A, index_rotor_minus_A,
              index_bead_plus_A,index_bead_plus_rotor_plus_A, index_bead_plus_rotor_minus_A,
              index_bead_minus_A,index_bead_minus_rotor_plus_A, index_bead_minus_rotor_minus_A,
              index_rotor_plus_APrime,index_rotor_minus_APrime, index_bead_plus_APrime,
              index_bead_plus_rotor_plus_APrime,index_bead_plus_rotor_minus_APrime,
              index_bead_minus_APrime,index_bead_minus_rotor_plus_APrime,
              index_bead_minus_rotor_minus_APrime);

            nextConfigsA.at(beadNum).at(j) = generated_samples.at(0);
            nextConfigsAPrime.at(beadNum).at(j) = generated_samples.at(1);
        }

        // n = N_A - 1
        index_bead_plus_A = nextConfigsAPrime.at(beadNum+1).at(numRotorsSectorA-1);
        index_bead_plus_rotor_plus_A = nextConfigsBPrime.at(beadNum+1).at(0);
        int index_bead_plus_rotor_minus_A = nextConfigsAPrime.at(beadNum+1).at(numRotorsSectorA-2);

        index_rotor_plus_A = nextConfigsB.at(beadNum).at(0);
        int index_rotor_plus_A_p = nextConfigsBPrime.at(beadNum).at(0);
        int index_rotor_minus_A = nextConfigsA.at(beadNum).at(numRotorsSectorA-2);

        index_bead_minus_A = nextConfigsA.at(beadNum-1).at(numRotorsSectorA-1);
        index_bead_minus_rotor_plus_A = nextConfigsB.at(beadNum-1).at(0);
        int index_bead_minus_rotor_minus_A = nextConfigsA.at(beadNum-1).at(numRotorsSectorA-2);

        index_bead_plus_APrime = nextConfigsA.at(beadNum+1).at(numRotorsSectorA-1);
        index_bead_plus_rotor_plus_APrime = nextConfigsB.at(beadNum+1).at(0);
        int index_bead_plus_rotor_minus_APrime = nextConfigsA.at(beadNum+1).at(numRotorsSectorA-2);

        index_rotor_plus_APrime = nextConfigsBPrime.at(beadNum).at(0);
        int index_rotor_plus_APrime_p = nextConfigsB.at(beadNum).at(0);
        int index_rotor_minus_APrime = nextConfigsAPrime.at(beadNum).at(numRotorsSectorA-2);

        index_bead_minus_APrime = nextConfigsAPrime.at(beadNum-1).at(numRotorsSectorA-1);
        index_bead_minus_rotor_plus_APrime = nextConfigsBPrime.at(beadNum-1).at(0);
        int index_bead_minus_rotor_minus_APrime = nextConfigsAPrime.at(beadNum-1).at(numRotorsSectorA-2);

        generated_samples.clear();
        generated_samples = generateSampleMidBeadMidRotorSwappedA(index_rotor_plus_A,
          index_rotor_plus_A_p, index_rotor_minus_A,
          index_bead_plus_A,index_bead_plus_rotor_plus_A,
          index_bead_plus_rotor_minus_A,index_bead_minus_A,
          index_bead_minus_rotor_plus_A,index_bead_minus_rotor_minus_A,
          index_rotor_plus_APrime,index_rotor_plus_APrime_p,
          index_rotor_minus_APrime,index_bead_plus_APrime,
          index_bead_plus_rotor_plus_APrime,
          index_bead_plus_rotor_minus_APrime,index_bead_minus_APrime,
          index_bead_minus_rotor_plus_APrime,
          index_bead_minus_rotor_minus_APrime);

        nextConfigsA.at(beadNum).at(numRotorsSectorA - 1) = generated_samples.at(0);
        nextConfigsAPrime.at(beadNum).at(numRotorsSectorA - 1) = generated_samples.at(1);

        // n = N_A (first rotor in sector B)
        int index_bead_plus_B = nextConfigsB.at(beadNum+1).at(0);
        int index_bead_plus_rotor_plus_B = nextConfigsB.at(beadNum+1).at(1);
        int index_bead_plus_rotor_minus_B = nextConfigsA.at(beadNum+1).at(numRotorsSectorA-1);

        int index_rotor_plus_B = nextConfigsB.at(beadNum).at(1);
        int index_rotor_minus_B = nextConfigsA.at(beadNum).at(numRotorsSectorA-1);
        int index_rotor_minus_B_p = nextConfigsAPrime.at(beadNum).at(numRotorsSectorA-1);

        int index_bead_minus_B = nextConfigsB.at(beadNum-1).at(0);
        int index_bead_minus_rotor_plus_B = nextConfigsB.at(beadNum-1).at(1);
        int index_bead_minus_rotor_minus_B = nextConfigsA.at(beadNum-1).at(numRotorsSectorA-1);

        int index_bead_plus_BPrime = nextConfigsBPrime.at(beadNum+1).at(0);
        int index_bead_plus_rotor_plus_BPrime = nextConfigsBPrime.at(beadNum+1).at(1);
        int index_bead_plus_rotor_minus_BPrime = nextConfigsAPrime.at(beadNum+1).at(numRotorsSectorA-1);

        int index_rotor_plus_BPrime = nextConfigsBPrime.at(beadNum).at(1);
        int index_rotor_minus_BPrime = nextConfigsAPrime.at(beadNum).at(numRotorsSectorA-1);
        int index_rotor_minus_BPrime_p = nextConfigsA.at(beadNum).at(numRotorsSectorA-1);

        int index_bead_minus_BPrime = nextConfigsBPrime.at(beadNum-1).at(0);
        int index_bead_minus_rotor_plus_BPrime = nextConfigsBPrime.at(beadNum-1).at(1);
        int index_bead_minus_rotor_minus_BPrime = nextConfigsAPrime.at(beadNum-1).at(numRotorsSectorA-1);

        generated_samples.clear();
        generated_samples = generateSampleMidBeadMidRotorSwappedB(index_rotor_plus_B, index_rotor_minus_B,
          index_rotor_minus_B_p,index_bead_plus_B,index_bead_plus_rotor_plus_B,
          index_bead_plus_rotor_minus_B,index_bead_minus_B,
          index_bead_minus_rotor_plus_B,index_bead_minus_rotor_minus_B,
          index_rotor_plus_BPrime, index_rotor_minus_BPrime, index_rotor_minus_BPrime_p,
          index_bead_plus_BPrime,index_bead_plus_rotor_plus_BPrime,
          index_bead_plus_rotor_minus_BPrime,index_bead_minus_BPrime,
          index_bead_minus_rotor_plus_BPrime,index_bead_minus_rotor_minus_BPrime);

        nextConfigsB.at(beadNum).at(0) = generated_samples.at(0);
        nextConfigsBPrime.at(beadNum).at(0) = generated_samples.at(1);

        // N_A < n < N_B - 1 = 2 * N_A - 1
        for (int j = 1; j < numRotorsSectorB-1; j++) {
            index_bead_plus_B = nextConfigsB.at(beadNum+1).at(j);
            index_bead_minus_B = nextConfigsB.at(beadNum-1).at(j);
            index_rotor_plus_B = nextConfigsB.at(beadNum).at(j+1);
            index_rotor_minus_B = nextConfigsB.at(beadNum).at(j-1);
            index_bead_plus_rotor_plus_B = nextConfigsB.at(beadNum+1).at(j+1);
            index_bead_plus_rotor_minus_B = nextConfigsB.at(beadNum+1).at(j-1);
            index_bead_minus_rotor_plus_B = nextConfigsB.at(beadNum-1).at(j+1);
            index_bead_minus_rotor_minus_B = nextConfigsB.at(beadNum-1).at(j-1);

            index_bead_plus_BPrime = nextConfigsBPrime.at(beadNum+1).at(j);
            index_bead_minus_BPrime = nextConfigsBPrime.at(beadNum-1).at(j);
            index_rotor_plus_BPrime = nextConfigsBPrime.at(beadNum).at(j+1);
            index_rotor_minus_BPrime = nextConfigsBPrime.at(beadNum).at(j-1);
            index_bead_plus_rotor_plus_BPrime = nextConfigsBPrime.at(beadNum+1).at(j+1);
            index_bead_plus_rotor_minus_BPrime = nextConfigsBPrime.at(beadNum+1).at(j-1);
            index_bead_minus_rotor_plus_BPrime = nextConfigsBPrime.at(beadNum-1).at(j+1);
            index_bead_minus_rotor_minus_BPrime = nextConfigsBPrime.at(beadNum-1).at(j-1);

            generated_samples.clear();
            generated_samples = generateSampleMidBeadMidRotor(index_rotor_plus_B, index_rotor_minus_B,
              index_bead_plus_B,index_bead_plus_rotor_plus_B, index_bead_plus_rotor_minus_B,
              index_bead_minus_B,index_bead_minus_rotor_plus_B, index_bead_minus_rotor_minus_B,
              index_rotor_plus_BPrime,index_rotor_minus_BPrime, index_bead_plus_BPrime,
              index_bead_plus_rotor_plus_BPrime,index_bead_plus_rotor_minus_BPrime,
              index_bead_minus_BPrime,index_bead_minus_rotor_plus_BPrime,
              index_bead_minus_rotor_minus_BPrime);

            nextConfigsB.at(beadNum).at(j) = generated_samples.at(0);
            nextConfigsBPrime.at(beadNum).at(j) = generated_samples.at(1);
        }

        // n = N_B - 1 = 2*N_A - 1
        index_bead_plus_B = nextConfigsB.at(beadNum+1).at(numRotorsSectorB-1);
        index_bead_minus_B = nextConfigsB.at(beadNum-1).at(numRotorsSectorB-1);
        index_rotor_minus_B = nextConfigsB.at(beadNum).at(numRotorsSectorB-2);
        index_bead_plus_rotor_minus_B = nextConfigsB.at(beadNum+1).at(numRotorsSectorB-2);
        index_bead_minus_rotor_minus_B = nextConfigsB.at(beadNum-1).at(numRotorsSectorB-2);
        
        index_bead_plus_BPrime = nextConfigsBPrime.at(beadNum+1).at(numRotorsSectorB-1);
        index_bead_minus_BPrime = nextConfigsBPrime.at(beadNum-1).at(numRotorsSectorB-1);
        index_rotor_minus_BPrime = nextConfigsBPrime.at(beadNum).at(numRotorsSectorB-2);
        index_bead_plus_rotor_minus_BPrime = nextConfigsBPrime.at(beadNum+1).at(numRotorsSectorB-2);
        index_bead_minus_rotor_minus_BPrime = nextConfigsBPrime.at(beadNum-1).at(numRotorsSectorB-2);

        generated_samples.clear();
        generated_samples = generateSampleMidBeadLastRotor(index_bead_plus_B, index_rotor_minus_B,
           index_bead_plus_rotor_minus_B, index_bead_minus_B,
           index_bead_minus_rotor_minus_B,index_bead_plus_BPrime,
           index_rotor_minus_BPrime,index_bead_plus_rotor_minus_BPrime,
           index_bead_minus_BPrime,index_bead_minus_rotor_minus_BPrime);

        nextConfigsB.at(beadNum).at(numRotorsSectorB-1) = generated_samples.at(0);
        nextConfigsBPrime.at(beadNum).at(numRotorsSectorB-1) = generated_samples.at(1);

        // p = P/2 + 1
        beadNum = numBeads/2 + 1;

        // n = 0
        index_bead_plus_A = nextConfigsA.at(beadNum+1).at(0);
        index_bead_plus_rotor_plus_A = nextConfigsA.at(beadNum+1).at(1);

        index_rotor_plus_A = nextConfigsA.at(beadNum).at(1);

        index_bead_minus_A = nextConfigsAPrime.at(beadNum-1).at(0);
        index_bead_minus_rotor_plus_A = nextConfigsAPrime.at(beadNum-1).at(1);

        index_bead_plus_APrime = nextConfigsAPrime.at(beadNum+1).at(0);
        index_bead_plus_rotor_plus_APrime = nextConfigsAPrime.at(beadNum+1).at(1);

        index_rotor_plus_APrime = nextConfigsAPrime.at(beadNum).at(1);

        index_bead_minus_APrime = nextConfigsA.at(beadNum-1).at(0);
        index_bead_minus_rotor_plus_APrime = nextConfigsA.at(beadNum-1).at(1);

        generated_samples = generateSampleMidBeadFirstRotor(index_bead_plus_A,
             index_bead_minus_A,index_rotor_plus_A,
             index_bead_plus_rotor_plus_A,index_bead_minus_rotor_plus_A,
             index_bead_plus_APrime, index_bead_minus_APrime,
             index_rotor_plus_APrime, index_bead_plus_rotor_plus_APrime,
             index_bead_minus_rotor_plus_APrime);

        nextConfigsA.at(beadNum).at(0) = generated_samples.at(0);
        nextConfigsAPrime.at(beadNum).at(0) = generated_samples.at(1);

        // 1 <= n <= N_A-2
        for (int j = 1; j < numRotorsSectorA-1; j++) {
            index_bead_plus_A = nextConfigsA.at(beadNum+1).at(j);
            index_bead_plus_rotor_minus_A = nextConfigsA.at(beadNum+1).at(j-1);
            index_bead_plus_rotor_plus_A = nextConfigsA.at(beadNum+1).at(j+1);

            index_rotor_plus_A = nextConfigsA.at(beadNum).at(j+1);
            index_rotor_minus_A = nextConfigsA.at(beadNum).at(j-1);

            index_bead_minus_A = nextConfigsAPrime.at(beadNum-1).at(j);
            index_bead_minus_rotor_plus_A = nextConfigsAPrime.at(beadNum-1).at(j+1);
            index_bead_minus_rotor_minus_A = nextConfigsAPrime.at(beadNum-1).at(j-1);

            index_bead_plus_APrime = nextConfigsAPrime.at(beadNum+1).at(j);
            index_bead_plus_rotor_minus_APrime = nextConfigsAPrime.at(beadNum+1).at(j-1);
            index_bead_plus_rotor_plus_APrime = nextConfigsAPrime.at(beadNum+1).at(j+1);

            index_rotor_plus_APrime = nextConfigsAPrime.at(beadNum).at(j+1);
            index_rotor_minus_APrime = nextConfigsAPrime.at(beadNum).at(j-1);

            index_bead_minus_APrime = nextConfigsA.at(beadNum-1).at(j);
            index_bead_minus_rotor_plus_APrime = nextConfigsA.at(beadNum-1).at(j+1);
            index_bead_minus_rotor_minus_APrime = nextConfigsA.at(beadNum-1).at(j-1);

            generated_samples.clear();
            generated_samples = generateSampleMidBeadMidRotor(index_rotor_plus_A, index_rotor_minus_A,
                                                              index_bead_plus_A,index_bead_plus_rotor_plus_A, index_bead_plus_rotor_minus_A,
                                                              index_bead_minus_A,index_bead_minus_rotor_plus_A, index_bead_minus_rotor_minus_A,
                                                              index_rotor_plus_APrime,index_rotor_minus_APrime, index_bead_plus_APrime,
                                                              index_bead_plus_rotor_plus_APrime,index_bead_plus_rotor_minus_APrime,
                                                              index_bead_minus_APrime,index_bead_minus_rotor_plus_APrime,
                                                              index_bead_minus_rotor_minus_APrime);

            nextConfigsA.at(beadNum).at(j) = generated_samples.at(0);
            nextConfigsAPrime.at(beadNum).at(j) = generated_samples.at(1);
        }

        // n = N_A - 1
        index_bead_plus_A = nextConfigsA.at(beadNum+1).at(numRotorsSectorA - 1);
        index_bead_plus_rotor_minus_A = nextConfigsA.at(beadNum+1).at(numRotorsSectorA - 2);
        index_bead_plus_rotor_plus_A = nextConfigsB.at(beadNum+1).at(0);

        index_rotor_plus_A = nextConfigsB.at(beadNum).at(0);
        index_rotor_minus_A = nextConfigsA.at(beadNum).at(numRotorsSectorA - 2);

        index_bead_minus_A = nextConfigsAPrime.at(beadNum-1).at(numRotorsSectorA - 1);
        index_bead_minus_rotor_plus_A = nextConfigsB.at(beadNum-1).at(0);
        index_bead_minus_rotor_minus_A = nextConfigsAPrime.at(beadNum-1).at(numRotorsSectorA - 2);

        index_bead_plus_APrime = nextConfigsAPrime.at(beadNum+1).at(numRotorsSectorA - 1);
        index_bead_plus_rotor_minus_APrime = nextConfigsAPrime.at(beadNum+1).at(numRotorsSectorA - 2);
        index_bead_plus_rotor_plus_APrime = nextConfigsBPrime.at(beadNum+1).at(0);

        index_rotor_plus_APrime = nextConfigsBPrime.at(beadNum).at(0);
        index_rotor_minus_APrime = nextConfigsAPrime.at(beadNum).at(numRotorsSectorA - 2);

        index_bead_minus_APrime = nextConfigsA.at(beadNum-1).at(numRotorsSectorA - 1);
        index_bead_minus_rotor_plus_APrime = nextConfigsBPrime.at(beadNum-1).at(0);
        index_bead_minus_rotor_minus_APrime = nextConfigsA.at(beadNum-1).at(numRotorsSectorA - 2);

        generated_samples.clear();
        generated_samples = generateSampleMidBeadMidRotor(index_rotor_plus_A, index_rotor_minus_A,
          index_bead_plus_A,index_bead_plus_rotor_plus_A, index_bead_plus_rotor_minus_A,
          index_bead_minus_A,index_bead_minus_rotor_plus_A, index_bead_minus_rotor_minus_A,
          index_rotor_plus_APrime,index_rotor_minus_APrime, index_bead_plus_APrime,
          index_bead_plus_rotor_plus_APrime,index_bead_plus_rotor_minus_APrime,
          index_bead_minus_APrime,index_bead_minus_rotor_plus_APrime,
          index_bead_minus_rotor_minus_APrime);

        nextConfigsA.at(beadNum).at(numRotorsSectorA-1) = generated_samples.at(0);
        nextConfigsAPrime.at(beadNum).at(numRotorsSectorA-1) = generated_samples.at(1);

        // n = N_A (first rotor of sector B)
        index_bead_plus_B = nextConfigsB.at(beadNum+1).at(0);
        index_bead_plus_rotor_minus_B = nextConfigsA.at(beadNum+1).at(numRotorsSectorA-1);
        index_bead_plus_rotor_plus_B = nextConfigsB.at(beadNum+1).at(1);

        index_rotor_plus_B = nextConfigsB.at(beadNum).at(1);
        index_rotor_minus_B = nextConfigsA.at(beadNum).at(numRotorsSectorA - 1);

        index_bead_minus_B = nextConfigsB.at(beadNum-1).at(0);
        index_bead_minus_rotor_plus_B = nextConfigsB.at(beadNum-1).at(1);
        index_bead_minus_rotor_minus_B = nextConfigsAPrime.at(beadNum-1).at(numRotorsSectorA - 1);

        index_bead_plus_BPrime = nextConfigsBPrime.at(beadNum+1).at(0);
        index_bead_plus_rotor_minus_BPrime = nextConfigsAPrime.at(beadNum+1).at(numRotorsSectorA-1);
        index_bead_plus_rotor_plus_BPrime = nextConfigsBPrime.at(beadNum+1).at(1);

        index_rotor_plus_BPrime = nextConfigsBPrime.at(beadNum).at(1);
        index_rotor_minus_BPrime = nextConfigsAPrime.at(beadNum).at(numRotorsSectorA - 1);

        index_bead_minus_BPrime = nextConfigsBPrime.at(beadNum-1).at(0);
        index_bead_minus_rotor_plus_BPrime = nextConfigsBPrime.at(beadNum-1).at(1);
        index_bead_minus_rotor_minus_BPrime = nextConfigsA.at(beadNum-1).at(numRotorsSectorA - 1);

        generated_samples.clear();
        generated_samples = generateSampleMidBeadMidRotor(index_rotor_plus_B, index_rotor_minus_B,
          index_bead_plus_B,index_bead_plus_rotor_plus_B, index_bead_plus_rotor_minus_B,
          index_bead_minus_B,index_bead_minus_rotor_plus_B, index_bead_minus_rotor_minus_B,
          index_rotor_plus_BPrime,index_rotor_minus_BPrime, index_bead_plus_BPrime,
          index_bead_plus_rotor_plus_BPrime,index_bead_plus_rotor_minus_BPrime,
          index_bead_minus_BPrime,index_bead_minus_rotor_plus_BPrime,
          index_bead_minus_rotor_minus_BPrime);

        nextConfigsB.at(beadNum).at(0) = generated_samples.at(0);
        nextConfigsBPrime.at(beadNum).at(0) = generated_samples.at(1);

        // N_A < n < N_B - 1 = 2 * N_A - 1
        for (int j = 1; j < numRotorsSectorB-1; j++) {
            index_bead_plus_B = nextConfigsB.at(beadNum+1).at(j);
            index_bead_minus_B = nextConfigsB.at(beadNum-1).at(j);
            index_rotor_plus_B = nextConfigsB.at(beadNum).at(j+1);
            index_rotor_minus_B = nextConfigsB.at(beadNum).at(j-1);
            index_bead_plus_rotor_plus_B = nextConfigsB.at(beadNum+1).at(j+1);
            index_bead_plus_rotor_minus_B = nextConfigsB.at(beadNum+1).at(j-1);
            index_bead_minus_rotor_plus_B = nextConfigsB.at(beadNum-1).at(j+1);
            index_bead_minus_rotor_minus_B = nextConfigsB.at(beadNum-1).at(j-1);
            index_bead_plus_BPrime = nextConfigsBPrime.at(beadNum+1).at(j);
            index_bead_minus_BPrime = nextConfigsBPrime.at(beadNum-1).at(j);
            index_rotor_plus_BPrime = nextConfigsBPrime.at(beadNum).at(j+1);
            index_rotor_minus_BPrime = nextConfigsBPrime.at(beadNum).at(j-1);
            index_bead_plus_rotor_plus_BPrime = nextConfigsBPrime.at(beadNum+1).at(j+1);
            index_bead_plus_rotor_minus_BPrime = nextConfigsBPrime.at(beadNum+1).at(j-1);
            index_bead_minus_rotor_plus_BPrime = nextConfigsBPrime.at(beadNum-1).at(j+1);
            index_bead_minus_rotor_minus_BPrime = nextConfigsBPrime.at(beadNum-1).at(j-1);

            generated_samples.clear();
            generated_samples = generateSampleMidBeadMidRotor(index_rotor_plus_B, index_rotor_minus_B,
              index_bead_plus_B,index_bead_plus_rotor_plus_B, index_bead_plus_rotor_minus_B,
              index_bead_minus_B,index_bead_minus_rotor_plus_B, index_bead_minus_rotor_minus_B,
              index_rotor_plus_BPrime,index_rotor_minus_BPrime, index_bead_plus_BPrime,
              index_bead_plus_rotor_plus_BPrime,index_bead_plus_rotor_minus_BPrime,
              index_bead_minus_BPrime,index_bead_minus_rotor_plus_BPrime,
              index_bead_minus_rotor_minus_BPrime);

            nextConfigsB.at(beadNum).at(j) = generated_samples.at(0);
            nextConfigsBPrime.at(beadNum).at(j) = generated_samples.at(1);
        }

        // n = N_B - 1 = 2*N_A - 1
        index_bead_plus_B = nextConfigsB.at(beadNum+1).at(numRotorsSectorB - 1);
        index_bead_minus_B = nextConfigsB.at(beadNum-1).at(numRotorsSectorB - 1);
        index_rotor_minus_B = nextConfigsB.at(beadNum).at(numRotorsSectorB - 2);
        index_bead_plus_rotor_minus_B = nextConfigsB.at(beadNum+1).at(numRotorsSectorB - 2);
        index_bead_minus_rotor_minus_B = nextConfigsB.at(beadNum-1).at(numRotorsSectorB - 2);
        index_bead_plus_BPrime = nextConfigsBPrime.at(beadNum+1).at(numRotorsSectorB - 1);
        index_bead_minus_BPrime = nextConfigsBPrime.at(beadNum-1).at(numRotorsSectorB - 1);
        index_rotor_minus_BPrime = nextConfigsBPrime.at(beadNum).at(numRotorsSectorB - 2);
        index_bead_plus_rotor_minus_BPrime = nextConfigsBPrime.at(beadNum+1).at(numRotorsSectorB - 2);
        index_bead_minus_rotor_minus_BPrime = nextConfigsBPrime.at(beadNum-1).at(numRotorsSectorB - 2);

        generated_samples.clear();
        generated_samples = generateSampleMidBeadLastRotor(index_bead_plus_B, index_rotor_minus_B,
           index_bead_plus_rotor_minus_B, index_bead_minus_B,
           index_bead_minus_rotor_minus_B,index_bead_plus_BPrime,
           index_rotor_minus_BPrime,index_bead_plus_rotor_minus_BPrime,
           index_bead_minus_BPrime,index_bead_minus_rotor_minus_BPrime);

        nextConfigsB.at(beadNum).at(numRotorsSectorB - 1) = generated_samples.at(0);
        nextConfigsBPrime.at(beadNum).at(numRotorsSectorB - 1) = generated_samples.at(1);

    }

    std::vector<std::vector<std::vector<int>>> gibbsSampleSwappedUnswapped(const std::vector<std::vector<std::vector<int>>> &step_configs){

        std::vector<std::vector<int>> nextConfigsA = step_configs.at(0);
        std::vector<std::vector<int>> nextConfigsB = step_configs.at(1);
        std::vector<std::vector<int>> nextConfigsAPrime = step_configs.at(2);
        std::vector<std::vector<int>> nextConfigsBPrime = step_configs.at(3);

        // First bead update
        updateFirstBead(nextConfigsA, nextConfigsB, nextConfigsAPrime, nextConfigsBPrime);

        // Middle beads up to P/2
        for (int i = 1; i < numBeads/2; i++){
            updateMidBeadUnswapped(nextConfigsA, nextConfigsB, nextConfigsAPrime, nextConfigsBPrime, i);
        }

        // TODO: Uncomment this
        sampleExchangeSwappedUnswapped(nextConfigsA, nextConfigsB, nextConfigsAPrime,nextConfigsBPrime);

        if (Swapped){
            // Middle beads (P/2 & P/2 + 1) swapped
            updateMidBeadsSwapped(nextConfigsA, nextConfigsB, nextConfigsAPrime, nextConfigsBPrime);
            CountSwappedMoves += 1;
        }
        else {
            // Middle beads (P/2 & P/2 + 1) unswapped
            updateMidBeadUnswapped(nextConfigsA, nextConfigsB, nextConfigsAPrime, nextConfigsBPrime, numBeads/2);
            updateMidBeadUnswapped(nextConfigsA, nextConfigsB, nextConfigsAPrime, nextConfigsBPrime, numBeads/2 + 1);
            CountUnswappedMoves += 1;
        }

        // Middle beads beyond P/2 + 1
        for (int i = numBeads/2 + 2; i < numBeads; i++){
            updateMidBeadUnswapped(nextConfigsA, nextConfigsB, nextConfigsAPrime, nextConfigsBPrime, i);
        }

        // Last bead
        updateFinalBead(nextConfigsA, nextConfigsB, nextConfigsAPrime, nextConfigsBPrime);

        std::vector<std::vector<std::vector<int>>> nextConfigs;
        nextConfigs.push_back(nextConfigsA);
        nextConfigs.push_back(nextConfigsB);
        nextConfigs.push_back(nextConfigsAPrime);
        nextConfigs.push_back(nextConfigsBPrime);

        CountSwappedMovesVector.push_back(CountSwappedMoves);
        CountUnswappedMovesVector.push_back(CountUnswappedMoves);

        return nextConfigs;
    }

    std::vector<std::vector<std::vector<int>>> gibbsSampleSwapped(const std::vector<std::vector<std::vector<int>>> &step_configs){

        std::vector<std::vector<int>> nextConfigsA = step_configs.at(0);
        std::vector<std::vector<int>> nextConfigsB = step_configs.at(1);
        std::vector<std::vector<int>> nextConfigsAPrime = step_configs.at(2);
        std::vector<std::vector<int>> nextConfigsBPrime = step_configs.at(3);

        updateFirstBead(nextConfigsA, nextConfigsB, nextConfigsAPrime, nextConfigsBPrime);

        // Middle beads up to P/2
        for (int i = 1; i < numBeads/2; i++){
            updateMidBeadUnswapped(nextConfigsA, nextConfigsB, nextConfigsAPrime, nextConfigsBPrime, i);
        }
        // Middle beads (P/2 & P/2 + 1) swapped
        updateMidBeadsSwapped(nextConfigsA, nextConfigsB, nextConfigsAPrime, nextConfigsBPrime);
        // Middle beads beyond P/2 + 1
        for (int i = numBeads/2 + 2; i < numBeads; i++){
            updateMidBeadUnswapped(nextConfigsA, nextConfigsB, nextConfigsAPrime, nextConfigsBPrime, i);
        }

        // Last bead
        updateFinalBead(nextConfigsA, nextConfigsB, nextConfigsAPrime, nextConfigsBPrime);

        std::vector<std::vector<std::vector<int>>> nextConfigs;
        nextConfigs.push_back(nextConfigsA);
        nextConfigs.push_back(nextConfigsB);
        nextConfigs.push_back(nextConfigsAPrime);
        nextConfigs.push_back(nextConfigsBPrime);

        /*
        if (numberSwappedRotors > numberSwappedRotorsInit){
            // Need to suggest unswapping the last rotor of A
            // if unswap is successful: numberSwappedRotors -= 1;
        }
        else {
            // Need to suggest swapping the first rotor of B
            // if swap is successful: numberSwappedRotors += 1;
        }

        return nextConfigs;
        */

        return nextConfigs;
    }

    void integrateSwappedUnswapped(){

        std::vector<std::vector<std::vector<int>>> currentStepConfigs = initialConfigurations;
        for (int i = 1; i < simulationSteps; i++){

            std::vector<std::vector<std::vector<int>>> prevStepConfigs = currentStepConfigs;

            currentStepConfigs = gibbsSampleSwappedUnswapped(prevStepConfigs);

            if (diagnostics){
                phiConfigurations.push_back(currentStepConfigs);
            }
        }

        if (diagnostics){
            OutWriter.outputDiagnosticData(phiConfigurations);
        }

        OutWriter.outputStepData(CountSwappedMovesVector, CountUnswappedMovesVector);

    }
};
