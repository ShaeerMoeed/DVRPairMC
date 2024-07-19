#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <stdexcept>
#include <sstream>

class ProbabilityDensities {
private:
    double couplingStrength;
    double simulationTemperature;
    double effectiveTemperature;
    int lMax;
    int numBeads;
    int numGridPts;

    std::string freeRhoPath;
    std::string pairRhoPath;

    void readDensity(const std::string &file_path, std::vector<std::vector<double> > &density_matrix) { // NOLINT(*-convert-member-functions-to-static)

        std::ifstream fileStream(file_path);
        if (!fileStream.good())
            throw std::runtime_error("Could not open file with path: " + file_path + "\n");

        std::string line, entry;
        std::vector<double> rowData;

        int lineCount = 0;
        while (std::getline(fileStream, line)) {

            lineCount += 1;
            if (lineCount < 3) {
                continue;
            }

            std::stringstream rowStream(line);

            rowData.clear();
            while (std::getline(rowStream, entry, ',')) {
                rowData.push_back(std::stod(entry));
                if (isnan(rowData.back()) || rowData.back() < 0.0) {
                    throw std::runtime_error("Value is NaN\n");
                }
            }

            density_matrix.push_back(rowData);
        }
        fileStream.close();
    }

public:
    std::vector<std::vector<double>> freeRho;
    std::vector<std::vector<double>> pairRho;

    ProbabilityDensities(const double &coupling_strength, const double &simulation_temperature, const int &max_states,
                         const int &num_beads, const std::string &density_directory, int density_file_digits) {

        couplingStrength = coupling_strength;
        simulationTemperature = simulation_temperature;
        numBeads = num_beads;
        effectiveTemperature = (simulationTemperature * numBeads);
        lMax = max_states;
        numGridPts = 2 * lMax + 1;

        std::string fileParameters = "_Propagator_DVR_l_" + std::to_string(lMax) + "_g_" +
                             std::to_string(couplingStrength).substr(0, density_file_digits)
                             +"_T_" + std::to_string(effectiveTemperature).substr(0, density_file_digits)
                             + ".csv";
        std::string freeRhoFileName = "Free" + fileParameters;
        std::string pairRhoFileName = "Pair" + fileParameters;

        freeRhoPath = density_directory + freeRhoFileName;
        pairRhoPath = density_directory + pairRhoFileName;

        readDensity(freeRhoPath, freeRho);
        readDensity(pairRhoPath, pairRho);

        std::cout << "Probability Densities Read" << "\n";

    }

    double pairProbBead1Rotor1(const int &index_current, const int &index_bead_plus, const int &index_rotor_plus,
                               const int &index_bead_plus_rotor_plus) {

        double free_term = freeRho.at(index_current).at(index_bead_plus);

        int index_current_bead_rotor_plus = index_current * numGridPts + index_rotor_plus;
        int index_next_bead_rotor_plus = index_bead_plus * numGridPts + index_bead_plus_rotor_plus;

        double potential_term = pairRho.at(index_current_bead_rotor_plus).at(index_next_bead_rotor_plus);

        return free_term * potential_term;
    }

    double pairProbBead1RotorMid(const int &index_current, const int &index_bead_plus, const int &index_rotor_plus,
                               const int &index_bead_plus_rotor_plus, const int &index_rotor_minus,
                               const int &index_bead_plus_rotor_minus) {

        double free_term = 1.0;
        double potential_term = 1.0;

        free_term *= freeRho.at(index_current).at(index_bead_plus);

        int index_current_bead_rotor_plus = index_current * numGridPts + index_rotor_plus;
        int index_next_bead_rotor_plus = index_bead_plus * numGridPts + index_bead_plus_rotor_plus;

        potential_term *= pairRho.at(index_current_bead_rotor_plus).at(index_next_bead_rotor_plus);

        int index_current_bead_rotor_minus = index_rotor_minus * numGridPts + index_current;
        int index_next_bead_rotor_minus = index_bead_plus_rotor_minus * numGridPts + index_bead_plus;

        potential_term *= pairRho.at(index_current_bead_rotor_minus).at(index_next_bead_rotor_minus);

        return free_term * potential_term;
    }

    double pairProbBead1RotorN(const int &index_current, const int &index_bead_plus, const int &index_rotor_minus,
                               const int &index_bead_plus_rotor_minus) {

        double free_term = freeRho.at(index_current).at(index_bead_plus);

        int index_current_bead_rotor_minus = index_rotor_minus * numGridPts + index_current;
        int index_next_bead_rotor_minus = index_bead_plus_rotor_minus * numGridPts + index_bead_plus;

        double potential_term = pairRho.at(index_current_bead_rotor_minus).at(index_next_bead_rotor_minus);

        return free_term * potential_term;
    }

    double pairProbBeadMidRotor1(const int &index_current, const int &index_bead_plus, const int &index_bead_minus,
                                 const int &index_rotor_plus, const int &index_bead_plus_rotor_plus,
                                 const int &index_bead_minus_rotor_plus){

        double free_term = 1.0;
        double potential_term = 1.0;

        free_term *= freeRho.at(index_current).at(index_bead_plus);
        free_term *= freeRho.at(index_bead_minus).at(index_current);

        int index_current_bead_rotor_plus = index_current * numGridPts + index_rotor_plus;
        int index_next_bead_rotor_plus = index_bead_plus * numGridPts + index_bead_plus_rotor_plus;
        int index_prev_bead_rotor_plus = index_bead_minus * numGridPts + index_bead_minus_rotor_plus;

        potential_term *= pairRho.at(index_current_bead_rotor_plus).at(index_next_bead_rotor_plus);
        potential_term *= pairRho.at(index_prev_bead_rotor_plus).at(index_current_bead_rotor_plus);

        return free_term * potential_term;

    }

    double pairProbBeadMidRotorN(const int &index_current, const int &index_bead_plus, const int &index_bead_minus,
                                 const int &index_rotor_minus, const int &index_bead_plus_rotor_minus,
                                 const int &index_bead_minus_rotor_minus){

        double free_term = 1.0;
        double potential_term = 1.0;

        free_term *= freeRho.at(index_current).at(index_bead_plus);
        free_term *= freeRho.at(index_bead_minus).at(index_current);

        int index_current_bead_rotor_minus = index_rotor_minus * numGridPts + index_current;
        int index_next_bead_rotor_minus = index_bead_plus_rotor_minus * numGridPts + index_bead_plus;
        int index_prev_bead_rotor_minus = index_bead_minus_rotor_minus * numGridPts + index_bead_minus;

        potential_term *= pairRho.at(index_current_bead_rotor_minus).at(index_next_bead_rotor_minus);
        potential_term *= pairRho.at(index_prev_bead_rotor_minus).at(index_current_bead_rotor_minus);

        return free_term * potential_term;

    }

    double pairProbBeadMidRotorMid(const int &index_current, const int &index_rotor_plus, const int &index_rotor_minus,
                                   const int &index_bead_plus, const int &index_bead_plus_rotor_plus,
                                   const int &index_bead_plus_rotor_minus, const int &index_bead_minus,
                                   const int &index_bead_minus_rotor_plus, const int &index_bead_minus_rotor_minus){

        double free_term = 1.0;
        double potential_term = 1.0;

        free_term *= freeRho.at(index_current).at(index_bead_plus);
        free_term *= freeRho.at(index_bead_minus).at(index_current);

        int index_current_bead_rotor_plus = index_current * numGridPts + index_rotor_plus;
        int index_next_bead_rotor_plus = index_bead_plus * numGridPts + index_bead_plus_rotor_plus;
        int index_prev_bead_rotor_plus = index_bead_minus * numGridPts + index_bead_minus_rotor_plus;

        potential_term *= pairRho.at(index_current_bead_rotor_plus).at(index_next_bead_rotor_plus);
        potential_term *= pairRho.at(index_prev_bead_rotor_plus).at(index_current_bead_rotor_plus);

        int index_current_bead_rotor_minus = index_rotor_minus * numGridPts + index_current;
        int index_next_bead_rotor_minus = index_bead_plus_rotor_minus * numGridPts + index_bead_plus;
        int index_prev_bead_rotor_minus = index_bead_minus_rotor_minus * numGridPts + index_bead_minus;

        potential_term *= pairRho.at(index_current_bead_rotor_minus).at(index_next_bead_rotor_minus);
        potential_term *= pairRho.at(index_prev_bead_rotor_minus).at(index_current_bead_rotor_minus);

        return free_term * potential_term;

    }

    double pairProbBeadMidRotorMidSwappedA(const int &index_current, const int &index_rotor_plus,
                      const int &index_rotor_plus_p, const int &index_rotor_minus, const int &index_bead_plus,
                      const int &index_bead_plus_rotor_plus, const int &index_bead_plus_rotor_minus,
                      const int &index_bead_minus, const int &index_bead_minus_rotor_plus,
                      const int &index_bead_minus_rotor_minus){

        double free_term = 1.0;
        double potential_term = 1.0;

        free_term *= freeRho.at(index_current).at(index_bead_plus);
        free_term *= freeRho.at(index_bead_minus).at(index_current);

        int index_current_bead_rotor_plus = index_current * numGridPts + index_rotor_plus;
        int index_current_bead_rotor_plus_p = index_current * numGridPts + index_rotor_plus_p;
        int index_next_bead_rotor_plus = index_bead_plus * numGridPts + index_bead_plus_rotor_plus;
        int index_prev_bead_rotor_plus = index_bead_minus * numGridPts + index_bead_minus_rotor_plus;

        potential_term *= pairRho.at(index_current_bead_rotor_plus_p).at(index_next_bead_rotor_plus);
        potential_term *= pairRho.at(index_prev_bead_rotor_plus).at(index_current_bead_rotor_plus);

        int index_current_bead_rotor_minus = index_rotor_minus * numGridPts + index_current;
        int index_next_bead_rotor_minus = index_bead_plus_rotor_minus * numGridPts + index_bead_plus;
        int index_prev_bead_rotor_minus = index_bead_minus_rotor_minus * numGridPts + index_bead_minus;

        potential_term *= pairRho.at(index_current_bead_rotor_minus).at(index_next_bead_rotor_minus);
        potential_term *= pairRho.at(index_prev_bead_rotor_minus).at(index_current_bead_rotor_minus);

        return free_term * potential_term;

    }

    double pairProbBeadMidRotorMidSwappedB(const int &index_current, const int &index_rotor_plus,
                                           const int &index_rotor_minus, const int &index_rotor_minus_p,
                                           const int &index_bead_plus, const int &index_bead_plus_rotor_plus,
                                           const int &index_bead_plus_rotor_minus, const int &index_bead_minus,
                                           const int &index_bead_minus_rotor_plus, const int &index_bead_minus_rotor_minus){

        double free_term = 1.0;
        double potential_term = 1.0;

        free_term *= freeRho.at(index_current).at(index_bead_plus);
        free_term *= freeRho.at(index_bead_minus).at(index_current);

        int index_current_bead_rotor_plus = index_current * numGridPts + index_rotor_plus;
        int index_next_bead_rotor_plus = index_bead_plus * numGridPts + index_bead_plus_rotor_plus;
        int index_prev_bead_rotor_plus = index_bead_minus * numGridPts + index_bead_minus_rotor_plus;

        potential_term *= pairRho.at(index_current_bead_rotor_plus).at(index_next_bead_rotor_plus);
        potential_term *= pairRho.at(index_prev_bead_rotor_plus).at(index_current_bead_rotor_plus);

        int index_current_bead_rotor_minus = index_rotor_minus * numGridPts + index_current;
        int index_current_bead_rotor_minus_p = index_rotor_minus_p * numGridPts + index_current;
        int index_next_bead_rotor_minus = index_bead_plus_rotor_minus * numGridPts + index_bead_plus;
        int index_prev_bead_rotor_minus = index_bead_minus_rotor_minus * numGridPts + index_bead_minus;

        potential_term *= pairRho.at(index_current_bead_rotor_minus_p).at(index_next_bead_rotor_minus);
        potential_term *= pairRho.at(index_prev_bead_rotor_minus).at(index_current_bead_rotor_minus);

        return free_term * potential_term;

    }

    double pairProbBeadLastRotor1(const int &index_current, const int &index_bead_minus, const int &index_rotor_plus,
                                  const int &index_bead_minus_rotor_plus){

        double free_term = freeRho.at(index_bead_minus).at(index_current);

        int index_current_bead_rotor_plus = index_current * numGridPts + index_rotor_plus;
        int index_prev_bead_rotor_plus = index_bead_minus * numGridPts + index_bead_minus_rotor_plus;

        double potential_term = pairRho.at(index_prev_bead_rotor_plus).at(index_current_bead_rotor_plus);

        return free_term * potential_term;

    }

    double pairProbBeadLastRotorN(const int &index_current, const int &index_bead_minus, const int &index_rotor_minus,
                                  const int &index_bead_minus_rotor_minus){

        double free_term = freeRho.at(index_bead_minus).at(index_current);

        int index_current_bead_rotor_minus = index_rotor_minus * numGridPts + index_current;
        int index_prev_bead_rotor_minus = index_bead_minus_rotor_minus * numGridPts + index_bead_minus;

        double potential_term = pairRho.at(index_prev_bead_rotor_minus).at(index_current_bead_rotor_minus);

        return free_term * potential_term;

    }

    double pairProbBeadLastRotorMid(const int &index_current, const int &index_bead_minus, const int &index_rotor_plus,
                                    const int &index_bead_minus_rotor_plus, const int &index_rotor_minus,
                                    const int &index_bead_minus_rotor_minus) {

        double free_term = freeRho.at(index_bead_minus).at(index_current);
        double potential_term = 1.0;

        int index_current_bead_rotor_plus = index_current * numGridPts + index_rotor_plus;
        int index_prev_bead_rotor_plus = index_bead_minus * numGridPts + index_bead_minus_rotor_plus;

        potential_term *= pairRho.at(index_prev_bead_rotor_plus).at(index_current_bead_rotor_plus);

        int index_current_bead_rotor_minus = index_rotor_minus * numGridPts + index_current;
        int index_prev_bead_rotor_minus = index_bead_minus_rotor_minus * numGridPts + index_bead_minus;

        potential_term *= pairRho.at(index_prev_bead_rotor_minus).at(index_current_bead_rotor_minus);

        return free_term * potential_term;

    }

    std::vector<double> swapUnswapProbs(const int &index_first_rotor_mid_bead, const int &index_second_rotor_mid_bead,
                           const int &index_third_rotor_mid_bead, const int &index_first_rotor_mid_bead_plus,
                           const int &index_second_rotor_mid_bead_plus, const int &index_third_rotor_mid_bead_plus,
                           const int &index_first_rotor_mid_bead_prime, const int &index_second_rotor_mid_bead_prime,
                           const int &index_third_rotor_mid_bead_prime, const int &index_first_rotor_mid_bead_plus_prime,
                           const int &index_second_rotor_mid_bead_plus_prime, const int &index_third_rotor_mid_bead_plus_prime){

        double free_term_swap = freeRho.at(index_first_rotor_mid_bead_prime).at(index_first_rotor_mid_bead_plus);
        free_term_swap *= freeRho.at(index_first_rotor_mid_bead).at(index_first_rotor_mid_bead_plus_prime);
        free_term_swap *= freeRho.at(index_second_rotor_mid_bead_prime).at(index_second_rotor_mid_bead_plus);
        free_term_swap *= freeRho.at(index_second_rotor_mid_bead).at(index_second_rotor_mid_bead_plus_prime);

        int pair_index_rotor_01_mid_bead_all_prime = index_first_rotor_mid_bead_prime * numGridPts + index_second_rotor_mid_bead_prime;
        int pair_index_rotor_01_mid_bead_plus = index_first_rotor_mid_bead_plus * numGridPts + index_second_rotor_mid_bead_plus;
        double pair_term_swap = pairRho.at(pair_index_rotor_01_mid_bead_all_prime).at(pair_index_rotor_01_mid_bead_plus);

        int pair_index_rotor_01_mid_bead = index_first_rotor_mid_bead * numGridPts + index_second_rotor_mid_bead;
        int pair_index_rotor_01_mid_bead_plus_all_prime = index_first_rotor_mid_bead_plus_prime * numGridPts +
                index_second_rotor_mid_bead_plus_prime;
        pair_term_swap *= pairRho.at(pair_index_rotor_01_mid_bead).at(pair_index_rotor_01_mid_bead_plus_all_prime);

        int pair_index_rotor_12_mid_bead_first_prime = index_second_rotor_mid_bead_prime * numGridPts + index_third_rotor_mid_bead;
        int pair_index_rotor_12_mid_bead_plus = index_second_rotor_mid_bead_plus * numGridPts + index_third_rotor_mid_bead_plus;
        pair_term_swap *= pairRho.at(pair_index_rotor_12_mid_bead_first_prime).at(pair_index_rotor_12_mid_bead_plus);

        int pair_index_rotor_12_mid_bead_second_prime = index_second_rotor_mid_bead * numGridPts + index_third_rotor_mid_bead_prime;
        int pair_index_rotor_12_mid_bead_plus_all_prime = index_second_rotor_mid_bead_plus_prime * numGridPts +
                index_third_rotor_mid_bead_plus_prime;
        pair_term_swap *= pairRho.at(pair_index_rotor_12_mid_bead_second_prime).at(pair_index_rotor_12_mid_bead_plus_all_prime);

        double free_term_unswap = freeRho.at(index_first_rotor_mid_bead).at(index_first_rotor_mid_bead_plus);
        free_term_unswap *= freeRho.at(index_first_rotor_mid_bead_prime).at(index_first_rotor_mid_bead_plus_prime);
        free_term_unswap *= freeRho.at(index_second_rotor_mid_bead).at(index_second_rotor_mid_bead_plus);
        free_term_unswap *= freeRho.at(index_second_rotor_mid_bead_prime).at(index_second_rotor_mid_bead_plus_prime);

        double pair_term_unswap = pairRho.at(pair_index_rotor_01_mid_bead).at(pair_index_rotor_01_mid_bead_plus);
        pair_term_unswap *= pairRho.at(pair_index_rotor_01_mid_bead_all_prime).at(pair_index_rotor_01_mid_bead_plus_all_prime);

        int pair_index_rotor_12_mid_bead = index_second_rotor_mid_bead * numGridPts + index_third_rotor_mid_bead;
        pair_term_unswap *= pairRho.at(pair_index_rotor_12_mid_bead).at(pair_index_rotor_12_mid_bead_plus);

        int pair_index_rotor_12_mid_bead_all_prime = index_second_rotor_mid_bead_prime * numGridPts + index_third_rotor_mid_bead_prime;
        pair_term_unswap *= pairRho.at(pair_index_rotor_12_mid_bead_all_prime).at(pair_index_rotor_12_mid_bead_plus_all_prime);

        double swap_prob = free_term_swap * pair_term_swap;
        double unswap_prob = free_term_unswap * pair_term_unswap;

        double normalization = swap_prob + unswap_prob;
        swap_prob /= normalization;
        unswap_prob /= normalization;

        std::vector<double> swap_unwap_probs;

        swap_unwap_probs.push_back(swap_prob);
        swap_unwap_probs.push_back(unswap_prob);

        return swap_unwap_probs;
    }

};
