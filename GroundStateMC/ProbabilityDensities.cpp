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
    int numRotors;
    int numGridPts;

    std::string freeRhoPath;
    std::string pairRhoPath;
    std::string sgnRhoPath;

    void readDensity(const std::string &file_path, std::vector<std::vector<double>> &density_matrix) { // NOLINT(*-convert-member-functions-to-static)

        std::ifstream fileStream(file_path);
        if (!fileStream.good())
            throw std::runtime_error("Could not open file with path: " + file_path);

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

    void readDensity(const std::string &file_path, std::vector<std::vector<int>> &density_matrix) { // NOLINT(*-convert-member-functions-to-static)

        std::ifstream fileStream(file_path);
        if (!fileStream.good())
            throw std::runtime_error("Could not open file with path: " + file_path);

        std::string line, entry;
        std::vector<int> rowData;

        int lineCount = 0;
        while (std::getline(fileStream, line)) {

            lineCount += 1;
            if (lineCount < 3) {
                continue;
            }

            std::stringstream rowStream(line);

            rowData.clear();
            while (std::getline(rowStream, entry, ',')) {
                rowData.push_back(std::stoi(entry));
                if (isnan(rowData.back())) {
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
    std::vector<std::vector<int>> sgnRho;

    ProbabilityDensities(const double &coupling_strength, const double &simulation_temperature, const int &max_states,
                         const int &num_beads, const int &num_rotors, const std::string &density_directory,
                         int density_file_digits) {

        couplingStrength = coupling_strength;
        simulationTemperature = simulation_temperature;
        numBeads = num_beads;
        numRotors = num_rotors;
        effectiveTemperature = (simulationTemperature * numBeads);
        lMax = max_states;
        numGridPts = 2 * lMax + 1;

        std::string effectiveTempString =  std::to_string(effectiveTemperature);
        std::size_t decimal_pos = effectiveTempString.find('.');
        effectiveTempString = effectiveTempString.substr(0, decimal_pos+2);

        std::cout << effectiveTempString << "\n";

        std::string fileParameters = "_Propagator_DVR_l_" + std::to_string(lMax) + "_g_" +
                             std::to_string(couplingStrength).substr(0, density_file_digits)
                             +"_T_" + effectiveTempString
                             + ".csv";
        std::string freeRhoFileName = "Free" + fileParameters;
        std::string pairRhoFileName = "Pair" + fileParameters;
        std::string pairSgnFileName = "Sign" + fileParameters;

        freeRhoPath = density_directory + freeRhoFileName;
        pairRhoPath = density_directory + pairRhoFileName;
        sgnRhoPath = density_directory + pairSgnFileName;

        readDensity(freeRhoPath, freeRho);
        readDensity(pairRhoPath, pairRho);
        readDensity(sgnRhoPath, sgnRho);

        std::cout << "Probability Densities Read" << "\n";

    }

    double pairProbabilityFirstBead(const std::vector<std::vector<int>> &step_configs, const int &rotor_num) {

        double free_term = 1.0;
        double potential_term = 1.0;

        int p = 0;
        int n = rotor_num;

        int index_current = step_configs.at(p).at(n);
        int index_bead_plus = step_configs.at(p + 1).at(n);

        free_term *= freeRho.at(index_current).at(index_bead_plus);

        if (n < numRotors - 1){
            int index_rotor_plus = step_configs.at(p).at(n+1);
            int index_bead_plus_rotor_plus = step_configs.at(p + 1).at(n+1);

            int index_current_bead_rotor_plus = index_current * numGridPts + index_rotor_plus;
            int index_next_bead_rotor_plus = index_bead_plus * numGridPts + index_bead_plus_rotor_plus;

            potential_term *= pairRho.at(index_current_bead_rotor_plus).at(index_next_bead_rotor_plus);
        }
        if (n > 0){
            int index_rotor_minus = step_configs.at(p).at(n-1);
            int index_bead_plus_rotor_minus = step_configs.at(p + 1).at(n-1);

            int index_current_bead_rotor_minus = index_rotor_minus * numGridPts + index_current;
            int index_next_bead_rotor_minus = index_bead_plus_rotor_minus * numGridPts + index_bead_plus;

            potential_term *= pairRho.at(index_current_bead_rotor_minus).at(index_next_bead_rotor_minus);
        }

        return free_term * potential_term;
    }
    double pairProbabilityLastBead(const std::vector<std::vector<int>> &step_configs, const int &rotor_num) {

        double free_term = 1.0;
        double potential_term = 1.0;

        int p = numBeads;
        int n = rotor_num;

        int index_current = step_configs.at(p).at(n);
        int index_bead_minus = step_configs.at(p - 1).at(n);

        free_term *= freeRho.at(index_bead_minus).at(index_current);

        if (n < numRotors - 1){
            int index_rotor_plus = step_configs.at(p).at(n+1);
            int index_bead_minus_rotor_plus = step_configs.at(p - 1).at(n+1);

            int index_current_bead_rotor_plus = index_current * numGridPts + index_rotor_plus;
            int index_prev_bead_rotor_plus = index_bead_minus * numGridPts + index_bead_minus_rotor_plus;

            potential_term *= pairRho.at(index_prev_bead_rotor_plus).at(index_current_bead_rotor_plus);
        }
        if (n > 0){
            int index_rotor_minus = step_configs.at(p).at(n-1);
            int index_bead_minus_rotor_minus = step_configs.at(p - 1).at(n-1);

            int index_current_bead_rotor_minus = index_rotor_minus * numGridPts + index_current;
            int index_prev_bead_rotor_minus = index_bead_minus_rotor_minus * numGridPts + index_bead_minus;

            potential_term *= pairRho.at(index_prev_bead_rotor_minus).at(index_current_bead_rotor_minus);
        }

        return free_term * potential_term;
    }

    double pairProbabilityMiddleBead(const std::vector<std::vector<int>> &step_configs, const int &bead_num,
                                    const int &rotor_num) {

        double free_term = 1.0;
        double potential_term = 1.0;

        int p = bead_num;
        int n = rotor_num;

        int index_current = step_configs.at(p).at(n);
        int index_bead_plus = step_configs.at(p + 1).at(n);
        int index_bead_minus = step_configs.at(p - 1).at(n);

        free_term *= freeRho.at(index_current).at(index_bead_plus);
        free_term *= freeRho.at(index_bead_minus).at(index_current);

        if (n < numRotors - 1){
            int index_rotor_plus = step_configs.at(p).at(n+1);
            int index_bead_plus_rotor_plus = step_configs.at(p + 1).at(n+1);
            int index_bead_minus_rotor_plus = step_configs.at(p - 1).at(n+1);

            int index_current_bead_rotor_plus = index_current * numGridPts + index_rotor_plus;
            int index_next_bead_rotor_plus = index_bead_plus * numGridPts + index_bead_plus_rotor_plus;
            int index_prev_bead_rotor_plus = index_bead_minus * numGridPts + index_bead_minus_rotor_plus;

            potential_term *= pairRho.at(index_current_bead_rotor_plus).at(index_next_bead_rotor_plus);
            potential_term *= pairRho.at(index_prev_bead_rotor_plus).at(index_current_bead_rotor_plus);
        }
        if (n > 0){
            int index_rotor_minus = step_configs.at(p).at(n-1);
            int index_bead_plus_rotor_minus = step_configs.at(p + 1).at(n-1);
            int index_bead_minus_rotor_minus = step_configs.at(p - 1).at(n-1);

            int index_current_bead_rotor_minus = index_rotor_minus * numGridPts + index_current;
            int index_next_bead_rotor_minus = index_bead_plus_rotor_minus * numGridPts + index_bead_plus;
            int index_prev_bead_rotor_minus = index_bead_minus_rotor_minus * numGridPts + index_bead_minus;

            potential_term *= pairRho.at(index_current_bead_rotor_minus).at(index_next_bead_rotor_minus);
            potential_term *= pairRho.at(index_prev_bead_rotor_minus).at(index_current_bead_rotor_minus);
        }
        return free_term * potential_term;
    }

    int extractConfigSign(const std::vector<std::vector<int>> &step_configs){

        int signConfig = 1;
        for (int i = 0; i < numBeads; i++){
            for (int j = 0; j < numRotors-1; j++){
                int stepBeadRotorIndex = step_configs.at(i).at(j);
                int stepBeadRotorPlusIndex = step_configs.at(i).at(j+1);
                int stepBeadPlusRotorIndex = step_configs.at(i+1).at(j);
                int stepBeadPlusRotorPlusIndex = step_configs.at(i+1).at(j+1);

                int indexCurrentBead = stepBeadRotorIndex * numGridPts + stepBeadRotorPlusIndex;
                int indexNextBead = stepBeadPlusRotorIndex * numGridPts + stepBeadPlusRotorPlusIndex;
                signConfig *= sgnRho.at(indexCurrentBead).at(indexNextBead);
            }
        }
        return signConfig;
    }
};
