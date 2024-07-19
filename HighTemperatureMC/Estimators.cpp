#include <vector>
#include <cmath>
#include "Grid.cpp"
#include <fstream>
#include <iostream>

class Estimators{
private:
    int simulationSteps;
    int numRotors;
    int numBeads;
    double couplingStrength;
    double simulationTemperature;
    int numBlocks;
    int blockNum;

    std::vector<double> stepwiseCorrelation;
    std::vector<double> stepwisePotential;
    std::vector<double> stepwiseBinder;

    int lMax;
    Grid phiGrid;
    std::vector<double> phiGridPts;

    std::string outFolderPath;

public:
    void updateOrientationalCorrelation(const std::vector<std::vector<int>> &stepConfigurations){

        std::vector<int> firstBeadConfigs = stepConfigurations.at(0);

        double stepCorrelation = 0.0;
        for (int i = 0; i < numRotors-1; i++){
            double phi1 = phiGridPts.at(firstBeadConfigs.at(i));
            double phi2 = phiGridPts.at(firstBeadConfigs.at(i+1));
            stepCorrelation += cos(phi1 - phi2);
        }

        stepwiseCorrelation.push_back(stepCorrelation);
    }

    void updatePotentialEnergy(const std::vector<std::vector<int>> &stepConfigurations){

        std::vector<int> firstBeadConfigs = stepConfigurations.at(0);

        double stepPotential = 0.0;
        for (int i = 0; i < numRotors-1; i++){
            double phi_i = phiGridPts.at(firstBeadConfigs.at(i));
            double phi_j = phiGridPts.at(firstBeadConfigs.at(i+1));
            stepPotential += (sin(phi_i)*sin(phi_j) - 2.0 * cos(phi_i)*cos(phi_j));
        }

        stepwisePotential.push_back(couplingStrength * stepPotential);
    }

    void updateAllProps(const std::vector<std::vector<int>> &stepConfigurations) {

        std::vector<int> firstBeadConfigs = stepConfigurations.at(0);

        double stepPotential = 0.0;
        double stepCorrelation = 0.0;
        double stepKinetic = 0.0;

        double phi_0 = phiGridPts.at(firstBeadConfigs.at(0));

        double stepBinder = cos(phi_0);

        for (int i = 0; i < numRotors-1; i++){
            double phi_i = phiGridPts.at(firstBeadConfigs.at(i));
            double phi_j = phiGridPts.at(firstBeadConfigs.at(i+1));
            stepPotential += (sin(phi_i) * sin(phi_j) - 2.0 * cos(phi_i) * cos(phi_j));
            stepCorrelation += cos(phi_i - phi_j);
            stepBinder += cos(phi_j);
        }
        double binderRatio = std::pow(stepBinder, 2);

        stepwisePotential.push_back(couplingStrength * stepPotential);
        stepwiseCorrelation.push_back(stepCorrelation);
        stepwiseBinder.push_back(binderRatio);
    }

    void outputStepData(){

        std::string filePath = outFolderPath + "/" + "MC Step Outputs.csv";
        std::string header = "N = " + std::to_string(numRotors) + ", P = " + std::to_string(numBeads) +
                ", g = " + std::to_string(couplingStrength) + ", T = " +
                std::to_string(simulationTemperature) + ", l = " + std::to_string(lMax) +
                ", Block = " + std::to_string(blockNum) + ", Number of Blocks = " +
                std::to_string(numBlocks) + ", Number of Steps = " + std::to_string(simulationSteps) +
                "\n";
        header += "MC Step, Potential Energy, Correlation, Binder Ratio\n";
        std::ofstream ofs(filePath);
        ofs << header;
        double potential_sum = 0.0;
        for (int i = 0; i < simulationSteps; i++){
            ofs << std::to_string(i+1) << "," << std::to_string(stepwisePotential.at(i))
            << "," << std::to_string(stepwiseCorrelation.at(i))
            << "," << std::to_string(stepwiseBinder.at(i)) << "\n";
            potential_sum += stepwisePotential.at(i);
        }
        double potential_avg = potential_sum/simulationSteps;
        ofs.close();

        std::cout << "Potential Average = " << potential_avg << "\n";
    }

    Estimators(const int &sim_steps, const int &num_rotors, const int &num_beads, const double &coupling_strength,
               const int &max_states, const double &simulation_temperature, const int &block_num, const int &num_blocks,
               const std::string &out_dir):phiGrid(max_states){

        simulationSteps = sim_steps;
        numRotors = num_rotors;
        numBeads = num_beads;
        couplingStrength = coupling_strength;
        lMax = max_states;
        simulationTemperature = simulation_temperature;
        blockNum = block_num;
        numBlocks = num_blocks;
        phiGridPts = phiGrid.gridPts;
        outFolderPath = out_dir;
    }

};

