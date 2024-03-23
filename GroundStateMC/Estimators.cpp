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
    int midBeadIndex;

    std::vector<double> stepwiseCorrelation;
    std::vector<double> stepwisePotential;
    std::vector<double> stepwiseBinder;
    std::vector<double> stepwiseEnergy;

    int lMax;
    Grid phiGrid;
    std::vector<double> phiGridPts;

    std::string outFolderPath;

public:

    void updateAllProps(const std::vector<std::vector<int>> &stepConfigurations) {

        std::vector<int> middleBeadConfigs = stepConfigurations.at(midBeadIndex);
        std::vector<int> finalBeadConfigs = stepConfigurations.at(numBeads);

        double stepPotential = 0.0;
        double stepCorrelation = 0.0;
        double stepEnergy = 0.0;
        //double stepBinderNum = std::pow(cos(phiGridPts.at(middleBeadConfigs.at(0))),4);
        //double stepBinderDenom = std::pow(cos(phiGridPts.at(middleBeadConfigs.at(0))),2);
        double stepBinder = cos(phiGridPts.at(middleBeadConfigs.at(0)));;
        for (int i = 0; i < numRotors-1; i++){
            double phi_i = phiGridPts.at(middleBeadConfigs.at(i));
            double phi_j = phiGridPts.at(middleBeadConfigs.at(i+1));
            stepPotential += (sin(phi_i) * sin(phi_j) - 2.0 * cos(phi_i) * cos(phi_j));
            stepCorrelation += cos(phi_i - phi_j);
            //stepBinderNum += std::pow(cos(phi_j),4);
            //stepBinderDenom += std::pow(cos(phi_j),2);
            stepBinder += cos(phi_j);

            double phi_1 = phiGridPts.at(finalBeadConfigs.at(i));
            double phi_2 = phiGridPts.at(finalBeadConfigs.at(i+1));
            stepEnergy += (sin(phi_1) * sin(phi_2) - 2.0 * cos(phi_1) * cos(phi_2));
        }
        double binderRatio = std::pow(stepBinder, 2);

        stepwisePotential.push_back(couplingStrength * stepPotential);
        stepwiseCorrelation.push_back(stepCorrelation);
        stepwiseBinder.push_back(binderRatio);
        stepwiseEnergy.push_back(stepEnergy);
    }

    void outputStepData(){

        std::string filePath = outFolderPath + "/" + "MC Step Outputs.csv";
        std::string header = "N = " + std::to_string(numRotors) + ", P = " + std::to_string(numBeads) +
                ", g = " + std::to_string(couplingStrength) + ", T = " +
                std::to_string(simulationTemperature) + ", l = " + std::to_string(lMax) +
                ", Block = " + std::to_string(blockNum) + ", Number of Blocks = " +
                std::to_string(numBlocks) + ", Number of Steps = " + std::to_string(simulationSteps) +
                "\n";
        header += "MC Step, Potential Energy, Correlation, Binder Ratio, Energy\n";
        std::ofstream ofs(filePath);
        ofs << header;
        for (int i = 0; i < simulationSteps; i++){
            ofs << std::to_string(i+1)
            << "," << std::to_string(stepwisePotential.at(i))
            << "," << std::to_string(stepwiseCorrelation.at(i))
            << "," << std::to_string(stepwiseBinder.at(i))
            << "," << std::to_string(stepwiseEnergy.at(i)) << "\n";
        }
        ofs.close();

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
        midBeadIndex = numBeads/2;
        if (numBeads % 2 != 0){
            throw std::runtime_error("Number of beads " + std::to_string(num_beads) + " needs to be even\n");
        }
    }

};

