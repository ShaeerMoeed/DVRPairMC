#include <vector>
#include "Grid.cpp"
#include <fstream>
#include <iostream>

class OutputWriter{
private:
    int simulationSteps;
    std::string numRotorsA;
    std::string numRotorsB;
    std::string exchangeType;
    std::string exchangeLength;
    int numBeads;
    double couplingStrength;
    double simulationTemperature;

    int lMax;
    Grid phiGrid;
    std::vector<double> phiGridPts;

    std::string outFolderPath;

public:

    void outputDiagnosticData(const std::vector<std::vector<std::vector<std::vector<int>>>> &phi_configs){

        std::string filePath = outFolderPath + "/" + "MC Diagnostic Outputs.csv";
        std::string header = "N_A = " + numRotorsA +
                ", N_B = " + numRotorsB +
                ", Exchange Type = " + exchangeType +
                ", Exchange Length = " + exchangeLength +
                ", P = " + std::to_string(numBeads) +
                ", g = " + std::to_string(couplingStrength) +
                ", T = " + std::to_string(simulationTemperature) +
                ", l = " + std::to_string(lMax) +
                ", Number of Steps = " + std::to_string(simulationSteps) +
                "\n";
        header += "MC Step, Bead Number, Sector, Rotor Number, Rotor Index\n";

        std::ofstream ofs(filePath);
        ofs << header;
        for (int i = 0; i < simulationSteps; i++){
            for (int j = 0; j < numBeads; j++){
                for (int k = 0; k < phi_configs.at(i).at(0).size(); k++){
                    ofs << std::to_string(i+1);
                    ofs << "," << std::to_string(j+1);
                    ofs << "," << "A";
                    ofs << "," << std::to_string(k+1);
                    ofs << "," << std::to_string(phi_configs.at(i).at(0).at(j).at(k));
                    ofs << "\n";
                }
                for (int k = 0; k < phi_configs.at(i).at(1).size(); k++){
                    ofs << std::to_string(i+1);
                    ofs << "," << std::to_string(j+1);
                    ofs << "," << "B";
                    ofs << "," << std::to_string(k+1);
                    ofs << "," << std::to_string(phi_configs.at(i).at(1).at(j).at(k));
                    ofs << "\n";
                }
                for (int k = 0; k < phi_configs.at(i).at(2).size(); k++){
                    ofs << std::to_string(i+1);
                    ofs << "," << std::to_string(j+1);
                    ofs << "," << "A Prime";
                    ofs << "," << std::to_string(k+1);
                    ofs << "," << std::to_string(phi_configs.at(i).at(2).at(j).at(k));
                    ofs << "\n";
                }
                for (int k = 0; k < phi_configs.at(i).at(3).size(); k++){
                    ofs << std::to_string(i+1);
                    ofs << "," << std::to_string(j+1);
                    ofs << "," << "B Prime";
                    ofs << "," << std::to_string(k+1);
                    ofs << "," << std::to_string(phi_configs.at(i).at(3).at(j).at(k));
                    ofs << "\n";
                }
            }
        }
        ofs.close();

    }

    void outputStepData(const std::vector<int> &swapped_counts, const std::vector<int> &unswapped_counts){

        std::string filePath = outFolderPath + "/" + "MC Step Outputs.csv";
        std::string header = "NA = " + numRotorsA + ", NB = " + numRotorsB +
                ", Exchange Type = " + exchangeType + ", Exchange Length = " + exchangeLength +
                ", P = " + std::to_string(numBeads) +
                ", g = " + std::to_string(couplingStrength) + ", T = " +
                std::to_string(simulationTemperature) + ", l = " + std::to_string(lMax) +
                ", Number of Steps = " + std::to_string(simulationSteps) +
                "\n";
        header += "MC Step, N, D \n";
        std::ofstream ofs(filePath);
        ofs << header;
        for (int i = 0; i < simulationSteps; i++){
            ofs << std::to_string(i+1)
                << "," << std::to_string(swapped_counts.at(i))
                << "," << std::to_string(unswapped_counts.at(i)) << "\n";
        }
        ofs.close();

    }

    OutputWriter(const int &sim_steps, const std::string &num_rotors_A, const std::string &num_rotors_B,
                 const int &num_beads, const double &coupling_strength, const int &max_states,
                 const double &simulation_temperature, const std::string &exchange_type,
                 const std::string &exchange_length, const std::string &out_dir):phiGrid(max_states){

        simulationSteps = sim_steps;
        numRotorsA = num_rotors_A;
        numRotorsB = num_rotors_B;
        exchangeType = exchange_type;
        exchangeLength = exchange_length;
        numBeads = num_beads;
        couplingStrength = coupling_strength;
        lMax = max_states;
        simulationTemperature = simulation_temperature;
        phiGridPts = phiGrid.gridPts;
        outFolderPath = out_dir;
        if (numBeads % 2 != 0){
            throw std::runtime_error("Number of beads " + std::to_string(num_beads) + " needs to be even\n");
        }
    }

};

