#include <iostream>
#include <vector>
#include <random>

class test_class{
public:
    std::vector<std::vector<std::vector<int>>> a;

    test_class(){

        for (int j = 0; j < 5; j++){
            std::vector<std::vector<int>> d;
            for (int i = 0; i < 3; i++){
                std::vector<int> b {i+1,i+2,i+3,i+4};
                d.push_back(b);
            }
            a.push_back(d);
        }
    }

    void test_sample(std::vector<std::vector<int>> &step_configs){

        std::vector<std::vector<int>> nextConfigs = step_configs;

        for (int i = 0; i < 3; i++){
            for (int j = 0; j < 4; j++){
                nextConfigs.at(i).at(j) = 0;
            }
        }
        a.push_back(nextConfigs);
    }
};


int main_1() {

    std::string rootFolder = "/Users/shaeermoeed/Github/DVRPairMC/ProbabilityDensities/";
    //ProbabilityDensities probTables = ProbabilityDensities(1.0, 4.0, 10, rootFolder, 3);
    //Grid phiGrid = Grid(10);
    //MonteCarloIntegrator mcSimulator = MonteCarloIntegrator(2, 2, 10, 1.0, 1.0, rootFolder, 3, 1000, false);

    int numBeads = 10;
    int p = 0;
    int test1 = ((p+1) + numBeads) % numBeads;
    int test2 = ((p-1) + numBeads) % numBeads;
    int test3 = (p + numBeads) % numBeads;
    std::cout << test1 << "\n";
    std::cout << test2 << "\n";
    std::cout << test3 << "\n";
    std::cout << "Finished" << "\n";

    test_class test_instance_class = test_class();
    std::vector<std::vector<int>> step_config = test_instance_class.a.at(2);
    test_instance_class.test_sample(step_config);
    std::vector<std::vector<std::vector<int>>> test_a = test_instance_class.a;

    for (int i = 0; i < test_a.size(); i++) {
        std::vector<std::vector<int>> test_b = test_a.at(i);
        for (int j = 0; j < test_b.size(); j++) {
            std::vector<int> test_c = test_b.at(j);
            for (int k = 0; k < test_c.size(); k++){
                std::cout << test_c.at(k) << ",";
            }
            std::cout << "\n";
        }
        std::cout << "\n";
    }

    return 0;
}

int main() {
    std::vector<double> probs;
    for (int i = 0; i < 2; i ++){
        probs.push_back(double(i+1));
    }
    std::discrete_distribution<int> distribution(probs.begin(),probs.end());
    std::mt19937_64 sampleGenerator;
    sampleGenerator.seed(123456789);
    double seedThrowaway = 0.0;
    for (int i = 0; i < 2000; i++){
        seedThrowaway += distribution(sampleGenerator);
    }

    int count_2 = 0;
    int count_1 = 0;
    for (int j = 0; j < 1000; j++){
        int sample = distribution(sampleGenerator);
        if (sample == 1){
            count_2 += 1;
        }
        else if (sample == 0){
            count_1 += 1;
        }
        else {
            std::cout << "Sample = " << sample << "\n";
        }
    }
    std::cout << count_1 << "\n";
    std::cout << count_2 << "\n";
}
