#include<vector>

class Grid{
private:
    int lMax;
    int numGridPts;
    double gridSpacing;

public:
    std::vector<double> gridPts;

    Grid(const int &max_states){ // NOLINT(*-explicit-constructor)

        lMax = max_states;
        numGridPts = 2 * lMax + 1;
        gridSpacing = 2*M_PI/numGridPts;

        for (int i = 0; i < numGridPts; i++){
            gridPts.push_back(i * gridSpacing);
        }

    }

};
