cmake_minimum_required(VERSION 3.27)
project(HighTemperatureMC)

set(CMAKE_CXX_STANDARD 20)

add_executable(HighTemperatureMC main.cpp
        ProbabilityDensities.cpp
        Estimators.cpp
        Grid.cpp
        MonteCarloIntegrator.cpp
)
