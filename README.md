# EWE_project

The repository is structured as such:
1) The notebook EWE_project_SEIRS_simulation tests the SEIR model built.
2) The Parallelized_SEIR_model.jl runs the stochastic simulation in a Julia environment. An analogous version of the code was devised in a C++ environment.
3) params.jl defines the SEIR model parameters for the Julia code.
4) Params_production.py instead produces an analogous .txt file with all the parameters to run the SEIR model simulation (specifically for the C++ simulation). It takes in input the parameters defined in params_values.py
