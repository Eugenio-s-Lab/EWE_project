#include <iostream>
#include <random>
#include <vector>
#include <fstream>
#include <string>
#include <sstream>
#include <iomanip>
//#include <omp.h> // Include the OpenMP header


typedef std::vector<std::vector<double>> colocation_type;  // NxN matrix
typedef std::vector<int> count_vector;  // N vector of integers
typedef std::vector<count_vector>  run_container; //  T vectors of N 



double compute_p_beta(int i, double beta, double dt, int N, const colocation_type& coloc, const count_vector& I) {
    double p = 0.;
    for (int j=0; j<N; j++) {
        p += coloc[i][j] * double(I[j]);
    }
    return 1. - exp(-beta*dt*p);
}




int main(int argc, char* argv[]) {
    if (argc != 7) {
        std::cerr << "Wrong command-line arguments provided." << std::endl;
        return 1;
    }

    // read the integer parameters
    std::vector<int> parameters_int(6);
    std::ifstream read_param_int(argv[5]);
    if (!read_param_int.is_open()) {
            std::cerr << "Error opening file: " << argv[5] << std::endl;
        }
    for (int i = 0; i < 6; ++i) {
        if (!(read_param_int >> parameters_int[i])) {
                    std::cerr << "Error reading from file: " << argv[5] << std::endl;
                }
    }  // end i
    read_param_int.close();

    // read the float parameters
    std::vector<double> parameters_double(2);
    std::ifstream read_param_double(argv[4]);
    if (!read_param_double.is_open()) {
            std::cerr << "Error opening file: " << argv[4] << std::endl;
        }
    for (int i = 0; i < 2; ++i) {
        if (!(read_param_double >> parameters_double[i])) {
                    std::cerr << "Error reading from file: " << argv[4] << std::endl;
                }
    }  // end i
    read_param_double.close();

    // read the address list
    std::vector<std::string> address_file(2);
    std::ifstream read_param_address(argv[3]);
    if (!read_param_address.is_open()) {
            std::cerr << "Error opening file: " << argv[3] << std::endl;
        }
    for (int i = 0; i < 2; ++i) {
        if (!(read_param_address >> address_file[i])) {
                    std::cerr << "Error reading from file: " << argv[3] << std::endl;
                }
    }  // end i
    read_param_address.close();
    
    // PARAMETERS
    double R_local = std::stod(argv[1]);
    double alpha = parameters_double[0];
    double mu = parameters_double[1];
    int N = parameters_int[0];  // number of patches
    int n_runs = std::stod(argv[2]);  // number of stochastic runs
    int T = parameters_int[1];  // max time (days)
    int T_coloc = parameters_int[2];  // max colocation
    int t_start = parameters_int[3];  // day to start the simulation
    int index_seed = parameters_int[4];  // index of the patch in which the infection starts
    int initial_infected = parameters_int[5];  // number of initial infected
    std::string file_coloc_template = address_file[0];  // directory of the colocation files up to the index. eg. "./baseline/" 
    std::string file_pop = address_file[1]; // filename and path of the pop file
    std::string output_address = argv[6]; //

    // FIXED PARAMETERS
    int colocation_periodicity = 7;   // new colocation every xx days
    double dt = 1.0;  // discretization time (=1 if time is in days)


    ////////
    // number of available threads
    //int n_cpu = omp_get_num_threads();
    //std::cout << "Number of threads: " << n_cpu << std::endl;

    ////////
    // read the colocation
    std::vector<colocation_type> colocation(T_coloc, colocation_type(N, std::vector<double>(N, 0.0)));   // coloc[t][i][j] where t=colocation time index and i,j patches 
    colocation_type colocation_baseline(N, std::vector<double>(N, 0.0));

    std::string file_coloc = file_coloc_template + "coloc_matrix_baseline.txt";
        std::ifstream file(file_coloc);
        if (!file.is_open()) {
            std::cerr << "Error opening file: " << file_coloc << std::endl;
            return 1;
	}

        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                if (!(file >> colocation_baseline[i][j])) {
                    std::cerr << "Error reading from file: " << file_coloc << std::endl;
                    return 1;
		}
            }
        }
    for (int c=0; c<T_coloc; c++) {
        std::string file_coloc = file_coloc_template + "coloc_matrix_" + std::to_string(c) + ".txt";  // eg. "./baseline/coloc_3.txt"
        std::ifstream file(file_coloc);
        if (!file.is_open()) {
            std::cerr << "Error opening file: " << file_coloc << std::endl;
            return 1;
	}

        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                if (!(file >> colocation[c][i][j])) {
                    std::cerr << "Error reading from file: " << file_coloc << " at position (" << i << ", " << j << ")" << std::endl;
                    return 1;
		}
            }  // end for j
        } // end for i

        file.close();

    }  // end for c

    ////////
    // read the population
    std::vector<int> population(N);
    std::ifstream rin(file_pop);
    if (!rin.is_open()) {
            std::cerr << "Error opening file: " << file_pop << std::endl;
        }
    for (int i = 0; i < N; ++i) {
        if (!(rin >> population[i])) {
                    std::cerr << "Error reading from file: " << file_pop << std::endl;
                }
    }  // end i
    rin.close();

    double beta = (R_local*mu)/(colocation[0][index_seed][index_seed]*population[index_seed]);

    ////////////
    // base random seed
    unsigned int seed0 = static_cast<unsigned int>(time(0));  


    ////////
    // compute the probas that I can compute
    double p_mu = 1.-exp(-mu*dt);
    double p_alpha = 1.-exp(-alpha*dt);

    std::vector<run_container> attack_rate(n_runs, run_container(T, count_vector(N, 0)) );


    #pragma omp parallel for schedule(dynamic) // Parallelize the loop
    for (int it=0; it<n_runs; it++) {

        std::mt19937 gen(seed0+it);  // initialize the random seed, making sure no run has the same seed

        std::vector<count_vector> lE(T, count_vector(N, 0));  // lE[t][i]
        std::vector<count_vector> lI(T, count_vector(N, 0));
        std::vector<count_vector> lR(T, count_vector(N, 0));
        int E_new, I_new, R_new, S;
        double p_beta;
        int n_active_people;

        // initialize outbreak
        lI[t_start][index_seed] = initial_infected;
        attack_rate[it][t_start][index_seed] = initial_infected;

        // loop over times
        for (int t=t_start+1; t<T; t++) {

            // loop over patches
            n_active_people = 0;
            for (int i=0; i<N; i++) {

                // E to I
                if (lE[t-1][i]>0) {
                    std::binomial_distribution<> dEI(lE[t-1][i], p_alpha);
                    I_new = dEI(gen); 
                } else {
                    I_new = 0;
                }

                // I to R
                if (lI[t-1][i]>0) {
                    std::binomial_distribution<> dIR(lI[t-1][i], p_mu);
                    R_new = dIR(gen); 
                } else {
                    R_new = 0;
                }

                // S to E
                S = population[i]-lE[t-1][i]-lI[t-1][i]-lR[t-1][i];
                if (S>0) {
                    p_beta = compute_p_beta(i, beta, dt, N, colocation[t/colocation_periodicity], lI[t-1]);
                    std::binomial_distribution<> dSE(S, p_beta);
                    E_new = dSE(gen);
                } else {
                    E_new = 0;
                }
                
                // update
                lE[t][i] = lE[t-1][i] - I_new + E_new;
                lI[t][i] = lI[t-1][i] - R_new + I_new;
                lR[t][i] = lR[t-1][i] + R_new;
                attack_rate[it][t][i] = lE[t][i]+lI[t][i]+lR[t][i];
                n_active_people += lE[t][i] + lI[t][i];
    
            }  // end for i

            // 
            if (n_active_people == 0) {
                // the outbreak has died out. Propagate the number of recovered and stop
                for (int t2=t+1; t2<T; t2++) {
                    for (int i=0; i<N; i++) {
                        attack_rate[it][t2][i] = attack_rate[it][t][i];
                    }
                }
                break;
            }


        }  // end for t
        
        int integer_part = static_cast<int>(std::floor(R_local));
        double decimal_part = std::round(100*(R_local - std::floor(R_local)));

        std::ostringstream file_path_stream;
        file_path_stream << "gzip - > " << output_address << "/I+E+R_output_" << it+1 << "_R_local_" << integer_part << "_"
                        << decimal_part << ".txt.gz";
        std::string file_path = file_path_stream.str();

        FILE *pipe = popen(file_path.c_str(), "w");

        for (int i=0; i<T; i++) { 
            for (int j=0; j<N-1; j++) {
                fprintf(pipe, "%d ", attack_rate[it][i][j]);
            }
            fprintf(pipe, "%d\n", attack_rate[it][i][N-1]);
        }

        pclose(pipe);

    }  // end for it


    return 0;
}
