#ifndef MCMC_h
#define MCMC_h


#include "Waveform.h"
#include <vector>
#include <random>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>


using namespace std;


//Dimensionanility of the parameter space
const string path_to_print_results = "Results/";
const int NUM_PAR = 4;

//Not constants anymore, since they can be modified using the MCMC_settings.txt file... it may degrade somewhat the performance, but it is more flexible
//I also have to declare them as extern in the header file, so that they can be used in main.cpp; yet I initialize them do a default value in MCMC.cpp
extern int CHAIN_LENGTH;
extern int BURN_IN_LENGTH;
extern vector<double> CHAIN_JUMP_SIZE; //The CHAIN_JUMP_SIZE should be chosen so that the acceptance rate is around 0.25-0.5


//Defining a type for the pointer to the posterior probability function
typedef double (*FunctionPointer)(const std::vector<double>&, const std::vector<double>&);

struct chain{
	vector<vector<double>> chain_points;
	int chain_index;
	bool removed_burn_in = false;
	
	vector<double> parameter_mean;
	vector<double> parameter_std_dev;

	void remove_burn_in(int burn_in_length){
		if(!removed_burn_in){
			cout << "Removing " << burn_in_length << " chain points due to burn-in from chain " << chain_index << endl;
			chain_points.erase(chain_points.begin(), chain_points.begin() + burn_in_length);
			removed_burn_in = true;
		}else{
			cout << "Burn-in already removed!" << endl;
		}
	}
	
	void set_chain_index(int index){
		chain_index = index;
	}

   	void evaluate_parameter_mean_and_std_dev() {
		if(!removed_burn_in){
			cout << "!!!!!!!! ERROR: Burn-in not removed before evaluating mean and standard deviation" << endl;
		}

        cout << "Evaluating mean and standard deviation of chain " << chain_index << endl;
		int num_points = chain_points.size();
		int num_parameters = chain_points[0].size();
		parameter_mean = vector<double>(num_parameters, 0.0);
		parameter_std_dev = vector<double>(num_parameters, 0.0);
		for (int i = 0; i < num_points; i++) {
			for (int j = 0; j < num_parameters; j++) {
				parameter_mean[j] += chain_points[i][j];
			}
		}

		for (int j = 0; j < num_parameters; j++) {
			parameter_mean[j] /= num_points;
		}

		for (int i = 0; i < num_points; i++) {
			for (int j = 0; j < num_parameters; j++) {
				parameter_std_dev[j] += pow(chain_points[i][j] - parameter_mean[j], 2);
			}
		}

		for (int j = 0; j < num_parameters; j++) {
			parameter_std_dev[j] = sqrt(parameter_std_dev[j] / num_points);
		}
        cout << "Evaluation completed" << endl;
    }

	void print_mean_and_std_dev() {
		cout << "Mean and standard deviation of chain " << chain_index << endl;
		for (int i = 0; i < parameter_mean.size(); i++) {
			cout << "Parameter " << i << " \t| mean:\t" << parameter_mean[i] << " \t| std_dev:\t" << parameter_std_dev[i] << endl;
		}
	}

	void save_chain_to_file(string def_path){
		ofstream file;
		file.open(def_path + path_to_print_results + "Points_chain_" + to_string(chain_index) + ".txt");
		file.precision(7);
		for(int i = 0; i < chain_points.size(); i++){
			for(int j = 0; j < chain_points[i].size(); j++){
				file << chain_points[i][j] << "\t";
			}
			file << endl;
		}
		file.close();
	}

	void print_mean_and_std_dev_to_file(string def_path) {
		ofstream file;
		file.open(def_path + path_to_print_results + "Results_chain_" + to_string(chain_index) + ".txt");
		file << "Chain " << chain_index << " results:" << endl;
		file << "Number of points used from chain: " << chain_points.size() << endl;
		file << "Mean and standard deviation of chain " << chain_index << endl;
		for (int i = 0; i < parameter_mean.size(); i++) {
			file << "Parameter " << i << " \t| mean:\t" << parameter_mean[i] << " \t| std_dev:\t" << parameter_std_dev[i] << endl;
		}
		file.close();
	}

	void process_chain(int burn_in_length, bool debug = false){
		remove_burn_in(burn_in_length);
		evaluate_parameter_mean_and_std_dev();
		if(debug) print_mean_and_std_dev();
	}

};

chain MCMC(FunctionPointer, vector<double> const &, vector<double> const &, int, vector<double> const &, int, bool);
double log_posterior_prob(vector<double> const &, vector<double> const &);
vector<double> load_data(string);
void initialize_random_number_generator();
vector<vector<double>> load_initial_paramater_values_MCMC(string);
vector<vector<double>> parseDoubleListFromFile(string);
vector<double> draw_point_from_multivariate_normal_distribution(vector<double> const &, vector<double> const &);
vector<double> draw_point_from_multivariate_normal_distribution(vector<double> const &, double);
double log_uniform(double, double, double);
double log_gaussian(double, double, double);
void load_initial_settings_MCMC(int, int, vector<double>);

#endif