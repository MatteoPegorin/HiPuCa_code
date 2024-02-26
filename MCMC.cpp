#include "MCMC.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <random>
#include <string>
using namespace std;


//Initialize variables that will be also in MCMC_settings.txt
int CHAIN_LENGTH = 200000;
int BURN_IN_LENGTH = CHAIN_LENGTH/5;
vector<double> CHAIN_JUMP_SIZE = {0.05, 0.01, 0.01, 0.00005}; //The CHAIN_JUMP_SIZE should be chosen so that the acceptance rate is around 0.25-0.5



//Random generator --- MAY NEED TO BE EXECUTES SEPARATELY BY EACH THREAD?
default_random_engine generator;


//We implement the Metropolis-Hastings algorithm to perform the MCMC
//Pg 63 Note Liguori and https://www.roma1.infn.it/~dagos/rpp/node42.html


//Takes as input a function P, the log of the Bayesian posterior probability distribution, and samples from it; from the starting point starting_point_parameters in parameter space, and with the given array of data. The return is a vector of vectors, with the points sampled from the MCMC
//You should remove the burn-in of the chain, run several chains, and perform some statistical tests to assess if the chain has converged. The MH_MCMC_jump_size is the standard deviation of the proposal distribution, which is assumed to be a multivariate normal distribution with diagonal covariance matrix
chain MCMC(FunctionPointer log_posterior_probability, vector<double> const & starting_point_parameters, vector<double> const & data, int chain_length, vector<double> const & MH_MCMC_jump_size, int chain_index, bool debug = false){

	cout << "=== Initializing MCMC" << "\n";
	
	chain chain_to_return;

	//starting_point_parameters is the first point in the chain
	chain_to_return.chain_points.push_back(starting_point_parameters);
	chain_to_return.chain_index = chain_index;
	
	vector<double> new_point;
	
	double mean_acceptance_ratio = 0.0;
	const int DEBUG_TO_SKIP = 100;
	for(int i = 0; i < chain_length; i++){
		if(debug && (i + 1) % DEBUG_TO_SKIP == 0){
			cout << "=== MCMC index " << chain_index << " -- iteration " << i << " of " << chain_length << " (" << (100.*i)/chain_length << "%) -- mean acceptance ratio: " << 100.*mean_acceptance_ratio/DEBUG_TO_SKIP << "%" << endl;
			mean_acceptance_ratio = 0.0;
			for(double last_point : chain_to_return.chain_points.back()){
				cout << last_point << "\t";
			}
			cout << endl;
		}

		new_point = draw_point_from_multivariate_normal_distribution(chain_to_return.chain_points.back(), MH_MCMC_jump_size);
		
		/*
			With the Metropolis-Hastings algorithm we have to evaluate the ratio of the two probabilities in these points (and we may save them afterwards for additional performance if we reject the proposal)
			I.e. the acceptance ratio was

			double acceptance_ratio = posterior_probability(new_point, data) / posterior_probability(chain_points.back(), data);
			acceptance_ratio = min(acceptance_ratio, 1.0);

			But, to increase the numerical stability and precision, we can work with the log of the posterior probability
			
			In this case the same formula becomes, from p(y)/p(x), to: exp(log(p(y))/exp(log(p(x)) = exp( log(p(y)) - log(p(x)) )
			Therefore:

			double acceptance_ratio = exp( log_posterior_probability(new_point, data) - log_posterior_probability(chain_points.back(), data) ) ;
			acceptance_ratio = min(acceptance_ratio, 1.0);

		*/

		double acceptance_ratio = exp( log_posterior_probability(new_point, data) - log_posterior_probability(chain_to_return.chain_points.back(), data) ) ;
		acceptance_ratio = min(acceptance_ratio, 1.0);

		//We have to accept the new point with probability acceptance_ratio, therefore we draw a random number from a uniform distribution
		double accept_probability = uniform_real_distribution<double>(0.0, 1.0)(generator);

		if(accept_probability <= acceptance_ratio){
			chain_to_return.chain_points.push_back(new_point);
			mean_acceptance_ratio += 1.0;
		}else{
			chain_to_return.chain_points.push_back(chain_to_return.chain_points.back());
		}

	
	
	}
	
  	return chain_to_return;
  
}

//Function that loads data from file
vector<double> load_data(string filename){

	cout << "=== Loading data from file: " << filename << endl;

	vector<double> data;
	ifstream file(filename);
	double number;
	while (file >> number) {
		data.push_back(number);
	}
	file.close();
	
	cout << "=== Data loaded\n";

	return data;	
}

void initialize_random_number_generator(){
    // Initialize and give seed to random number generator
    cout << "=== IMPLEMENT THIS FUNCTION FOR EACH THREAD SEPARATELY?" << endl;
}

vector<vector<double>> load_initial_paramater_values_MCMC(string def_path){
	//I read the initial points from a text file
	const string filename = "Settings/initial_points.txt";
	vector<vector<double>> initial_values = parseDoubleListFromFile(def_path + filename);
	
	cout << "=== Initial parameter values for MCMC chains:\n";
	for (const auto& line : initial_values) {
        for (const auto& num : line) {
            cout << num << "\t";
        }
        cout << endl;
    }
	
	cout << "=== Initial parameter values for " << initial_values.size() << " chains loaded\n";
	return initial_values;
}

vector<vector<double>> parseDoubleListFromFile(string filename) {
	vector<vector<double>> result;
	ifstream file(filename);
	string line;

	//Throw error if file is not found
	if (!file) {
		cout << "ERROR: File " << filename << " not found!" << endl;
		return {{}};
	}

	while (getline(file, line)) {
        	vector<double> numbers;
        	istringstream iss(line);
        double num;

        while (iss >> num) {
            numbers.push_back(num);
        }

        result.push_back(numbers);
    }

	cout << "File loaded!" << endl;

    return result;
}

//Draw a point from a multivariate normal distribution, assuming diagonal covariance matrix
vector<double> draw_point_from_multivariate_normal_distribution(vector<double> const & initial_point, vector<double> const & std_dev){
	
	vector<double> new_point;
	
	for(unsigned long int i = 0; i < initial_point.size(); i++){
		normal_distribution<double> normal_dist(initial_point[i], std_dev[i]);
		new_point.push_back(normal_dist(generator));
	}
	
	return new_point;
}

//Draw a point from a multivariate normal distribution, assuming covariance matrix proportional to the identity matrix
vector<double> draw_point_from_multivariate_normal_distribution(vector<double> const & initial_point, double std_dev){
	
	//I create a vector with dimension initial_point.size() and all elements equal to std_dev
	vector<double> std_dev_vector(initial_point.size(), std_dev);
	
	return draw_point_from_multivariate_normal_distribution(initial_point, std_dev_vector);
}

//Define log of commonly used probability distribution
double log_uniform(double x, double min, double max){
	if(x < min || x > max){
		return -INFINITY;
	}else{
		return -log(max - min);
	}
}

double log_normal(double x, double mean, double sigma){
	return -0.5 * log(2. * M_PI) - log(sigma) - 0.5 * pow(((x - mean)/sigma), 2.);
}

void load_initial_settings_MCMC(int number_points_per_chain, int chain_burn_in_length, vector<double> chain_jump_sizes){
	CHAIN_LENGTH = number_points_per_chain;
	BURN_IN_LENGTH = chain_burn_in_length;
	CHAIN_JUMP_SIZE = chain_jump_sizes;
}