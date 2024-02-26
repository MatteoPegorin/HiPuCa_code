#include "MCMC.h"
#include "Waveform.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <random>
#include <string>

using namespace std;

//Default path where to write temporary data and results
string def_path = "./";



//To_FIX
int execute_MCMC();
int main_waveform();
int create_fiducial_waveform();
double log_posterior_prob(vector<double> const & parameters, vector<double> const & data);
double log_prior_prob(vector<double> const & parameters);
double log_likelihood_prob(vector<double> const & parameters, vector<double> const & data);

int main(int argc, char* argv[]){

 	cout << "=== Code for the PhD course: Advanced topics on scientific and parallel programming with practical application on the CAPRI HPC infrastructure" << endl;
	cout << "=== Giovanna Saleh and Matteo Pegorin" << endl;
	
	cout << "=== Initializing program" << endl;
	// Check if path argument is provided, if not, use default path (first argument is the name of the program)
    if (argc == 2){
		def_path = argv[1];
    }

	cout << "Default path: " << def_path << endl;

	cout << "=== Creating fiducial waveform" << endl;

	string command = "mkdir -p " + def_path + "Data";
	int return_value = system(command.c_str());

	create_fiducial_waveform();

	cout << "=== Executing MCMC" << endl;
	command = "mkdir -p " + def_path + "Results";
	return_value = system(command.c_str());

	execute_MCMC();

  return 0;
}

int execute_MCMC(){	

	cout << "=== Initializing random number generator" << endl;
	initialize_random_number_generator();

	cout << "=== Loading experimental (mock) dataset -- fiducial waveform with gaussian noise" << endl;
	vector<double> data = load_data(def_path+"Data/h_measured_strain_with_noise.txt");
	
	cout << "=== Loading initial parameters values for MCMCs" << endl;
	vector<vector<double>> initial_parameter_values = load_initial_paramater_values_MCMC();

	cout << "=== Calling MCMC function" << endl;
	chain executed_chain = MCMC(log_posterior_prob, initial_parameter_values[0], data, CHAIN_LENGTH, CHAIN_JUMP_SIZE, 0, true);

	cout << "=== DEBUG: print MCMC result" << endl;
	executed_chain.process_chain(BURN_IN_LENGTH, true);

	cout << "=== Writing MCMC results to file" << endl;
	executed_chain.print_mean_and_std_dev_to_file(def_path);

	cout << "=== Writing MCMC chain to file" << endl;
	executed_chain.save_chain_to_file(def_path);

	cout << "=== MCMC executed" << endl;

	return 0;
}


const double sigma_noise = 1e-21;

int create_fiducial_waveform(){

	cout << "=== Creating fiducial post-Newtonian waveform for the inspiral, using the TaylorT3 approximant" << endl;


	// Fiducial Parameters
	double m1 = 80.; //solar masses
	double m2 = 60.; //solar masses
	double t_c = t_max*0.7; //s
	double phi_0 = 0.; //rad, phase at coalescence, could be degenerate with other parameters, from 0 to 2 pi
	double r_Mpc = 600.; //Mpc
	double iota = 0.5; //rad, from 0 to pi
	double obs_theta = 0.4; //rad, from 0 to pi
	double obs_phi = 1.2; //rad, from 0 to 2 pi

	cout << "=== Evaluating waveform" << endl;
	vector<double> h_strain = h_measured_strain(obs_theta, obs_phi, iota, phi_0, t_c, r_Mpc, m1, m2);

	cout << "=== Waveform evaluated" << endl;

	cout << "=== Writing waveform to file" << endl;
	write_vector_to_file(h_strain, def_path+"Data/h_measured_strain.txt");

	cout << "=== Adding noise to the waveform" << endl;

	cout << "=== sigma_noise: " << sigma_noise << endl;


	default_random_engine generator_temp;
	vector<double> h_strain_with_noise;

	for(unsigned long int i = 0; i < h_strain.size(); i++){
		normal_distribution<double> normal_dist(h_strain[i], sigma_noise);
		h_strain_with_noise.push_back(normal_dist(generator_temp));
	}
	
	cout << "=== Writing waveform with noise to file" << endl;
	write_vector_to_file(h_strain_with_noise, def_path+"Data/h_measured_strain_with_noise.txt");

  return 0;
}

// log_posterior_prob(vector<double> parameters, vector<double> data) return the unnormalized value of the log of posterior, given a specic array of values for the parameters (in parameters), and the data for the likelihood in the vector data
/*
	parameters to be estimated are:
	0: r_Mpc
	1: m1
	2: m2
	3: t_c

*/

double log_prior_prob(vector<double> const & parameters){
	double log_prior = 0.;

	log_prior += log_uniform(parameters[0], 100, 2000); //r_Mpc ... you should use loglog_uniform
	log_prior += log_uniform(parameters[1], 10., 100.); //m1
	log_prior += log_uniform(parameters[2], 10., 100.); //m2
	log_prior += log_uniform(parameters[3], 0., t_max); //t_c

	return log_prior;
}

double log_likelihood_prob(vector<double> const & parameters, vector<double> const & data){
	double log_likelihood = 0.;
	
	vector<double> theoretical_waveform = h_measured_strain(0.4, 1.2, 0.5, 0., parameters[3], parameters[0], parameters[1], parameters[2]);

	for(unsigned long int i = 0; i < data.size(); i++){
		log_likelihood += pow( (data[i] - theoretical_waveform[i])/(sigma_noise) , 2.0);
	}

	log_likelihood = -0.5 * log_likelihood - data.size() * log(sigma_noise);

	return log_likelihood;
}

double log_posterior_prob(vector<double> const & parameters, vector<double> const & data){

	double log_prior = log_prior_prob(parameters);
	double log_likelihood = log_likelihood_prob(parameters, data);

	//Multiply log likelihood by log prior (so sum them) to get log posterior
	double log_posterior = log_likelihood + log_prior;

	return log_posterior;
}
