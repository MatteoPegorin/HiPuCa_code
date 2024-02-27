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

/*
	Parameters for the MCMC
*/
int number_of_chains_to_run = 1;

void execute_MCMC(bool);
void create_fiducial_waveform();
double log_posterior_prob(vector<double> const & parameters, vector<double> const & data);
double log_prior_prob(vector<double> const & parameters);
double log_likelihood_prob(vector<double> const & parameters, vector<double> const & data);
void load_initial_settings();


int main(int argc, char* argv[]){

 	cout << "=== Code for the PhD course: Advanced topics on scientific and parallel programming with practical application on the CAPRI HPC infrastructure" << endl;
	cout << "=== Giovanna Saleh and Matteo Pegorin" << endl;
	
	cout << "=== Initializing program" << endl;
	// Check if path argument is provided, if not, use default path (first argument is the name of the program)
    if (argc >= 3){
		def_path = argv[2];
    }

	cout << "=== Loading MCMC settings" << endl;

	load_initial_settings();

	//If parameter is given at runtime, can override number_of_chains from settings file
	if (argc >= 2){
		number_of_chains_to_run = stoi(argv[1]);
		cout << "=== Overriding number of chains to run with parameter from command line: " << number_of_chains_to_run << endl;
    }

	cout << "Default path: " << def_path << endl;

	cout << "=== Creating fiducial waveform" << endl;

	string command = "mkdir -p " + def_path + "Data";
	int return_value = system(command.c_str());

	create_fiducial_waveform();


	cout << "=== Executing MCMC" << endl;
	command = "mkdir -p " + def_path + "Results";
	return_value = system(command.c_str());

	execute_MCMC(true);

  return 0;
}

void execute_MCMC(bool debug = true){	

	cout << "=== Loading experimental (mock) dataset -- fiducial waveform with gaussian noise" << endl;
	vector<double> data = load_data(def_path+"Data/h_measured_strain_with_noise.txt");
	
	cout << "=== Loading initial parameters values for MCMCs" << endl;
	vector<vector<double>> initial_parameter_values = load_initial_paramater_values_MCMC(def_path);
	if(initial_parameter_values.size() < number_of_chains_to_run){
		cout << "ERROR: Not enough initial parameter values (" << initial_parameter_values.size() << ") for the number of chains to run (" << number_of_chains_to_run << ")!" << endl;
		return;
	}

	vector<chain*> list_MCMC_chains;

	cout << "=== Initializing " << number_of_chains_to_run << " MCMC chains" << endl;
	for(int i = 0; i < number_of_chains_to_run; i++){
		chain* chain_pointer = new chain();
		chain_pointer->initialize_chain(i, initial_parameter_values[i]);
		list_MCMC_chains.push_back(chain_pointer);

		cout << "=== Chain " << i << " initialized" << endl;
	}

	cout << "=== Executing MCMC" << endl;
	
	#pragma omp parallel for
	for(int i = 0; i < list_MCMC_chains.size(); i++){
		// (//Let us process each chain and write the results to file indipendently, as soon as each chain has finished, to improve parallelism)
		// In the following, for each chain we:
		//     - Run the MCMC
		//     - Process the chain (remove burn-in and evaluate statistics)
		//     - Write all the chain points in the output file
		//     - Write the mean and standard deviation of the parameters in the output file
		//     - Deallocate the memory for the chain (as soon as it is not needed anymore, to reduce memory footprint)


		cout << "== Start running chain " << i + 1 << " out of " << list_MCMC_chains.size() << " (with index " << list_MCMC_chains[i]->chain_index << ")" << endl;
		MCMC(list_MCMC_chains[i], log_posterior_prob, data, CHAIN_LENGTH, CHAIN_JUMP_SIZE, debug);
		cout << "== Finished running MCMC chain with index " << list_MCMC_chains[i]->chain_index << endl;

		cout << "== Process MCMC chain index " << list_MCMC_chains[i]->chain_index << " (removing burn-in and evaluating statistics)" << endl;
		list_MCMC_chains[i]->process_chain(BURN_IN_LENGTH, debug);

		cout << "== Writing MCMC chain to file for chain index " << list_MCMC_chains[i]->chain_index << endl;
		list_MCMC_chains[i]->save_chain_to_file(def_path);

		cout << "== Writing MCMC results to file for chain index " << list_MCMC_chains[i]->chain_index << endl;
		list_MCMC_chains[i]->print_mean_and_std_dev_to_file(def_path);

		cout << "== Deallocating memory for chain index " << list_MCMC_chains[i]->chain_index << endl;
		delete list_MCMC_chains[i];

	}

	cout << "=== MCMC execution finished!" << endl;

}

//Not a const anymore as it can be modified in the MCMC_settings.txt file!
double sigma_noise = 1e-21;

void create_fiducial_waveform(){

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

void load_initial_settings(){
	//Check if file exist, and open it
	ifstream file(def_path + "Settings/MCMC_settings.txt");
	string line;

	if(file.is_open()){
		//Number of chains
		getline(file, line);
		getline(file, line);
		number_of_chains_to_run = stoi(line);
		cout << "Number of chains to run: " << number_of_chains_to_run << endl;

		//Number of points for each MCMC chain
		getline(file, line);
		getline(file, line);
		int number_points_per_chain = stoi(line);
		cout << "Number of points per chain: " << number_points_per_chain << endl;

		//Number of burnin points to remove
		getline(file, line);
		getline(file, line);
		int burn_in_length = stoi(line);
		cout << "Burn in length to remove: " << burn_in_length << endl;

		//Chain jump size
		getline(file, line);
		vector<double> jump_sizes = {};
		for(int i = 0; i < 4; i++){
			getline(file, line);
			jump_sizes.push_back(stod(line));
			cout << "Jump size for parameter " << i << ": " << jump_sizes[i] << endl;
		}
		
		//Sigma_noise
		getline(file, line);
		getline(file, line);
		sigma_noise = stod(line);
		cout << "Sigma noise: " << sigma_noise << endl;

		load_initial_settings_MCMC(number_points_per_chain, burn_in_length, jump_sizes);

		cout << "=== MCMC successfully loaded" << endl;

	}else{
		cout << "ERROR: File MCMC_settings.txt not found!" << endl;
	}

	file.close();
}
