#include "Waveform.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <random>
#include <string>
using namespace std;

/*
    Using master thesis formulae:
  	1.91, 1.92, 1.109, 1.110, 6.39
*/

vector<double> h_measured_strain(double obs_theta, double obs_phi, double iota, double phi_0, double t_c, double r_Mpc, double m1, double m2)
{
  //phi_0 may be degenerate with otgher parameters

  vector<double> h_measured_strain;


  //Observation direction defined in page 35 of master thesis
  double Fplus = 0.5 * (1 + pow(cos(obs_theta), 2.)) * cos(2 * obs_phi);
  double Fcross = cos(obs_theta) * sin(2 * obs_phi);
  double nu = nu_source(m1, m2);
  double m = m_source(m1, m2);
  double mu = mu_source(m1, m2);

  for(double t = 0.; t < t_max; t+=deltat){
    double tau = t_c - t;
    double R_over_Rs = 0.;
    if(tau > 0.){
      R_over_Rs = pow( (16./5. * speed_of_light_cubed_over_G_over_solar_mass * mu / (m*m) * tau) , 1./4.);
    } 
    const double R_ISCO_over_Rs = 3.;

    double h_plus = 0.;
    double h_cross = 0.;

    //Need to implement cutoff at Risco!!! Then I simply return vanishing amplitude and do not compute the phase
    if(R_over_Rs > R_ISCO_over_Rs){
      double Theta = Theta_val(t, t_c, nu, m);
      double phase = phase_2PN(Theta, phi_0, nu);

      
      double amplitude = A_amplitude(r_Mpc, m1, m2, t_c - t);

      //I try to smoothly transition to zero amplitude
      const double multiple_of_R_ISCO = 1.1;
      if(R_over_Rs < multiple_of_R_ISCO* R_ISCO_over_Rs){
        double xvar = ((R_over_Rs - R_ISCO_over_Rs) /((multiple_of_R_ISCO-1.)*R_ISCO_over_Rs));
        amplitude *= (cos((1. - xvar)*M_PI) * 0.5 + 0.5);
      }

      h_plus = amplitude * 0.5 * (1 + cos(iota) * cos(iota)) * cos(phase); 
      h_cross = amplitude * (cos(iota)) * sin(phase);  
    }
    
    double h_measured_strain_temp = Fplus * h_plus + Fcross * h_cross;
    h_measured_strain.push_back(h_measured_strain_temp);


  }


  return h_measured_strain;
}

double A_amplitude(double r_Mpc, double m1, double m2, double tau){
  //r_Mpc is the distance to the source in Mpc
  //m1 and m2 are the masses of the binary compact objects in solar masses
  //Tau is time to coalescence in seconds
  return 3.4e-23 * ( 100. / r_Mpc ) * pow(chirp_mass_source(m1, m2), 5./4.) * pow(tau, -1./4.);
}


double phase_2PN(double Theta, double phi0, double nu){
  return phi0 - pow(Theta, 5./8.) / nu * (1. + (55./96. * nu + 3715./8064.) * pow(Theta, - 1./4.) + (-3./4. * M_PI) * pow(Theta, - 3./8.) + (1855./2048. * nu * nu + 284875./258048. * nu + 9275495./14450688.) * pow(Theta, -1./2.)  );
}

double Theta_val(double t, double tc, double nu, double m){
  return (speed_of_light_cubed_over_G_over_solar_mass * nu / 5.) * (tc - t);
}

double mu_source(double m1, double m2){
  return m1 * m2 / (m1 + m2);
}

double nu_source(double m1, double m2){
  return m1 * m2 / pow(m1 + m2, 2.);
}

double m_source(double m1, double m2){
  return m1 + m2;
}

double chirp_mass_source(double m1, double m2){
  return pow(m1 * m2, 3./5.) / pow(m1 + m2, 1./5.);
}

//Write vector of doubles to file
void write_vector_to_file(vector<double> v, string filename){
  ofstream file;
  file.open(filename);
  for(int i = 0; i < v.size(); i++){
    file << v[i] << endl;
  }
  file.close();
}