#ifndef Waveform_h
#define Waveform_h

#include "Waveform.h"
#include <vector>
#include <string>

using namespace std;


const double t_max = 0.2; //s
//const double deltat = 1./4096.; //s, sampling time
const double deltat = 1./(4096*32); //s, sampling time


const double speed_of_light_cubed_over_G_over_solar_mass = 203000.; //1/s, with masses expressed in solar masses

double A_amplitude(double, double, double, double);
double phase_2PN(double, double, double);
double Theta_val(double, double, double, double);
double mu_source(double, double);
double nu_source(double, double);
double m_source(double, double);
double chirp_mass_source(double, double);
void write_vector_to_file(vector<double>, string);

vector<double> h_measured_strain(double, double, double, double, double, double, double, double);

#endif
