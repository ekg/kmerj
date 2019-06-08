#pragma once

//#include <memory>
#include <cassert>
#include <iostream>
#include <fstream>
//#include <strstream>
#include <sstream>
#include <cmath>

namespace kmerj {

namespace melting {

// adapted from https://github.com/Sely85/dna_melting

double wallace_rule(double sequence_length, double a_count , double c_count, double g_count, double t_count);
double salt(double salt_conc, double a_count , double c_count, double g_count, double t_count);
double khandelwal(int sequence_length, const std::string& specie, double salt_conc, double dna_conc);
double bre_nearest_neighbor(int sequence_length, const std::string& specie, double salt_conc, double dna_conc, double a_count, double c_count, double g_count, double t_count);
double san_nearest_neighbor(int sequence_length, const std::string& specie, double salt_conc, double dna_conc, double a_count, double c_count, double g_count, double t_count);
double sug_nearest_neighbor(int sequence_length, const std::string& specie, double salt_conc, double dna_conc, double a_count, double c_count, double g_count, double t_count);
double consensus(int sequence_length, const std::string& specie, double salt_conc, double dna_conc, double a_count, double c_count, double g_count, double t_count);
double bre_enthalpy(int sequence_length, const std::string& specie);
double bre_entropy(int sequence_length, const std::string& specie);
double bre_melting_curve(int sequence_length, const std::string& specie, double dna_conc);
double san_enthalpy(int sequence_length, const std::string& specie);
double san_entropy(int sequence_length, const std::string& specie);
double san_melting_curve(int sequence_length, const std::string& specie, double dna_conc);
double sug_enthalpy(int sequence_length, const std::string& specie);
double sug_entropy(int sequence_length, const std::string& specie);
double sug_melting_curve(int sequence_length, const std::string& specie, double dna_conc);

}

}
