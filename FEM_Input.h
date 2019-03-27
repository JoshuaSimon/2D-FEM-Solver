// FEM_Input.h
#pragma once

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
using namespace std;

// Function for reading mesh data(GiD export file) and creating solver input
void generateSolverInput(string mesh_data_file_name, string solver_input_file_name);
