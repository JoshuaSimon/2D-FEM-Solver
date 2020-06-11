#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <Rcpp.h>
#include <D:\Dokumente\C++\Includes\eigen\Eigen/Dense>
#include <D:\Dokumente\C++\Includes\eigen\Eigen/Sparse>

//Element data type
struct Element
{
	void CalculateStiffnessMatrix(const Eigen::Matrix3f& D, std::vector<Eigen::Triplet<float> >& triplets);

	Eigen::Matrix<float, 3, 6> B;
	int nodesIds[3];
};

//Boundary constraint data type
struct Constraint
{
	enum Type
	{
		UX = 1 << 0,
		UY = 1 << 1,
		UXY = UX | UY
	};
	int node;
	Type type;
};

//Globals
int				            nodesCount;
Eigen::VectorXf			    nodesX;
Eigen::VectorXf			    nodesY;
Eigen::VectorXf			    loads;
std::vector< Element >		elements;
std::vector< Constraint >	constraints;

//Function for calculating the element stiffness matrix.
void Element::CalculateStiffnessMatrix(const Eigen::Matrix3f& D, std::vector<Eigen::Triplet<float> >& triplets)
{
	Eigen::Vector3f x, y;
	x << nodesX[nodesIds[0]], nodesX[nodesIds[1]], nodesX[nodesIds[2]];
	y << nodesY[nodesIds[0]], nodesY[nodesIds[1]], nodesY[nodesIds[2]];

	Eigen::Matrix3f C;
	C << Eigen::Vector3f(1.0f, 1.0f, 1.0f), x, y;

	//Calculating coefficients for shape functions (a1, a2, a3). 
	//These are relevant for interpolation.
	Eigen::Matrix3f IC = C.inverse();

	//Assemble B matrix
	for (int i = 0; i < 3; i++)
	{
		B(0, 2 * i + 0) = IC(1, i);
		B(0, 2 * i + 1) = 0.0f;
		B(1, 2 * i + 0) = 0.0f;
		B(1, 2 * i + 1) = IC(2, i);
		B(2, 2 * i + 0) = IC(2, i);
		B(2, 2 * i + 1) = IC(1, i);
	}

	//Calculate element stiffness (det(C)/2 = area of triangle).
	Eigen::Matrix<float, 6, 6> K = B.transpose() * D * B * C.determinant() / 2.0f;

	//Store values of element stiffness matrix with corresponding indices in global stiffness matrix in triplets.
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			Eigen::Triplet<float> trplt11(2 * nodesIds[i] + 0, 2 * nodesIds[j] + 0, K(2 * i + 0, 2 * j + 0));
			Eigen::Triplet<float> trplt12(2 * nodesIds[i] + 0, 2 * nodesIds[j] + 1, K(2 * i + 0, 2 * j + 1));
			Eigen::Triplet<float> trplt21(2 * nodesIds[i] + 1, 2 * nodesIds[j] + 0, K(2 * i + 1, 2 * j + 0));
			Eigen::Triplet<float> trplt22(2 * nodesIds[i] + 1, 2 * nodesIds[j] + 1, K(2 * i + 1, 2 * j + 1));

			triplets.push_back(trplt11);
			triplets.push_back(trplt12);
			triplets.push_back(trplt21);
			triplets.push_back(trplt22);
		}
	}
}

//Function for setting constraints. 
void SetConstraints(Eigen::SparseMatrix<float>::InnerIterator& it, int index)
{
	if (it.row() == index || it.col() == index)
	{
		it.valueRef() = it.row() == it.col() ? 1.0f : 0.0f;
	}
}

//Function for applying constraints in stiffness matrix.
void ApplyConstraints(Eigen::SparseMatrix<float>& K, const std::vector<Constraint>& constraints)
{
	std::vector<int> indicesToConstraint;

	for (std::vector<Constraint>::const_iterator it = constraints.begin(); it != constraints.end(); ++it)
	{
		if (it->type & Constraint::UX)
		{
			indicesToConstraint.push_back(2 * it->node + 0);
		}
		if (it->type & Constraint::UY)
		{
			indicesToConstraint.push_back(2 * it->node + 1);
		}
	}

	for (int k = 0; k < K.outerSize(); ++k)
	{
		for (Eigen::SparseMatrix<float>::InnerIterator it(K, k); it; ++it)
		{
			for (std::vector<int>::iterator idit = indicesToConstraint.begin(); idit != indicesToConstraint.end(); ++idit)
			{
				SetConstraints(it, *idit);
			}
		}
	}
}

// Functions
void generateSolverInput(std::string mesh_data_file_name, std::string solver_input_file_name);

// [[Rcpp::export]]
//std::vector<float> solver(std::string mesh_data_file) {
std::vector<std::vector<float>> solver(std::string mesh_data_file) {

    

    Rcpp::Rcout << "---------------------- 2D FEM SOLVER ---------------------" << std::endl;
	Rcpp::Rcout << "FEM software for solving elastic 2D plain stress problems." << std::endl;
	Rcpp::Rcout << "Created by Joshua Simon. Date: 25.02.2019." << std::endl;
	Rcpp::Rcout << std::endl << std::endl;


	//1. Paths and filenames
							
	std::string solver_input_file = "Solver_Input.txt";
	std::string displacement_plot_data = "GNUPlot_Input_displacement.txt";
	std::string stress_plot_data = "GNUPlot_Input_contour.txt";

	//2. Pre Processing:
	//Read GiD mesh Data and write solver input.
	Rcpp::Rcout << "Pre Processor: Define boundary conditions and loads." << std::endl << std::endl;
	generateSolverInput(mesh_data_file, solver_input_file);
	Rcpp::Rcout << "Pre Processor: Solver input generated!" << std::endl << std::endl;

	std::ifstream infile(solver_input_file);                    // This filename is read from user input.
	//std::ofstream outfile("Solution_Data.txt");					// This file contains solution data.
	//std::ofstream outfile_gnuplot(displacement_plot_data);		// This file contains plot data.
	//std::ofstream outfile_gnuplot_contur(stress_plot_data);		// This file contains plot data.

	//3. Solution:
	Rcpp::Rcout << "Solver: Creating mathematical model..." << std::endl;

	//Read material specifications
	float poissonRatio, youngModulus;
	infile >> poissonRatio >> youngModulus;

	//Assemble elasticity matrix D
	Eigen::Matrix3f D;
	D <<
		1.0f, poissonRatio, 0.0f,
		poissonRatio, 1.0, 0.0f,
		0.0f, 0.0f, (1.0f - poissonRatio) / 2.0f;

	D *= youngModulus / (1.0f - pow(poissonRatio, 2.0f));

	//Read number of nodes and their coordinates
	infile >> nodesCount;
	nodesX.resize(nodesCount);
	nodesY.resize(nodesCount);

	for (int i = 0; i < nodesCount; ++i)
	{
		infile >> nodesX[i] >> nodesY[i];
	}

	//Read number of elements and their nodes
	int elementCount;
	infile >> elementCount;

	for (int i = 0; i < elementCount; ++i)
	{
		Element element;
		infile >> element.nodesIds[0] >> element.nodesIds[1] >> element.nodesIds[2];
		elements.push_back(element);
	}

	//Read number of constraints and their node settings
	int constraintCount;
	infile >> constraintCount;

	for (int i = 0; i < constraintCount; ++i)
	{
		Constraint constraint;
		int type;
		infile >> constraint.node >> type;
		constraint.type = static_cast<Constraint::Type>(type);
		constraints.push_back(constraint);
	}

	loads.resize(2 * nodesCount);
	loads.setZero();

	//Read number of nodal loads and nodal forces
	int loadsCount;
	infile >> loadsCount;

	for (int i = 0; i < loadsCount; ++i)
	{
		int node;
		float x, y;
		infile >> node >> x >> y;
		loads[2 * node + 0] = x;
		loads[2 * node + 1] = y;
	}

	//Calculate stiffness matrix for each element
	std::vector<Eigen::Triplet<float> > triplets;
	for (std::vector<Element>::iterator it = elements.begin(); it != elements.end(); ++it)
	{
		it->CalculateStiffnessMatrix(D, triplets);
	}

	//Assemble global stiffness matirx
	Eigen::SparseMatrix<float> globalK(2 * nodesCount, 2 * nodesCount);
	globalK.setFromTriplets(triplets.begin(), triplets.end());

	//Apply Constraints
	ApplyConstraints(globalK, constraints);

	Rcpp::Rcout << "Solver: Mathematical model created!" << std::endl;

	//Solving
	Rcpp::Rcout << "Solver: Solving in progress..." << std::endl;
	Eigen::SimplicialLDLT<Eigen::SparseMatrix<float> > solver(globalK);
	Eigen::VectorXf displacements = solver.solve(loads);
	Rcpp::Rcout << "Solver: Solving done!" << std::endl << std::endl;

	//Writing output and display on console
	//std::cout << "Loads vector:" << std::endl << loads << std::endl << std::endl;			//Loads
	//std::cout << "Displacements vector:" << std::endl << displacements << std::endl;		//Displaysments


	//std::cout << "Stresses:" << std::endl;							//Von Mises Stress

	int m = 0;
	float sigma_max = 0.0;
	//float *sigma_mises = new float[elementCount];
    std::vector<float> output_sigma_mises (elementCount); 

	for (std::vector<Element>::iterator it = elements.begin(); it != elements.end(); ++it)
	{
		Eigen::Matrix<float, 6, 1> delta;
		delta << displacements.segment<2>(2 * it->nodesIds[0]),
			displacements.segment<2>(2 * it->nodesIds[1]),
			displacements.segment<2>(2 * it->nodesIds[2]);

		Eigen::Vector3f sigma = D * it->B * delta;
		output_sigma_mises[m] = sqrt(sigma[0] * sigma[0] - sigma[0] * sigma[1] + sigma[1] * sigma[1] + 3.0f * sigma[2] * sigma[2]);

		//Search for maximum stress
		if (output_sigma_mises[m] > sigma_max) {
			sigma_max = output_sigma_mises[m];
		}

		m++;
	}

	//4. Post Processing:
    //std::vector<std::vector<float>> output_displacement(elementCount*4, std::vector<float> (4));
	std::vector<std::vector<float>> output_displacement(elementCount*3, std::vector<float> (4));
    int output_index = 0;
    //4.1 Writing GNUPlot output file for ploting mesh and mesh + displacements
	for (std::vector<Element>::iterator it = elements.begin(); it != elements.end(); ++it)
	{
		//Prints x,y,dis-x,dis-y for every node of element in one line
		for (int i = 0; i < 3; i++) {
			output_displacement[output_index][0] = nodesX(it->nodesIds[i]);
            output_displacement[output_index][1] = nodesY(it->nodesIds[i]);
		    output_displacement[output_index][2] = displacements(it->nodesIds[i] * 2);
            output_displacement[output_index][3] = displacements(it->nodesIds[i] * 2 + 1);

            output_index++;
		}

		// First node of element has to appear twice for plotting purpose
        // to close the triangle. 
		/*output_displacement[output_index][0] = nodesX(it->nodesIds[0]);
        output_displacement[output_index][1] = nodesY(it->nodesIds[0]);
	    output_displacement[output_index][2] = displacements(it->nodesIds[0] * 2);
        output_displacement[output_index][3] = displacements(it->nodesIds[0] * 2 + 1);
		
        output_index++;*/
	}

    Rcpp::Rcout << std::endl;
    //return output_sigma_mises;
    return output_displacement;
}


// Read mesh data (GiD export file) and create solver input
void generateSolverInput(std::string mesh_data_file_name, std::string solver_input_file_name)
{
	std::ifstream infile(mesh_data_file_name);
	std::ofstream outfile(solver_input_file_name);
	
	//Read general information
	int problemDimension, nodesCount, elementCount, materialCount, BoundaryEdgesCount, BoundaryBlocksCount;
	infile >> problemDimension;
	infile >> nodesCount >> elementCount >> materialCount >> BoundaryEdgesCount >> BoundaryBlocksCount;

	//Wirte matiral specifcations (Possion's Ratio and E-Modulus)
	outfile << 0.3 << " " << 2000.0 << std::endl;

	//Write total number of nodes
	outfile << nodesCount << std::endl;

	//Read and write x,y coordinates of the nodes
	int nodeNumber;
	double x, y;

	for (int i = 0; i < nodesCount; i++) {
		infile >> nodeNumber >> x >> y;
		outfile << x << "  " << y << std::endl;
	}

	//Read type and number of Elements
	int type, elementNumber;
	infile >> type >> elementNumber;
	outfile << elementNumber << std::endl;

	//Read, format and wirte nodes to corresponding elements
	int elementMaterial, node_1, node_2, node_3;

	for (int i = 0; i < elementNumber; i++) {
		infile >> elementMaterial >> node_1 >> node_2 >> node_3;
		node_1 -= 1;
		node_2 -= 1;
		node_3 -= 1;
		outfile << node_1 << " " << node_2 << " " << node_3 << std::endl;
	}

	//Read boundary conditions
	int boundType, boundEdgeNumbers, lineSpecifier;
	int fixedNodesCount = 0;
	int forceNodesCount = 0;
	int k = 0;
	double fx, fy;
	int *nodeSpecifier = new int[BoundaryEdgesCount];
	int *lines = new int[BoundaryEdgesCount];
	int *start_nodes = new int[BoundaryEdgesCount];
	int *end_nodes = new int[BoundaryEdgesCount];
	
	for (int i = 0; i < BoundaryBlocksCount; i++) {

		infile >> boundType >> boundEdgeNumbers;

		//User input to specify boundarys
		Rcpp::Rcout << "Choose condition for boundary line " << i+1 << " [(1): Fixed Support or (2): Force]:" << std::endl;
		Rcpp::Rcout << ">> "; 
        std::cin >> lineSpecifier;
		Rcpp::Rcout << std::endl;
		
		//Input nodal force for x and y components
		if (lineSpecifier == 2) {
			Rcpp::Rcout << "Nodal force fx = ";
			std::cin >> fx;
			Rcpp::Rcout << "Nodal force fy = ";
			std::cin >> fy;
			Rcpp::Rcout << std::endl;
		}

		for (int j = 0; j < boundEdgeNumbers; j++) {
			infile >> lines[k] >> start_nodes[k] >> end_nodes[k];
			start_nodes[k] -= 1;
			end_nodes[k] -= 1;
			nodeSpecifier[k] = lineSpecifier;

			if (lineSpecifier == 1) {
				fixedNodesCount++;
			}

			if (lineSpecifier == 2) {
				forceNodesCount++;
			}

			k++;
		}
	}

	//Write number of fixed nodes to solver input file
	outfile << fixedNodesCount << std::endl;
	
	//Write boundary node and condition type to solver input file (1 = UX fixed, 2 = UY fixed, 3 = UX and UY fixed).
	for (int i = 0; i < BoundaryEdgesCount; i++) {
		if (nodeSpecifier[i] == 1) {
			outfile << start_nodes[i] << " " << 3 << std::endl;
		}
	}
	
	//Write number of force nodes to solver input file
	outfile << forceNodesCount << std::endl;

	//Write force node and nodal force componets to solver input file
	for (int i = 0; i < BoundaryEdgesCount; i++) {
		if (nodeSpecifier[i] == 2) {
			outfile << start_nodes[i] << " " << fx << " " << fy << std::endl;
		}
	}

	delete[] nodeSpecifier;
	delete[] lines;
	delete[] start_nodes;
	delete[] end_nodes;

	return;
}