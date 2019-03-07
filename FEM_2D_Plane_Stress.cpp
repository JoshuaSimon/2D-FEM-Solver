// FEM_2D_Plane_Stress.cpp (Main)

#include "pch.h"
#include "FEM_Input.h"
#include "FEM_GNUPlot.h"
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <string>
#include <vector>
#include <iostream>
#include <fstream>

struct Element
{
	void CalculateStiffnessMatrix(const Eigen::Matrix3f& D, std::vector<Eigen::Triplet<float> >& triplets);

	Eigen::Matrix<float, 3, 6> B;
	int nodesIds[3];
};

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

int				nodesCount;
Eigen::VectorXf			nodesX;
Eigen::VectorXf			nodesY;
Eigen::VectorXf			loads;
std::vector< Element >		elements;
std::vector< Constraint >	constraints;

void Element::CalculateStiffnessMatrix(const Eigen::Matrix3f& D, std::vector<Eigen::Triplet<float> >& triplets)
{
	Eigen::Vector3f x, y;
	x << nodesX[nodesIds[0]], nodesX[nodesIds[1]], nodesX[nodesIds[2]];
	y << nodesY[nodesIds[0]], nodesY[nodesIds[1]], nodesY[nodesIds[2]];

	Eigen::Matrix3f C;
	C << Eigen::Vector3f(1.0f, 1.0f, 1.0f), x, y;

	Eigen::Matrix3f IC = C.inverse();

	for (int i = 0; i < 3; i++)
	{
		B(0, 2 * i + 0) = IC(1, i);
		B(0, 2 * i + 1) = 0.0f;
		B(1, 2 * i + 0) = 0.0f;
		B(1, 2 * i + 1) = IC(2, i);
		B(2, 2 * i + 0) = IC(2, i);
		B(2, 2 * i + 1) = IC(1, i);
	}
	Eigen::Matrix<float, 6, 6> K = B.transpose() * D * B * C.determinant() / 2.0f;

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

void SetConstraints(Eigen::SparseMatrix<float>::InnerIterator& it, int index)
{
	if (it.row() == index || it.col() == index)
	{
		it.valueRef() = it.row() == it.col() ? 1.0f : 0.0f;
	}
}

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

int main(void)
{
	cout << "---------------------- 2D FEM SOLVER ---------------------" << std::endl;
	cout << "FEM software for solving elastic 2D plain stress problems." << std::endl;
	cout << "Created by Joshua Simon, Date: 25.05.2019." << std::endl;
	cout << std::endl;

	//1. Pre Processing
	//Read GiD Mesh Data and write solver input.
	cout << "Pre Processor: Define boundary conditions and loads." << std::endl << std::endl;
	generateSolverInput("Mesh.txt");
	cout << "Pre Processor: Solver input generated!" << endl << endl;

	std::ifstream infile("Solver_Input.txt");
	std::ofstream outfile("Outputfile_2.txt");
	std::ofstream outfile_gnuplot("GNUPlot_Input.txt");
	std::ofstream outfile_gnuplot_contur("GNUPlot_Input_contur.txt");

	//2. Solution
	std::cout << "Solver: Creating mathematical model..." << std::endl;

	//Read material specifications
	float poissonRatio, youngModulus;
	infile >> poissonRatio >> youngModulus;


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

	//Calculate stiffness matrix
	std::vector<Eigen::Triplet<float> > triplets;
	for (std::vector<Element>::iterator it = elements.begin(); it != elements.end(); ++it)
	{
		it->CalculateStiffnessMatrix(D, triplets);
	}

	Eigen::SparseMatrix<float> globalK(2 * nodesCount, 2 * nodesCount);
	globalK.setFromTriplets(triplets.begin(), triplets.end());

	//Apply Constraints
	ApplyConstraints(globalK, constraints);

	std::cout << "Solver: Mathematical model created!" << std::endl;

	//Solving
	std::cout << "Solver: Solving in progress..." << std::endl;
	Eigen::SimplicialLDLT<Eigen::SparseMatrix<float> > solver(globalK);
	Eigen::VectorXf displacements = solver.solve(loads);
	std::cout << "Solver: Solving done!" << std::endl << std::endl;

	//Writing output and display on consol
	//std::cout << "Loads vector:" << std::endl << loads << std::endl << std::endl;			//Loads
	//std::cout << "Displacements vector:" << std::endl << displacements << std::endl;		//Displaysments

	outfile << displacements << std::endl;

	//std::cout << "Stresses:" << std::endl;													//Von Mises Stress

	int m = 0;
	float sigma_max = 0.0;
	float *sigma_mises = new float[elementCount];

	for (std::vector<Element>::iterator it = elements.begin(); it != elements.end(); ++it)
	{
		Eigen::Matrix<float, 6, 1> delta;
		delta << displacements.segment<2>(2 * it->nodesIds[0]),
			displacements.segment<2>(2 * it->nodesIds[1]),
			displacements.segment<2>(2 * it->nodesIds[2]);

		Eigen::Vector3f sigma = D * it->B * delta;
		//float sigma_mises = sqrt(sigma[0] * sigma[0] - sigma[0] * sigma[1] + sigma[1] * sigma[1] + 3.0f * sigma[2] * sigma[2]);
		sigma_mises[m] = sqrt(sigma[0] * sigma[0] - sigma[0] * sigma[1] + sigma[1] * sigma[1] + 3.0f * sigma[2] * sigma[2]);

		//Search for maximum stress
		if (sigma_mises[m] > sigma_max) {
			sigma_max = sigma_mises[m];
		}

		//std::cout << sigma_mises << std::endl;
		outfile << sigma_mises[m] << std::endl;

		m++;
	}

	//3. Post Processing
	//3.1 Wirting GNUPlot output file for ploting mesh and mesh + displacements
	for (std::vector<Element>::iterator it = elements.begin(); it != elements.end(); ++it)
	{
		//Prints x,y,dis-x,dis-y for every node of element in one line
		for (int i = 0; i < 3; i++) {
			outfile_gnuplot << nodesX(it->nodesIds[i]) << " " << nodesY(it->nodesIds[i]) \
				<< " " << displacements(it->nodesIds[i] * 2) << " " << displacements(it->nodesIds[i] * 2 + 1) << std::endl;
		}
		//First node of element has to appear twice for plotting purpose
		outfile_gnuplot << nodesX(it->nodesIds[0]) << " " << nodesY(it->nodesIds[0]) \
			<< " " << displacements(it->nodesIds[0] * 2) << " " << displacements(it->nodesIds[0] * 2 + 1) << std::endl;

		//Empty line to sperate between elements
		outfile_gnuplot << std::endl;
	}



	//3.2 Wirting GNUPlot output file for stress conturplot
	outfile_gnuplot_contur << "unset xtics" << std::endl;
	outfile_gnuplot_contur << "unset ytics" << std::endl;
	outfile_gnuplot_contur << "set cbrange [0:1]" << std::endl << std::endl;
	outfile_gnuplot_contur << "plot[-15:15][-15:15] \\" << std::endl;

	//Wirte color information for every element
	int mm = 0;
	for (std::vector<Element>::iterator it = elements.begin(); it != (elements.end()-1); ++it) {
		outfile_gnuplot_contur << "\"-\" title \"\" with filledcurve lt palette cb " \
			                   << sigma_mises[mm] / sigma_max << " \\" << std::endl;
		outfile_gnuplot_contur << "fillstyle transparent solid 1.000000 ,\\" << std::endl;
		mm++;
	}

	//Wirte color information for last element
	outfile_gnuplot_contur << "\"-\" title \"\" with filledcurve lt palette cb " \
						   << sigma_mises[elementCount-1] / sigma_max << " \\" << std::endl;
	outfile_gnuplot_contur << "fillstyle transparent solid 1.000000 ;" << std::endl;


	for (std::vector<Element>::iterator it = elements.begin(); it != elements.end(); ++it)
	{
		//Elements and thier nodes:
		//Prints x,y for every node of element in one line
		for (int i = 0; i < 3; i++) {
			outfile_gnuplot_contur << nodesX(it->nodesIds[i]) << " " << nodesY(it->nodesIds[i]) << std::endl;
		}
		//First node of element has to appear twice for plotting purpose
		outfile_gnuplot_contur << nodesX(it->nodesIds[0]) << " " << nodesY(it->nodesIds[0]) << std::endl;

		//'e' to sperate between elements
		outfile_gnuplot_contur << "e" << std::endl;
	}

	std::cout << "Post Processor: Plotting solution..." << std::endl << std::endl;

	//3.3 Plot results
	plot("D:\\Dokumente\\C++\\FEM_2D_Plane_Stress\\Debug\\GNUPlot_Input.txt");

	delete[] sigma_mises;

	return 0;
}
