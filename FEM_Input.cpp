// FEM_Input.cpp
#include "pch.h"
#include "FEM_Input.h"

void generateSolverInput(string filename)
{
	ifstream infile(filename);
	ofstream outfile("D:\\Dokumente\\C++\\FEM_2D_Plane_Stress\\Debug\\Solver_Input.txt");
	
	//Read general information
	int problemDimension, nodesCount, elementCount, materialCount, BoundaryEdgesCount, BoundaryBlocksCount;
	infile >> problemDimension;
	infile >> nodesCount >> elementCount >> materialCount >> BoundaryEdgesCount >> BoundaryBlocksCount;

	//Wirte matiral specifcations (Possion's Ratio and E-Modulus)
	outfile << 0.3 << " " << 2000.0 << endl;

	//Write total number of nodes
	outfile << nodesCount << endl;

	//Read and write x,y coordinates of the nodes
	int nodeNumber;
	double x, y;

	for (int i = 0; i < nodesCount; i++) {
		infile >> nodeNumber >> x >> y;
		outfile << x << "  " << y << endl;
	}

	//Read type and number of Elements
	int type, elementNumber;
	infile >> type >> elementNumber;
	outfile << elementNumber << endl;

	//Read, format and wirte nodes to coresbonding elements
	int elementMaterial, node_1, node_2, node_3;

	for (int i = 0; i < elementNumber; i++) {
		infile >> elementMaterial >> node_1 >> node_2 >> node_3;
		node_1 -= 1;
		node_2 -= 1;
		node_3 -= 1;
		outfile << node_1 << " " << node_2 << " " << node_3 << endl;
	}

	//Read Boundary Conditions
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
		cout << "Choose condition for boundary line " << i+1 << " [(1): Fixed Support or (2): Force]:" << endl;
		cout << ">> "; cin >> lineSpecifier;
		cout << endl;
		
		//Input nodal force for x and y componets
		if (lineSpecifier == 2) {
			cout << "Nodal force fx = ";
			cin >> fx;
			cout << "Nodal force fy = ";
			cin >> fy;
			cout << endl;
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
	outfile << fixedNodesCount << endl;
	
	//Write boundary node and condition type to solver input file (1 = UX fixed, 2 = UY fixed, 3 = UX and UY fixed).
	for (int i = 0; i < BoundaryEdgesCount; i++) {
		if (nodeSpecifier[i] == 1) {
			outfile << start_nodes[i] << " " << 3 << endl;
		}
	}
	
	//Write number of force nodes to solver input file
	outfile << forceNodesCount << endl;

	//Write force node and nodal force componets to solver input file
	for (int i = 0; i < BoundaryEdgesCount; i++) {
		if (nodeSpecifier[i] == 2) {
			outfile << start_nodes[i] << " " << fx << " " << fy << endl;
		}
	}

	delete[] nodeSpecifier;
	delete[] lines;
	delete[] start_nodes;
	delete[] end_nodes;

	return;
}
