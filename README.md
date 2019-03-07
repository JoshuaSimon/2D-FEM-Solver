# 2D-FEM-Solver
Implementation of a 2D FEM solver for elastic plain stress problems. 


This FEM solver is able to read a specific input file, which contains all the mesh data. There is no limit to number of nodes or anything. The only thing to be ensured is that the input file has the right structure. An example for an input file is given within this repository. Form there the mathematical model is created and all calculations and solving is done. The 2D solvers comes with a lightweight post processor. Visualizing the mesh, displacement and stress output is done with "gnuplot". This program does not include a mesh generator. Creating a geometry and a meshed geometry must be done seperatly. 



The C++ template library for linear algebra "Eigen" is used in this project. You can download it for free here:
http://eigen.tuxfamily.org/index.php?title=Main_Page

All plots were created with "gnuplot". For more information see:
http://www.gnuplot.info/

Geometry and input mesh data were created with the "GiD" software. For more information see:
https://www.gidhome.com/


[![License: CC BY-SA 4.0](https://licensebuttons.net/l/by-sa/4.0/80x15.png)](https://creativecommons.org/licenses/by-sa/4.0/)

All files in this git are published under a Creative Commons Attribution-ShareAlike 4.0 International license (CC BY-SA 4.0).

https://creativecommons.org/licenses/by-sa/4.0/legalcode
