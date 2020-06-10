library(Rcpp)
library(shiny)
library(shinythemes)
library(dplyr)
library(ggplot2)

# Set system variable for the compiler, so that it
# can find external librarys and header files. 
Sys.setenv("PKG_CXXFLAGS"="-I/D:\\Dokumente\\C++\\Includes")

print("Compling files...")

# Start the clock.
ptm <- proc.time()

# Compile the C++ Source file with the main function. 
sourceCpp("D:/Dokumente/GitHub/2D-FEM-Solver/R Shiny Webapp/solver_main.cpp")

# Stop the clock and print it to console. 
print(proc.time() - ptm)
print("Done compling!")

# Call C++ function form the compiled source to calculate all 
# the displacements from the provided input file which contains
# the mesh data. 
displacements = solver("D:\\Dokumente\\GitHub\\2D-FEM-Solver\\R Shiny Webapp\\Mesh_Balken_V2.txt")

# Reshape data.
x <- unname(sapply(displacements, `[[`, 1))
y <- unname(sapply(displacements, `[[`, 2))
ux <- unname(sapply(displacements, `[[`, 3))
uy <- unname(sapply(displacements, `[[`, 4))

plot(x, y, "p")

displacements_data <- data.frame(x, y, ux, uy)

# Plot the mesh of the provided input data. 
#node_plt <- ggplot(displacements_data, aes(x = x, y = y)) + 
#                geom_line(color="blue") +
#                geom_point(color="blue") +
#                xlim(-10, 10) +
#                ylim(-2, 2)
#print(node_plt)

