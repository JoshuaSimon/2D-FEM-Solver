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

# Create empty plot.
plot(-5:5,-5:5,type='n')

# Add polygons to the plot for each element.
node_index <- 1
element_nodes_global <- character()
scale_factor <- 1

for (element in 1:(length(x)/3) ) {
    # Put the node coordinates (xx, yy) to one element (= triangle) in one list.
    # Mesh data.  
    element_nodes_local_xy <- list(xx=c(x[node_index], x[node_index + 1], x[node_index + 2]),
                                yy=c(y[node_index], y[node_index + 1], y[node_index + 2]))

    # Mesh data + scale_factor * displacement.
    element_nodes_local_uxuy <- list(xx=c(x[node_index] + scale_factor * ux[node_index],
                                        x[node_index + 1] + scale_factor * ux[node_index + 1],
                                        x[node_index + 2] + scale_factor * ux[node_index + 2]),
                                    yy=c(y[node_index] + scale_factor * uy[node_index],
                                        y[node_index + 1] + scale_factor * uy[node_index + 1],
                                        y[node_index + 2] + scale_factor * uy[node_index + 2]))

    node_index <- node_index + 3

    # Draw pologons for each element. First mesh, then mesh + displacement. 
    polygon(element_nodes_local_xy$xx, element_nodes_local_xy$yy, col="gray", border="blue")
    polygon(element_nodes_local_uxuy$xx, element_nodes_local_uxuy$yy, col="gray", border="red")
}



#displacements_data <- data.frame(x, y, ux, uy)
#print(head(displacements_data))

# Plot the mesh of the provided input data. 
#node_plt <- ggplot(displacements_data, aes(x = x, y = y)) + 
#                geom_line(color="blue") +
#                geom_point(color="blue") +
#                xlim(-10, 10) +
#                ylim(-2, 2)
#print(node_plt)

