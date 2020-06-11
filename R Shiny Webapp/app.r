library(Rcpp)
library(shiny)
library(shinythemes)
library(dplyr)
library(ggplot2)

mesh_rplot <- function(x, y, ux, uy, scale_factor) {
    # Create empty plot.
    plot(-5:5,-5:5,type='n')

    # Add polygons to the plot for each element.
    node_index <- 1
    element_nodes_global <- character()

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
}

mesh_ggplot <- function(x, y, ux, uy, stress) {
    # See https://ggplot2.tidyverse.org/reference/geom_polygon.html for info
    # on polygon plots with ggpolt2. 
    element_count <- length(x)/3

    # Create an unique ID for every element to match the data later.
    ids <- 1:element_count

    # Stress values to color the elements. 
    values <- data.frame(
        id = ids,
        value = stress
    )

    # Coordinates of the nodes. 
    positions <- data.frame(
        id = rep(ids, each = 3),
        xx = x,
        yy = y,
        uxx = x + ux,
        uyy = y + uy
    )

    # Match the data of the stress values and the positions using the ID.
    poly_data <- merge(values, positions, by = c("id"))

    # Plot the mesh data. 
    p <- ggplot() +
        geom_polygon(data = poly_data, aes(x = xx, y = yy, fill = value, group = id), colour = "black") +
        geom_polygon(data = poly_data, aes(x = uxx, y = uyy, fill = value, group = id), colour = "black") +
        xlim(-5, 5) +
        ylim(-5, 5) +
        scale_fill_gradient(low="blue", high="red")

    return(p)
}

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

# Call C++ functions form the compiled source to calculate all 
# the displacements from the provided input file which contains
# the mesh data. 
displacements = solver("D:\\Dokumente\\GitHub\\2D-FEM-Solver\\R Shiny Webapp\\Mesh_Balken_V2.txt")
mises_stresses = output_mises_stress()

# Reshape data.
x <- unname(sapply(displacements, `[[`, 1))
y <- unname(sapply(displacements, `[[`, 2))
ux <- unname(sapply(displacements, `[[`, 3))
uy <- unname(sapply(displacements, `[[`, 4))

# Plot mesh and mesh + displacement with standard r plot. 
mesh_rplot(x, y, ux, uy, 1)

# Plot contur of mesh + stress wwith ggpolt2.
print(mesh_ggplot(x, y, ux, uy, mises_stresses))
