# app.r
# This file contains the server and ui functions of the shiny web app.
#
# For running the app form a R terminal use:
# shiny::runApp("D:\\Dokumente\\GitHub\\2D-FEM-Solver\\R Shiny Webapp\\app.r")

library(shiny)
library(shinythemes)
source("D:\\Dokumente\\GitHub\\2D-FEM-Solver\\R Shiny Webapp\\2D_fem_pre_processor.r", chdir = TRUE)

# Compile the C++ source and run the solver.
init_solver()

# Call C++ functions form the compiled source to calculate all 
# the displacements from the provided input file which contains
# the mesh data. 
input_filename <- "D:\\Dokumente\\GitHub\\2D-FEM-Solver\\R Shiny Webapp\\Mesh_Balken_V2.txt"
displacements = solve_linear_elastic(input_filename)
mises_stresses = output_mises_stress()

# Reshape data.
x <- unname(sapply(displacements, `[[`, 1))
y <- unname(sapply(displacements, `[[`, 2))
ux <- unname(sapply(displacements, `[[`, 3))
uy <- unname(sapply(displacements, `[[`, 4))


# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("Mechanical Linear Elastic 2D FEM Solver"),

    # Color theme of the page.
    theme = shinytheme("slate"),

    fluidRow(
        # Left sidebar panel for inputs.
        column(2,
            fluidRow(
                wellPanel(
                    h4("Solver Input Filename"),
                    textOutput("text_filename")
                )
            ),

            fluidRow(
                wellPanel(
                # Input: Select a timestep.
                sliderInput("bins",
                        "Time:",
                        min = 0,
                        max = 1,
                        value = 1),
            
                # Button for updating plots.
                actionButton("update", "Update Plot"),
                width = 2
                )
            )

        ),

        # Main panel.
        column(10,
            wellPanel(
                h4("Stress and Displacment Plot"),

                # Plot mesh with stresses and displacements.
                plotOutput("distPlot")
            )
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {

    output$text_filename <- renderText({paste(input_filename)})

    output$distPlot <- renderPlot({
        # Generate bins based on input$bins from ui.R
        bins <- input$bins

        # Draw plot in app. 
        mesh_ggplot(x, y, ux*bins, uy*bins, mises_stresses)
    })
}

# Run the application 
shinyApp(ui = ui, server = server)
