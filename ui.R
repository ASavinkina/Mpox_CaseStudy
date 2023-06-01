#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(shinydashboard)

# Define UI for application that draws a histogram
shinyUI(fluidPage(

    # Application title
    titlePanel("University Monkeypox Model"),

    # Sidebar with a slider input for number of bins
    sidebarLayout(
        sidebarPanel(
            
          h3("Disease inputs"),  
          sliderInput("initial_inf",
                      "Infected students at start of semester",
                      min = 0,
                      max = 10,
                      value = 1),
          
          sliderInput("exograte",
                      "Average number of infections coming in from outside of university community in 100 days",
                      min = 0,
                      max = 5,
                      value = 1),
          
          sliderInput("R0_h",
                      "R0 in high-risk subgroup",
                      min = 1,
                      max = 5,
                      value = 2.4, 
                      step=0.2),
          
          sliderInput("R0_l",
                      "R0 in low-risk subgroup",
                      min = 0.5,
                      max = 2.0,
                      value = 0.8,
                      step=0.1),
          
          h3("University inputs"),
            sliderInput("Yale_undegrad_pop",
                        "Undergraduate student population:",
                        min = 1000,
                        max = 20000,
                        value = 6500),
          
          sliderInput("HR_prop",
                      "Proportion of university population in high-risk subgroup:",
                      min = 0.01,
                      max = 1,
                      value = 0.1),
          
          sliderInput("propVax",
                      "Proportion of high-risk sub-group already vaccinated:",
                      min = 0,
                      max = 1,
                      value = 0),
      
            
            sliderInput("studentcontacts",
                        "Vaccinated contacts per diagnosed student",
                        min = 0,
                        max = 20,
                        value = 0),
            
            sliderInput("diagrate",
                        "Proportion of infections diagnosed and isolated",
                        min = 0,
                        max = 1,
                        value = 0,
                        step=0.05),
          
          sliderInput("vaxefficacy",
                      "Efficacy of vaccine",
                      min = 0,
                      max = 1,
                      value = 0.8),
          
          sliderInput("isoduration",
                      "Duration of isolation",
                      min = 21,
                      max = 39,
                      value = 28,
                      step=1),
          
          # sliderInput("vaxduration",
          #             "Duration of quarantine",
          #             min = 8,
          #             max = 21,
          #             value = 14,
          #             step=1),
            
          h3("Cost inputs"), 
          sliderInput("costvax",
                      "Cost of vaccination",
                      min = 0,
                      max = 1000,
                      value = 50,
                      step=50),
          
          sliderInput("costdetect",
                      "Cost of isolation, per student per day",
                      min = 0,
                      max = 5000,
                      value = 200,
                      step=100),
          

            
            submitButton(text = "Apply Changes", icon = NULL, width = NULL)),
        
        
        


        # Show a plot of the generated distribution
        mainPanel(
            #tabsetPanel(type="tabs",
                       # tabPanel("Summary results",
                                    h3("Summary results over 100 days"),
                                    fluidRow(splitLayout(cellWidths = c("33%", "33%","33%"),
                                                         plotOutput("maxinfectionsplot"),
                                                         plotOutput("maxisoplot"),
                                                         plotOutput("maxvaxplot"), width=8)),
                                    h3("Outbreak summary results over 100 days"),
                                    fluidRow(valueBoxOutput("nonewinfections_hr"), valueBoxOutput("meaninfections_hr"),
                                             valueBoxOutput("meaninfections")),
                                    fluidRow(valueBoxOutput("nonewinfections_lr"), valueBoxOutput("meaninfections_lr")),
                                    # fluidRow(valueBoxOutput("isocaplikelihood"),
                                    #          plotOutput("maxisoplot")),
                                    # fluidRow(valueBoxOutput("quarcaplikelihood"),
                                    #          plotOutput("maxquarplot")),
                                    h3("Cost output"),
                                    fluidRow(valueBoxOutput("meanvaxcost"), valueBoxOutput("meandetectcost")),
                                    h3("Stochastic results over time"),
                                    fluidRow(splitLayout(cellWidths = c("33%", "33%","33%"), 
                                                         plotOutput("DPlot"),
                                                         plotOutput("QPlot"), 
                                                         plotOutput("IPlot"), width=8)),
                        #tabPanel("Stochastic results",fluidRow(splitLayout(cellWidths = c("50%", "50%"), plotOutput("D1Plot"),plotOutput("Q1Plot"))),
                                    #fluidRow(splitLayout(cellWidths = c("50%", "50%"), plotOutput("I1Plot"),plotOutput("R1Plot"))),

                                    plotOutput("I_plot")
        ))
  #  )
))


