# ui.R

shinyUI(fluidPage(
  titlePanel("LR Power"),
  
  sidebarLayout(
    sidebarPanel(
      helpText("Estimate power associated with a logistic regression interaction term in the 2 x 2 x 2 case."),
      
      numericInput("n", 
                  label = "Choose a sample size",
                  value = 2600,
                  min = 1),
      
      numericInput("nsims", 
                   label = "Choose the number of simulations",
                   value = 1000,
                   min = 1),
      
      numericInput("rho", 
                   label = "Choose correlation",
                   value = 0,
                   min = -1,
                   max = 1),
      
      numericInput("a", 
                   label = "Enter cell count a",
                   value = 170,
                   min = 0),
      
      numericInput("b", 
                   label = "Enter cell count b",
                   value = 170,
                   min = 0),
      
      numericInput("c", 
                   label = "Enter cell count c",
                   value = 170,
                   min = 0),
      
      numericInput("d", 
                   label = "Enter cell count d",
                   value = 340,
                   min = 0),
      
      numericInput("A", 
                   label = "Enter cell count A",
                   value = 30,
                   min = 0),
      
      numericInput("B", 
                   label = "Enter cell count B",
                   value = 18,
                   min = 0),
      
      numericInput("C", 
                   label = "Enter cell count C",
                   value = 30,
                   min = 0),
      
      numericInput("D", 
                   label = "Enter cell count D",
                   value = 72,
                   min = 0),
            
      numericInput("or4", 
                   label = "Choose OR4",
                   value = 1,
                   min = 0),
      
      numericInput("or5", 
                   label = "Choose OR5",
                   value = 1,
                   min = 0),
      
      numericInput("or6", 
                   label = "Choose OR6",
                   value = 1,
                   min = 0),
      
      numericInput("or7", 
                   label = "Choose OR7",
                   value = 1,
                   min = 0),
      
      numericInput("or8", 
                   label = "Choose OR8",
                   value = 1,
                   min = 0),
      
      numericInput("or9", 
                   label = "Choose OR9",
                   value = 1,
                   min = 0),
      
      numericInput("or10", 
                   label = "Choose OR10",
                   value = 1,
                   min = 0),
      
      textInput("fullmodel", 
                label = "Enter full model",
                value = "y ~ cx1 + cx2 + cx3"),
      
      textInput("reducedmodel", 
                label = "Enter reduced model",
                value = "y ~ cx1 + cx2"),
      
      numericInput("alpha", 
                   label = "Specify alpha",
                   value = 0.05,
                   min = 0,
                   max = 1),
      
      numericInput("df", 
                   label = "Specify degrees of freedom for LR test",
                   value = 1,
                   min = 1),
      
      numericInput("pcx1", 
                   label = "Specify proportion for cx1",
                   value = 0.6,
                   min = 0,
                   max = 1),
      
      numericInput("pcx2", 
                   label = "Specify proportion for cx2",
                   value = 0.15,
                   min = 0,
                   max = 1)
      ),
    
    mainPanel(textOutput("poweroutput"))
  )
))
