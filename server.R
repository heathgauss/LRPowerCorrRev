# server.R

source("helpers.R")

shinyServer(
  function(input, output) {
    
    output$poweroutput <- renderPrint({
      LRPowerCorr(sampsize = input$n, nsims = input$nsims, p = input$rho, 
                   a = input$a, b = input$b, c = input$c, d = input$d,
                   A = input$A, B = input$B, C = input$C, D = input$D,
                   or4 = input$or4, or5 = input$or5, or6 = input$or6, 
                   or7 = input$or7, or8 = input$or8, or9 = input$or9, or10 = input$or10,
                   fullmodel = input$fullmodel, 
                   reducedmodel = input$reducedmodel, 
                   alpha = input$alpha, 
                   dftest = input$df, pcx1 = input$pcx1, pcx2 = input$pcx2)
    })
    
  }
)
