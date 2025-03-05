################################################################
################################################################
## FUNCTIONS FOR PLOTTING MATRICES:

#################################################
#################################################
##  plot.pdm:
##    Function that generates a plot of a square
##    matrix.
##
##  INPUTS
##    pd          : square matrix to plot. 
##    main        : title of the plot.
##    breaks      : number of colors in the gradient
##                    used for plotting.
##    diagonal    : if FALSE, plot only displays 
##                    non-diagonal entries.
##    cex         : size of points in the plot.
##    entryrange  : if NULL, calculated automatically.
##                    Otherwise, bound on the 
##                    absolute values of entries plotted.
##    
##  OUTPUT:
##    Color plot of matrix with gradient color from blue to red
##    indicating sign and magnitude of entries.
##
plot.pdm <- function(pd,
                     main,
                     breaks = 20,
                     diagonal = FALSE,
                     cex = 0.05,
                     entryrange = NULL) {
  .d <- dim(pd)[1]
  
  if (is.null(entryrange)) {
    entryrange <- 0
    if (diagonal) {
      entryrange <- max(abs(pd))
    } else {
      for (.i in 1:.d) {
        for (.j in (1:.d)[-.i]) {
          entryrange <- max(entryrange, abs(pd[.i, .j]))
        }
      }
    }
  }
  
  
  rbPal <- colorRampPalette(c("red", "white", "blue"))
  color <- rbPal(breaks)
  plot(1, 1,
       col = "white",
       xlim = c(1, .d),
       ylim = c(1, .d),
       xlab = "Row Index",
       ylab = "Column Index",
       main = main)
  
  for (.i in 1:.d) {
    if (diagonal) {
      for (.j in 1:.d) {
        colorindex <- floor(
          breaks * (pd[.i, .j] + entryrange) / (2 * entryrange) + 0.99)
        points(x = c(.j), y = c(.d - .i + 1), col = color[colorindex],
               pch = 19, cex = cex)
      }
    } else {
      for (.j in (1:.d)[-.i]) {
        colorindex <- floor(
          breaks * (pd[.i, .j] + entryrange) / (2 * entryrange) + 0.99)
        points(x = c(.j), y = c(.d - .i + 1), col = color[colorindex],
               pch = 19, cex = cex)
      }
    }
  }
}