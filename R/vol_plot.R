#' @title stovol
#'
#' @description compare the absolute values of the returns with the estimated volatilizes
#'
#' @param y the training time series
#'
#' @param z The estimated volatility
#'
#' @examples
#' vol_plot(z)
#' @export

vol_plot = function(y, z){

  m = length(z)


    op=par(mfrow = c(2,1),
           #omg=c(2,2,0,0)+0.05,
           mar=c(2.5,2.5,3.0,3.0) +0.05)

    nn=1:m
    plot(nn, abs(y), type="l", lwd=2, col="blue",xlab="", ylab="",
         main=expression(paste("Time series of the absolute returns")), cex.main=1.5,
         cex.lab=1.9, cex.axis = 1.5)


    plot(nn, exp(z/2), type="l", lwd=2, col="red",
         main=expression( "Time series of the estimated volatilities"), cex.main=1.5,
         cex.lab=1.9, cex.axis = 1.5)

}




