# Example of THIeF (vs ARIMA)
# Nikolaos Kourentzes, 2025
# nikolaos (at) kourentzes.com

# Load the MAPA package for the temporal
# aggregation function and the example data
library(MAPA)
# Load the tsutils package to get the S-matrix easily
library(tsutils)

# Get the sampling rate of the time series
f <- frequency(admissions)
# Find the permitted aggregation levels
k <- rev(f/(1:f))
# Allow only factors of f
k <- k[k %% 1 == 0]
# Find how many aggregation levels are to be used
p <- length(k)

# Create the temporally aggregate series
# I use the 58 first observations for the example.
Y <- tsaggr(head(MAPA::admissions,58),k,FALSE,FALSE)

# This is the desired forecast horizon
h <- 12
# Calculate the horizons for all termpoal aggregation levels
H <- 12/k

# Caonstruct the base forecast
# An independently produced forecast per level
frc <- list()
for (j in 1:p){
  
  # Get the series of a single level
  temp <- Y[[1]][[j]]
  # I will keep a test set, just for the plot that follows
  # This is not necessary, and you can just forecast right away
  n <- length(temp)
  tempTrn <- head(temp,n-H[j])
  # Fit an ARIMA for that level
  fit <- auto.arima(tempTrn)
  # Forecast and store the result
  frc[[j]] <- forecast(fit, h=H[j])$mean
  
} # Close loop of base forecasts

# Collect forecasts in a vector as needed for THieF
fbase <- cbind(unlist(rev(frc)))
# Get the S matrix
S <- tsutils::Sthief(Y[[1]]$AL1)
# Get the structural approximation W
W <- diag(1/rowSums(S))
# Estimate G
G <- solve(t(S) %*% W %*% S) %*% t(S) %*% W
# Reconcile base forecasts into THieF forecasts
freco <- S %*% G %*% fbase

# The last f values correspond to the forecasts
# on the original time series
# See also the plot below. 

# Get some colours for plotting.
cmp <- RColorBrewer::brewer.pal(7,"Set1")

# Create a plto fo the THieF and base forecasts
par(mfrow=c(3,2),mar=c(3.5,2,1.5,1))
for (j in 1:p){
  temp <- Y[[1]][[j]]  
  # Extract the appropriate level forecasts
  tempF <- cbind(head(tail(fbase,sum(H[1:j])),H[j]),
                 head(tail(freco,sum(H[1:j])),H[j]))
  tempF <- rbind(rbind(rep(head(tail(temp,H[j]+1),1),2)),tempF)
  # Do the plotting
  n <- length(temp)
  yy <- range(temp)
  yy[2] <- max(yy[2],max(tempF))
  yy <- yy + 0.04*c(-1,1)*diff(yy)
  plot(1:n,temp,xlim=c(1,n),xlab="",ylab="",type="l",ylim=yy,lwd=2)
  lines((n-H[j]):n,#n:(n+H[j]),
        tempF[,1],col=cmp[1])
  points((n-H[j]+1):n,#n:(n+H[j]),
        tempF[-1,1],col=cmp[1],pch=20)
  lines((n-H[j]):n,#n:(n+H[j]),
        tempF[,2],col=cmp[2])
  points((n-H[j]+1):n,#n:(n+H[j]),
         tempF[-1,2],col=cmp[2],pch=20)
  abline(v=n-H[j]+0.5,lty=2,col="black")
  
  title(xlab=c("Months","2-Months","Quarters","4-Months","Half-Years","Years")[j],line=2)
  text(1-(n)*0.05,yy[2]+diff(yy)*0.1,paste0("Aggregation level ",k[j]),pos=4,xpd=NA)
  
  legend("topleft",c("Data","ARIMA","THieF"),col=c("black",cmp[1:2]),lty=1,lwd=2,cex=0.9,bty="n")
}
