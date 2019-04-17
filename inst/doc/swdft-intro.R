## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>", 
  fig.width = 7, 
  fig.height = 5
)

## ---- include=FALSE------------------------------------------------------
if (as.character(Sys.info()[8]) == "lee") {
  devtools::load_all("/home/lee/Dropbox/swdft/r/swdft")
} else {
  library(swdft)
}

## ---- eval=FALSE---------------------------------------------------------
#  library(swdft)

## ------------------------------------------------------------------------
N <- 100
window_size <- 32
x_random <- rnorm(n=N, mean=0, sd=1)
a <- swdft::swdft(x=x_random, n=window_size)

## ------------------------------------------------------------------------
# --- Generate a local cosine signal + Gaussian noise ---
set.seed(999)

## Specify the parameters for the local cosine signal 
signal_length <- 96 # The length of the local cosine signal
freq <- 3 / 32 # Frequency, interpreted as "3 complete cycles every 32 data-points"
amplitude <- 1 # The amplitude of the local cosine signal 
phase <- 0 # Phase of the local cosine signal 
periodic_start <- 20 # When the local cosine signal starts 
periodic_length <- 50 # How long the local cosine signal lasts 
signal <- swdft::local_signal(N=signal_length, A=amplitude, 
                              Fr=freq, phase=phase, S=periodic_start, 
                              L=periodic_length)

## Generate the Gaussian noise factor, and add to the local cosine signal 
sigma <- .3 # Standard deviation of the Gaussian noise factor 
noise <- rnorm(n=signal_length, mean=0, sd=sigma)
x <- signal + noise 

## Plot what the local cosine signal + noise looks like 
plot(x, pch=19, cex=1.4, xlab="", ylab="Signal", 
     main="Local Cosine Signal plus Gaussian Noise")
lines(signal, lwd=2, col="red")

## ------------------------------------------------------------------------
window_size <- 32
a <- swdft::swdft(x=x, n=window_size)

## ------------------------------------------------------------------------
length(x)
dim(a$a)

## ------------------------------------------------------------------------
a_nopad <- swdft::swdft(x=x, n=window_size, pad=FALSE)
dim(a_nopad$a)

## ------------------------------------------------------------------------
amod <- Mod(a$a)^2

## ---- eval=FALSE---------------------------------------------------------
#  swdft_fft <- function(x, n, taper) {
#    N <- length(x)
#    P <- N - n + 1
#    a <- array(data = NA, dim = c(n, P))
#  
#    for (p in n:N) {
#      a[, p - n + 1] <- stats::fft(z = x[(p - n + 1):p] * taper)
#    }
#  
#    return(a)
#  }

## ------------------------------------------------------------------------
plot(a)

## ------------------------------------------------------------------------
plot(a, col="tim.colors", title="SWDFT using the 'tim.colors' colorscale")

## ------------------------------------------------------------------------
sampling_frequency <- 100
plot(a, freq_type="hertz", fs=sampling_frequency)

## ------------------------------------------------------------------------
plot(a, freq_type="fraction")

## ------------------------------------------------------------------------
years <- 1900:1995
plot(a, freq_type="fraction", col="tim.colors", custom_xaxis=years, xlab="Years")

## ------------------------------------------------------------------------
a_taper <- swdft::swdft(x=x, n=window_size, taper_type="cosine", p=.2)
plot(a_taper)

## ------------------------------------------------------------------------
a_smooth <- swdft::swdft(x=x, n=window_size, smooth="daniell", m=1, num_convs=0)
plot(a_smooth)

## ------------------------------------------------------------------------
cosreg_fit <- swdft::cosreg(x=x, f=freq)
plot(cosreg_fit)

## ------------------------------------------------------------------------
coefficients(cosreg_fit)

## ------------------------------------------------------------------------
periodogram <- Mod(fft(x))^2
freqs <- (0:(length(x) - 1)) / length(x)
max_freq <- freqs[which.max(periodogram)]
cat("Estimated Frequency: ", max_freq, " True Frequency: ", freq, " \n")

## ------------------------------------------------------------------------
swdft::get_max_freq(x=x) == max_freq

## ------------------------------------------------------------------------
local_cosreg_fit <- swdft::local_cosreg(x=x)
plot(local_cosreg_fit)
coefficients(local_cosreg_fit)

## ------------------------------------------------------------------------
set.seed(666)

## Set the frequency and length of the signal 
N <- 128 
f0 <- 10/N

## Generate a time-varying amplitude 
amplitude <- rep(0, N)
inds11 <- 10:20
inds12 <- 21:50
inds13 <- 51:70
amplitude[inds11] <- seq(0, 1, length=length(inds11))
amplitude[inds12] <- seq(1, 1, length=length(inds12))
amplitude[inds13] <- seq(1, 0, length=length(inds13))

## Generate the cosine signal with time-varying amplitude plus noise
signal <- swdft::cosine(N=N, A=amplitude, Fr=f0, phase=0)
noise <- rnorm(n=N, mean=0, sd=sigma)
x_demod <- signal + noise

## Plot what the signal looks like 
plot(x_demod, pch=19, cex=1.4, main="Cosine signal with time-varying amplitude plus Gaussian noise", 
     xlab="", ylab="")
lines(signal, lwd=2, col="red")

## ------------------------------------------------------------------------
complex_demod_fit <- swdft::complex_demod(x=x_demod, f0=f0)
plot(complex_demod_fit)
lines(signal, col="blue")
legend("topright", col=c("black", "red", "blue"), 
       c("Signal + Noise", "Complex Demodulation", "True Signal"), lwd=1)

## ------------------------------------------------------------------------
complex_demod_fit_smoother <- swdft::complex_demod(x=x_demod, f0=f0, passfreq=.05)
plot(complex_demod_fit_smoother)
lines(signal, col="blue")
legend("topright", col=c("black", "red", "blue"), 
       c("Signal + Noise", "Complex Demodulation", "True Signal"), lwd=1)

## ------------------------------------------------------------------------
## Generate the time-varying amplitude for the second periodic component 
amplitude2 <- rep(0, N)
inds21 <- 50:70
inds22 <- 71:100
inds23 <- 101:120
amplitude2[inds21] <- seq(0, 1, length=length(inds21))
amplitude2[inds22] <- seq(1, 1, length=length(inds22))
amplitude2[inds23] <- seq(1, 0, length=length(inds23))

## Set the frequency for the second periodic component
f1 <- 30 / N

## Generate the second signal and add it to the first 
signal2 <- swdft::cosine(N=N, A=amplitude2, Fr=f1, phase=phase)
x_demod <- x_demod + signal2

## Plot what the two signals plus noise look like
plot(x_demod, pch=19, cex=1, main="Two cosine signals with time-varying amplitudes plus noise", 
     xlab="")
lines(signal, col="red", lwd=2)
lines(signal2, col="blue", lwd=2)

## ------------------------------------------------------------------------
matching_demod_fit <- swdft::matching_demod(x=x_demod, n=70)
plot(matching_demod_fit)

## ------------------------------------------------------------------------
plot(x_demod, cex=1, pch=19, ylim=c(-2, 3), xlab="", ylab="", 
     main="Fitted values from each iteration of matching demodulation")
lines(signal, col="red")
lines(signal2, col="blue")

lines(matching_demod_fit$iterations$iter_fits[1, ], col="green", lwd=2, lty=2)
lines(matching_demod_fit$iterations$iter_fits[2, ], col="purple", lwd=2, lty=2)

legend("topright", c("Signal + Noise", "Signal 1", "Signal 2", "Iteration 1", "Iteration 2"), 
       col=c("black", "red", "blue", "green", "purple"), lwd=2, cex=.8)

## ------------------------------------------------------------------------
plot(matching_demod_fit$coefficients$inst_amp[1, ], type="l", lwd=2, col="red", ylim=c(0, 1.5), 
     main="Time-varying amplitude for matching demodulation iterations", xlab="", ylab="")
lines(matching_demod_fit$coefficients$inst_amp[2, ], lwd=2, col="blue")
legend("topright", col=c("red", "blue"), lwd=2, 
       paste0("Iteration ", 1:2, ", Frequency: ", round(matching_demod_fit$coefficients$f0, digits=3)))

