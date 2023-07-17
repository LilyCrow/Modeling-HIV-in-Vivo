require(deSolve)
require(ggplot2)
require(gridExtra)
require(dplyr)

# parameters
v_HIV_parm1 = c(s =6 , p = 2, tmax = 54, dt = 1/7, k =.3 , delta = .49, N = 10.3e9, c = 3.07)

# initial values
v_HIV_init_val <- c(T = 0.999999, I = 1e-6, V = 0)

t <- seq(0, 200, by = 5)

model.HIV <- function(t, x, v_HIV_parm1) {
    with(as.list(c(v_HIV_parm1, x)), {
        dT <- s + p*T-dt*T-k*V*T
        dI <- k*V*T - delta*I
        dV <- N*delta*I - c*V

        dx <- c(dT, dI, dV)
        list(dx)
    })
}

#solve system of equations
out <- lsoda(v_HIV_init_val, t, model.HIV, v_HIV_parm1)

#plot data
mf <- par("mfrow")
plot(out, main = c("T", "I", "V"))
plot(out[,"T"], out[,"I"],  type = "l", xlab = "time", ylab = "percent of cell population")
par(mfrow = mf)
