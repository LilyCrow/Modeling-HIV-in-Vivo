# more learning to do differential equations in R.  You can run this
# by pasting it into an interactive R interpreter, or you can run
## R -f ode_1st_coupled.R

# I have put comments in detailing what's needed for R's deSolve() and
# then to make a simple plot.  This is a "sequel" to the
# od1_1st_single.R file.

# remember to do install.packages("ggplot2") (and same for deSolve) if
# needed.
library(deSolve)
library(ggplot2)
library(tidyr)

# tiv_func is the function on the right hand side of the
# differential equation. 
tiv_func <- function(t, TIV, params) {
    # this takes several lines to define, but really it's just the
    # two right hand side expressions in Lotka Volterra
    T = TIV["T"]
    I = TIV["I"]
    V = TIV["V"]
    # using the notation from the Wikipedia version of the equations
                                        # we have alpha, beta, gamma, and delta in the params structure
    s = params["s"]
    p = params ["p"]
    tmax = params["tmax"]
    d = params["d"]
    k = params["k"]
    delta = params["delta"]
    N = params["N"]
    c = params["c"]

    # now that we have extracted the parameters, we can evaluate the
    # differential equation's right hand side

    dT_dt = s + p*T - d*T - k*V*T
    dI_dt = k*V*T - delta*I
    dV_dt = N*delta*I - c*V
    list(c(dT_dt, dI_dt, dV_dt))
}

# set the constant parameters of the equation.  for this one we have
# alpha (max prey per capita growth rate), beta (effect of predators
# on prey growth rate), gamma (predator per capita death rate), and
# delta (effect of presence of prey on predator growth rate) as used
# in the wikipedia page.  
parms <- c(s = 34,
           p = .7 ,
           tmax = .8,
           d = 4,
           k = .222,
           delta = 0.49,
           N = 10e10,
           c = 3.07)
# set the times (full collection of x axis points)
times <- seq(from=0, to=200, by=1/100)
# initial values for both predator and prey
Pstart <- c(T = 0.99999999,
            I = 11,
            V = 0)

# invoke ode - it requires the four things we just defined: the right
# hand side function, the initial value, the times of operation, and
# the constant parameters
diffeq_result <- ode(
    func=tiv_func,
    y=Pstart,
    times=times,
    parms=parms
)

# R now likes to force the result into a data frame, which is a
# pervasive way of organizing data in R.  this as.data.frame()
# function enforces that diffeq_result will be usable as input to
# ggplot2's plotting functions.
diffeq_result <- as.data.frame(diffeq_result)

# here we use a simple invocation of ggplot2 to make the plot.  if
# this is invoked interactively then we will get a pop-up plot as well
# as saving it to png and pdf files; if it's invoked in batch mode
# then it just saves the output files.
##
# A particular note here: we have to extract the two columns of the
# data frame.  People who use R's tidyverse have a way of doing that
# with the gather() function, which I try to imitate here.

# gather (from tidyr) seems to be used widely - it seems to do the
# massaging of a data frame to make it ready for ggplot()
gathered_result <- gather(diffeq_result, variable, value, -time)
plot_result <- ggplot(data=gathered_result,
                      mapping=aes(x=time, y=value, color=variable)) +
    geom_line(linewidth=2) + theme_classic()

# now save to file
fname_base <- "result_1_tiv"
extensions = c("png", "pdf")
for (ext in extensions) {
    fname <- paste(fname_base, ".", ext, sep="")
    ggsave(fname, plot_result)
    print(paste("wrote output to file ", fname))
}
