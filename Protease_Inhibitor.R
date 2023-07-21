library(deSolve)
library(ggplot2)
library(tidyr)

# hiv_func is the function on the right hand side of the
# differential equation. 
pi_func <- function(t, PI, params) {

    TW = PI["TW"]
    TR = PI["TR"]
    I = PI["I"]
    VI = PI["VI"]
    VNI = PI["VNI"]
    Z = PI["Z"]
    ZA = PI["ZA"]

    lambda_T_W = params["lambda_T_W"]
    lambda_T_R = params["lambda_T_R"]
    X_W= params["X_W"]
    X_R = params["X_R"]
    mu_T_W = params["mu_T_W"]
    mu_T_R = params["mu_T_R"]
    mu_I = params["mu_I"]
#
    N_PI = params["N_PI"]
    epsilon_V = params["epsilon_V"]
    mu_V_I = params["mu_V_I"]
    mu_V_NI = params["mu_V_NI"]
#
    alpha = params["alpha"]
    c = params["c"]
    lambda_Z = params["lambda_Z"]
    mu_Z = params["mu_Z"]
    beta = params["beta"]
    mu_Z_A = params["mu_Z_A"]

    # now that we have extracted the parameters, we can evaluate the
    # differential equation's right hand side

    dTW = lambda_T_W - mu_T_W*TW - X_W*TW*V
    dTR = lambda_T_R - mu_T_R*TR - X_R*TR*V
    dI = V*(X_W*TW + X_R*TR) - mu_I*I - alpha*I*ZA
    #
    dVI = (1 - N_PI)*epsilon_V*mu_V_I*I - mu_V_I*VI
    dVNI = N_PI*epsilon_V*mu_V_NI*I - mu_V_NI*VNI
    #
    dZ = lambda_Z + c*Z*I - mu_Z*Z - beta*Z*I
    dZA = beta*Z*I - mu_Z_A*ZA

    list(c(dTW, dTR, dI, dVI, dVNI, dZ, dZA))
}

#parameters
parms <- c(lambda_T_W = 10,
      	   lambda_T_R = 0.03198,
      	   X_W= 0.000024,
      	   X_R = 0.01,
      	   mu_T_W = 0.01,
      	   mu_T_R = 0.01,
      	   mu_I = 0.5,
      	   alpha = 0.02,
      	   c = 0.000005,
      	   lambda_Z = 20,
      	   mu_Z = 0.06,
      	   beta = 0.004,
      	   mu_Z_A = 0.004
           #
           N_PI = ,
           epsilon_V = ,
           mu_V_I = ,
           mu_V_NI = )

	   

# initial values
Pstart <- c(TW = 1000,
       	    TR = 10,
       	    I = 10,
            #
            VI = 10,
            VNI = 0,
            #
            Z = 500,
	    ZA = 30)

# set the times (full collection of x axis points)
times <- seq(from=0, to=2000, by= 25)


diffeq_result <- ode(
    func=pi_func,
    y=Pstart,
    times=times,
    parms=parms)


diffeq_result <- as.data.frame(diffeq_result)

gathered_result <- gather(diffeq_result, variable, value, -time)

plot_result <- ggplot(data = gathered_result,#subset(gathered_result, variable == "TR"), #plot declared equation(s) only
                      mapping=aes(x=time, y=value, color = variable)) +
                      geom_line(linewidth=1) + theme_classic() #+

                      #plot each equation on its own graph
                      facet_wrap(~variable, scales = "free_y")


		

#save to file
fname_base <- "result_pi"
extensions = c("png", "pdf")
for (ext in extensions) {
    fname <- paste(fname_base, ".", ext, sep="")
    ggsave(fname, plot_result, height = 5, width = 9) #height,width of png/ pdf
    print(paste("wrote output to file ", fname))
    }
