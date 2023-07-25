library(deSolve)
library(ggplot2)
library(tidyr)

############
#Amprenavir#
############

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

    dTW = lambda_T_W - mu_T_W*TW - X_W*TW*VI
    dTR = lambda_T_R - mu_T_R*TR - X_R*TR*VI
    dI = VI*(X_W*TW + X_R*TR) - mu_I*I - alpha*I*ZA
    #
    dVI = (1 - N_PI)*epsilon_V*mu_V_I*I - mu_V_I*VI
    dVNI = N_PI*epsilon_V*mu_V_NI*I - mu_V_NI*VNI
    #
    dZ = lambda_Z + c*Z*I - mu_Z*Z - beta*Z*I
    dZA = beta*Z*I - mu_Z_A*ZA

    list(c(dTW, dTR, dI, dVI, dVNI, dZ, dZA))
}

C_min <- .28
C_max <- 4.63
EC_50_min <- .012
EC_50_max <- .08
Hill_Coefficient <- 3

N_PI_min = (C_min ^ Hill_Coefficient) / ((C_min ^ Hill_Coefficient) + (EC_50_min ^ Hill_Coefficient))
N_PI_max = (C_max ^ Hill_Coefficient) / ((C_max ^ Hill_Coefficient) + (EC_50_max ^ Hill_Coefficient))

#parameters
parms_min <- c(lambda_T_W = 10,
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
      	   mu_Z_A = 0.004,
           #
           N_PI = N_PI_min,
           epsilon_V = .3,
           mu_V_I = 4,
           mu_V_NI = 4)

parms_max <- c(lambda_T_W = 10,
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
      	   mu_Z_A = 0.004,
           #
           N_PI = N_PI_max,
           epsilon_V = 100,
           mu_V_I = 3,
           mu_V_NI = 3)

	   

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




diffeq_result_min <- ode(
    func=pi_func,
    y=Pstart,
    times=times,
    parms=parms_min)

diffeq_result_max <- ode(
    func=pi_func,
    y=Pstart,
    times=times,
    parms=parms_max)


diffeq_result_min <- as.data.frame(diffeq_result_min)
diffeq_result_max <- as.data.frame(diffeq_result_max)

gathered_result_min <- gather(diffeq_result_min, variable, value, -time)
gathered_result_min$type <- "min"

gathered_result_max <- gather(diffeq_result_max, variable, value, -time)
gathered_result_max$type <- "max"

gathered_result <- rbind(gathered_result_min, gathered_result_max)

head(gathered_result)

plot_result <- ggplot(data = gathered_result_min,#subset(gathered_result, variable == "TR"), #plot declared equation(s) only
                      mapping=aes(x=time, y=value, color = variable)) +
    geom_line(linewidth=1) +
    geom_line(data = gathered_result_max, linewidth = 1) +
    geom_ribbon(aes(ymin = N_RTI_min, ymax = N_RTI_max), fill = "grey70") +
    theme_classic() +

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
