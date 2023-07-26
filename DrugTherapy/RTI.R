library(deSolve)
library(ggplot2)
library(tidyr)
library(dplyr)


###############
#Emtricitabine#
###############


# hiv_func is the function on the right hand side of the
# differential equation. 
rti_func <- function(t, RTI, params) {

    TW = RTI["TW"]
    TR = RTI["TR"]
    I = RTI["I"]
    V = RTI["V"]
    Z = RTI["Z"]
    ZA = RTI["ZA"]

    lambda_T_W = params["lambda_T_W"]
    lambda_T_R = params["lambda_T_R"]
    X_W= params["X_W"]
    X_R = params["X_R"]
    mu_T_W = params["mu_T_W"]
    mu_T_R = params["mu_T_R"]
    #
    N_RTI = params["N_RTI"]
    #
    mu_I = params["mu_I"]
    epsilon_V = params["epsilon_V"]
    mu_V = params["mu_V"]
    alpha = params["alpha"]
    c = params["c"]
    lambda_Z = params["lambda_Z"]
    mu_Z = params["mu_Z"]
    beta = params["beta"]
    mu_Z_A = params["mu_Z_A"]

    # now that we have extracted the parameters, we can evaluate the
    # differential equation's right hand side

    #
    dTW = lambda_T_W - mu_T_W*TW -  N_RTI*X_W*TW*V
    dTR = lambda_T_R - mu_T_R*TR -  N_RTI*X_R*TR*V
    #dTW = lambda_T_W - mu_T_W*TW - (1 - N_RTI)*X_W*TW*V
    #dTR = lambda_T_R - mu_T_R*TR - (1 - N_RTI)*X_R*TR*V
    #
    dI = V*(X_W*TW + X_R*TR) - mu_I*I - alpha*I*ZA
    dV = epsilon_V*mu_I*I - mu_V*V
    dZ = lambda_Z + c*Z*I - mu_Z*Z - beta*Z*I
    dZA = beta*Z*I - mu_Z_A*ZA

    list(c(dTW, dTR, dI, dV, dZ, dZA))
}

C_min <- .064
C_max <- 1.77
EC_50_min <- .0013
EC_50_max <- .64
Hill_Coefficient <- 1

N_RTI_min = (C_min ^ Hill_Coefficient) / ((C_min ^ Hill_Coefficient) + (EC_50_min ^ Hill_Coefficient))
N_RTI_max = (C_max ^ Hill_Coefficient) / ((C_max ^ Hill_Coefficient) + (EC_50_max ^ Hill_Coefficient))


#parameters
parms_min <- c(lambda_T_W = 10,
      	   lambda_T_R = 0.03198,
      	   X_W= 0.000024,
      	   X_R = 0.01,
      	   mu_T_W = 0.01,
      	   mu_T_R = 0.01,
           #
           N_RTI = N_RTI_min,
           #
      	   mu_I = 0.5,
      	   epsilon_V = 100,
      	   mu_V = 3,
      	   alpha = 0.02,
      	   c = 0.000005,
      	   lambda_Z = 20,
      	   mu_Z = 0.06,
      	   beta = 0.004,
      	   mu_Z_A = 0.004)

parms_max <- c(lambda_T_W = 10,
      	   lambda_T_R = 0.03198,
      	   X_W= 0.000024,
      	   X_R = 0.01,
      	   mu_T_W = 0.01,
      	   mu_T_R = 0.01,
           #
           N_RTI = N_RTI_max,
           #
      	   mu_I = 0.5,
      	   epsilon_V = 100,
      	   mu_V = 3,
      	   alpha = 0.02,
      	   c = 0.000005,
      	   lambda_Z = 20,
      	   mu_Z = 0.06,
      	   beta = 0.004,
      	   mu_Z_A = 0.004)
	   

# initial values
Pstart <- c(TW = 1000,
       	    TR = 10,
       	    I = 100,
            V = 100,
	    Z = 500,
	    ZA = 30)

# set the times (full collection of x axis points)
times <- seq(from=0, to=1000, by= 15)


diffeq_result_min <- ode(
    func=rti_func,
    y=Pstart,
    times=times,
    parms=parms_min)

#max
diffeq_result_max <- ode(
    func=rti_func,
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


gathered_result_wider <- gathered_result %>%
    pivot_wider(
        id_cols = c("time", "variable"),
        names_from = type,
        values_from = value     
        )


plot_result <- ggplot(data = gathered_result_wider,#subset(gathered_result, variable == "TR"), #plot declared equation(s) only
                      mapping=aes(x=time, ymin = min, ymax = max, color = variable, fill = variable)) +
                      geom_ribbon(alpha = .7) + theme_classic() +

                      #plot each equation on its own graph
                      facet_wrap(~variable, scales = "free_y")


		

#save to file
fname_base <- "result_rti"
extensions = c("png", "pdf")
for (ext in extensions) {
    fname <- paste(fname_base, ".", ext, sep="")
    ggsave(fname, plot_result, height = 5, width = 9) #height,width of png/ pdf
    print(paste("wrote output to file ", fname))
    }


