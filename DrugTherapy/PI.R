library(deSolve)
library(ggplot2)
library(tidyr)

############
#Amprenavir#
############

# hiv_func is the function on the right hand side of the
# differential equation. 
pi_func <- function(t, HIV, params) {

    T = HIV["T"]
    I = HIV["I"]
    V = HIV["V"]
    VNI = HIV["VNI"]
    Z = HIV["Z"]
    ZA = HIV["ZA"]

    mu_x = params["mu_x"]
    mu_y = params["mu_y"]
    mu_v = params["mu_v"]
    mu_vni = params["mu_vni"]
    N_PI = params["N_PI"]
    mu_z = params["mu_z"]
    k_v = params["k_v"]
    beta_z = params["beta_z"]
    beta_v = params["beta_v"]
    p_y = params["p_y"]
    p_v = params["p_v"]
    lambda_x = params["lambda_x"]
    lambda_z= params["lambda_z"]


    # now that we have extracted the parameters, we can evaluate the
    # differential equation's right hand side

    dT = lambda_x - mu_x*T - beta_v*T*V
    dI = beta_v*T*V - mu_y*I - p_y*I*ZA
    #
    dV = (1 - N_PI)*k_v*mu_y*I - mu_v*V - p_v*V*ZA
    dVNI = N_PI*k_v*mu_y*I - mu_v*V - p_v*V*ZA
    #
    dZ = lambda_z - mu_z*Z - beta_z*Z*V
    dZA = beta_z*Z*I -mu_z*ZA

    list(c(dT, dI, dV, dVNI, dZ, dZA))
}

C_min <- 0.28
C_max <- 4.63
EC_50_min <- 0.006
EC_50_max <- 0.04
Hill_Coefficient <- 3

N_PI_min = (C_min ^ Hill_Coefficient) / ((C_min ^ Hill_Coefficient) + (EC_50_min ^ Hill_Coefficient))
N_PI_max = (C_max ^ Hill_Coefficient) / ((C_max ^ Hill_Coefficient) + (EC_50_max ^ Hill_Coefficient))

#N_PI_min = .2
#N_PI_max = .8

#parameters
parms_min <- c(mu_x = 0.02,
           mu_y = 0.24,
           #
           mu_v = 2.4,
           mu_vni =2.4,
           N_PI = N_PI_min,
           #
           mu_z = 0.04,
           k_v = 360,
           beta_z = 0.000005,
           beta_v = 0.000024,
           p_y = 0.02,
           p_v = 0.02,
           lambda_x = 20,
           lambda_z = 20)

#parameters
parms_max <- c(mu_x = 0.02,
           mu_y = 0.24,
           #
           mu_v = 2.4,
           mu_vni =2.4,
           N_PI = N_PI_max,
           #
           mu_z = 0.04,
           k_v = 360,
           beta_z = 0.000005,
           beta_v = 0.000024,
           p_y = 0.02,
           p_v = 0.02,
           lambda_x = 20,
           lambda_z = 20)
	   

# initial values
Pstart <- c(T = 1000,
       	    I = 0,
            V = 0.001,
            VNI = 0,
	    Z = 500,
	    ZA = 0)

# set the times (full collection of x axis points)
times <- seq(from=0, to=365, by= 1)

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

gathered_result_wider <- gathered_result %>%
    pivot_wider(
        id_cols = c("time", "variable"),
        names_from = type,
        values_from = value     
        )



plot_result <- ggplot(data = gathered_result,#subset(gathered_result, variable == "TR"), #plot declared equation(s) only
                      mapping=aes(x=time, y= value, color = type))+ #, fill = variable)) +
                      #geom_ribbon(alpha = .7) + theme_classic() +
                      geom_line(linewidth = 1)+
                      #plot each equation on its own graph
                      facet_wrap(~variable, scales = "free_y")
		

#save to file
fname_base <- "result_ocm_pi"
extensions = c("png", "pdf")
for (ext in extensions) {
    fname <- paste(fname_base, ".", ext, sep="")
    ggsave(fname, plot_result, height = 5, width = 9) #height,width of png/ pdf
    print(paste("wrote output to file ", fname))
    }


