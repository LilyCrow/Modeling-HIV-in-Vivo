library(deSolve)
library(ggplot2)
library(tidyr)

# hiv_func is the function on the right hand side of the
# differential equation. 
hiv_func <- function(t, HIV, params) {

    T = HIV["T"]
    I = HIV["I"]
    V = HIV["V"]
    Z = HIV["Z"]
    ZA = HIV["ZA"]

    mu_x = params["mu_x"]
    mu_y = params["mu_y"]
    mu_v = params["mu_v"]
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
    dV = k_v*mu_y*I - mu_v*V - p_v*V*ZA
    dZ = lambda_z - mu_z*Z - beta_z*Z*V
    dZA = beta_z*Z*I -mu_z*ZA

    list(c(dT, dI, dV, dZ, dZA))
}

#parameters
parms <- c(mu_x = 0.02,
           mu_y = 0.24,
           mu_v = 2.4,
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
	    Z = 500,
	    ZA = 0)

# set the times (full collection of x axis points)
times <- seq(from=0, to=365, by= 0.01)


diffeq_result <- ode(
    func=hiv_func,
    y=Pstart,
    times=times,
    parms=parms)


diffeq_result <- as.data.frame(diffeq_result)

gathered_result <- gather(diffeq_result, variable, value, -time)

gathered_result$variable[gathered_result$variable == "T"] <- "CD4+ T-Cells"
gathered_result$variable[gathered_result$variable == "I"] <- "HIV Infected T-Cells"
gathered_result$variable[gathered_result$variable == "V"] <- "Virions"
gathered_result$variable[gathered_result$variable == "Z"] <- "CD8+ T-Cells"
gathered_result$variable[gathered_result$variable == "ZA"] <- "Activated CD8+ T-Cells"

gathered_result$variable <- factor(gathered_result$variable, levels = c("CD4+ T-Cells", "HIV Infected T-Cells", "Virions",
                                                                        "CD8+ T-Cells", "Activated CD8+ T-Cells"))

plot_result <- ggplot(data = gathered_result,#subset(gathered_result, variable == "TR"), #plot declared equation(s) only
                      mapping=aes(x=time, y=value, color = variable)) +
    geom_line(linewidth = 0.5) +
    theme_classic() +
    #scale_color_manual("Type of Cell")+
    ylab("Concentration (mm3)")+
    xlab("Time(days)")+

                      #plot each equation on its own graph
    facet_wrap(~variable, scales = "free_y")+
    coord_cartesian(
        xlim = NULL,
        ylim = c(0, NA))


		

#save to file
fname_base <- "result_ocm"
extensions = c("png", "pdf")
for (ext in extensions) {
    fname <- paste(fname_base, ".", ext, sep="")
    ggsave(fname, plot_result, height = 5, width = 9) #height,width of png/ pdf
    print(paste("wrote output to file ", fname))
    }


