library(deSolve)
library(ggplot2)
library(tidyr)

###############
#Emtricitabine#
###############

# hiv_func is the function on the right hand side of the
# differential equation. 
rti_func <- function(t, HIV, params) {

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
    #
    N_RTI = params["N_RTI"]
    #
    beta_z = params["beta_z"]
    beta_v = params["beta_v"]
    p_y = params["p_y"]
    p_v = params["p_v"]
    lambda_x = params["lambda_x"]
    lambda_z= params["lambda_z"]


    # now that we have extracted the parameters, we can evaluate the
    # differential equation's right hand side

    dT = lambda_x - mu_x*T - N_RTI*beta_v*T*V
    dI = N_RTI*beta_v*T*V - mu_y*I - p_y*I*ZA
    dV = k_v*mu_y*I - mu_v*V - p_v*V*ZA
    dZ = lambda_z - mu_z*Z - beta_z*Z*V
    dZA = beta_z*Z*I -mu_z*ZA

    list(c(dT, dI, dV, dZ, dZA))
}

C_min <- .064
C_max <- 1.77
EC_50_min <- .0013
EC_50_max <- .64
Hill_Coefficient <- 1

N_RTI_min = (C_min ^ Hill_Coefficient) / ((C_min ^ Hill_Coefficient) + (EC_50_min ^ Hill_Coefficient))
N_RTI_max = (C_max ^ Hill_Coefficient) / ((C_max ^ Hill_Coefficient) + (EC_50_max ^ Hill_Coefficient))

#parameters
parms_min <- c(mu_x = 0.02,
           mu_y = 0.24,
           mu_v = 2.4,
           mu_z = 0.04,
           #
           N_RTI = N_RTI_min,
           #
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
           mu_v = 2.4,
           mu_z = 0.04,
           #
           N_RTI = N_RTI_max,
           #
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
times <- seq(from=0, to=365, by= 1)

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
gathered_result$variable[gathered_result$variable == "T"] <- "CD4+ T-Cells"
gathered_result$variable[gathered_result$variable == "I"] <- "HIV Infected T-Cells"
gathered_result$variable[gathered_result$variable == "V"] <- "Viral Load"
gathered_result$variable[gathered_result$variable == "Z"] <- "CD8+ T-Cells"
gathered_result$variable[gathered_result$variable == "ZA"] <- "Activated CD8+ T-Cells"

gathered_result$variable <- factor(gathered_result$variable, levels = c("CD4+ T-Cells", "HIV Infected T-Cells", "Viral Load",
                                                                        "CD8+ T-Cells", "Activated CD8+ T-Cells"))

gathered_result_wider <- gathered_result %>%
    pivot_wider(
        id_cols = c("time", "variable"),
        names_from = type,
        values_from = value     
        )



#not plotted
plot_result <- ggplot(data = gathered_result_wider,#subset(gathered_result, variable == "TR"), #plot declared equation(s) only
                      mapping=aes(x=time, ymin = max, ymax = min, color = variable, fill = variable)) +
                      #geom_line(linewidth = 1)+
                      geom_ribbon(alpha = .7) +
                      theme_classic() +

                      #plot each equation on its own graph
                      facet_wrap(~variable, scales = "free_y")+
                      coord_cartesian(
                      xlim = NULL,
                      ylim = c(0, NA))

plot_result_test <- ggplot(data = gathered_result, #subset(gathered_result, variable == "CD4+ T-Cells"), #plot declared equation(s) only
                      mapping=aes(x=time, y = value, color = type)) +
    geom_line(linewidth = 1)+
    scale_color_manual("Drug Efficacy", values = c(min = "red", max = "blue"))+
    ylab("Concentration (mm3)")+
    xlab("Time (days)")+
    
                      #geom_ribbon(alpha = .7) +
    theme_classic() +
    theme(legend.position = "bottom") +

                      #plot each equation on its own graph
                      facet_wrap(~variable, scales = "free_y")+#, ncol = 1)+
                      coord_cartesian(
                      xlim = NULL,
                      ylim = c(0, NA))

plot_result_cd4 <- ggplot(data = subset(gathered_result, variable == "CD4+ T-Cells"), #plot declared equation(s) only
                      mapping=aes(x=time, y = value, color = type)) +
    geom_line(linewidth = 1)+
    scale_color_manual("Drug Efficacy", values = c(min = "red", max = "blue"))+
    ylab("Concentration (mm3)")+
    xlab("Time (days)")+
    
                      #geom_ribbon(alpha = .7) +
                      theme_classic() +

                      #plot each equation on its own graph
                      #facet_wrap(~variable, scales = "free_y")+
                      coord_cartesian(
                      xlim = NULL,
                      ylim = c(0, NA))

plot_result_i <- ggplot(data = subset(gathered_result, variable == "Infected T-Cells"), #plot declared equation(s) only
                      mapping=aes(x=time, y = value, color = type)) +
    geom_line(linewidth = 1)+
    scale_color_manual("Drug Efficacy", values = c(min = "red", max = "blue"))+
    ylab("Concentration (mm3)")+
    xlab("Time (days)")+
    
                      #geom_ribbon(alpha = .7) +
                      theme_classic() +

                      #plot each equation on its own graph
                      facet_wrap(~variable, scales = "free_y")+
                      coord_cartesian(
                      xlim = NULL,
                      ylim = c(0, NA))

plot_result_v <- ggplot(data = subset(gathered_result, variable == "Virions"), #plot declared equation(s) only
                      mapping=aes(x=time, y = value, color = type)) +
    geom_line(linewidth = 1)+
    scale_color_manual("Drug Efficacy", values = c(min = "red", max = "blue"))+
    ylab("Concentration (mm3)")+
    xlab("Time (days)")+
    
                      #geom_ribbon(alpha = .7) +
                      theme_classic() +

                      #plot each equation on its own graph
                      facet_wrap(~variable, scales = "free_y")+
                      coord_cartesian(
                      xlim = NULL,
                      ylim = c(0, NA))

plot_result_z <- ggplot(data = subset(gathered_result, variable == "CD8+ T-Cells"), #plot declared equation(s) only
                      mapping=aes(x=time, y = value, color = type)) +
    geom_line(linewidth = 1)+
    scale_color_manual("Drug Efficacy", values = c(min = "red", max = "blue"))+
    ylab("Concentration (mm3)")+
    xlab("Time (days)")+
    
                      #geom_ribbon(alpha = .7) +
                      theme_classic() +

                      #plot each equation on its own graph
                      facet_wrap(~variable, scales = "free_y")+
                      coord_cartesian(
                      xlim = NULL,
                      ylim = c(0, NA))

plot_result_za <- ggplot(data = subset(gathered_result, variable == "Activated CD8+ T-Cells"), #plot declared equation(s) only
                      mapping=aes(x=time, y = value, color = type)) +
    geom_line(linewidth = 1)+
    scale_color_manual("Drug Efficacy", values = c(min = "red", max = "blue"))+
    ylab("Concentration (mm3)")+
    xlab("Time (days)")+
    
                      #geom_ribbon(alpha = .7) +
                      theme_classic() +

                      #plot each equation on its own graph
                      facet_wrap(~variable, scales = "free_y")+
                      coord_cartesian(
                      xlim = NULL,
                      ylim = c(0, NA))



		

#save to file
fname_base <- "result_ocm_rti"
extensions = c("png", "pdf")

fname_base_cd4 <- "result_ocm_cd4"
extensions = c("png", "pdf")

fname_base_i <- "result_ocm_i"
extensions = c("png", "pdf")

fname_base_v <- "result_ocm_v"
extensions = c("png", "pdf")

fname_base_z <- "result_ocm_z"
extensions = c("png", "pdf")

fname_base_za <- "result_ocm_za"
extensions = c("png", "pdf")

for (ext in extensions) {
    fname <- paste(fname_base, ".", ext, sep="")
    ggsave(fname, plot_result_test, height = 5, width = 9) #height,width of png/ pdf
    print(paste("wrote output to file ", fname))

        fname <- paste(fname_base_cd4, ".", ext, sep="")
    ggsave(fname, plot_result_cd4, height = 5, width = 5) #height,width of png/ pdf
    print(paste("wrote output to file ", fname))

        fname <- paste(fname_base_i, ".", ext, sep="")
    ggsave(fname, plot_result_i, height = 5, width = 5) #height,width of png/ pdf
    print(paste("wrote output to file ", fname))

        fname <- paste(fname_base_v, ".", ext, sep="")
    ggsave(fname, plot_result_v, height = 5, width = 5) #height,width of png/ pdf
    print(paste("wrote output to file ", fname))

        fname <- paste(fname_base_z, ".", ext, sep="")
    ggsave(fname, plot_result_z, height = 5, width = 5) #height,width of png/ pdf
    print(paste("wrote output to file ", fname))

        fname <- paste(fname_base_za, ".", ext, sep="")
    ggsave(fname, plot_result_za, height = 5, width = 5) #height,width of png/ pdf
    print(paste("wrote output to file ", fname))
    }
