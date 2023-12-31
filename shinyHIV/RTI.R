library(deSolve)
library(ggplot2)
library(tidyr)

###############
#Emtricitabine#
###############

# rti_func is the function on the right hand side of the
# differential equation. 
rti_func <- function(t, RTI, params) {

    T = RTI["T"]
    I = RTI["I"]
    V = RTI["V"]
    Z = RTI["Z"]
    ZA = RTI["ZA"]

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

#Minimum efficacy parameters where N_RTI = N_RTI_min
parms_min <- c(mu_x = 0.02,    #Natural death rate of uninfected CD4+ T-cells
           mu_y = 0.24,        #Natural death rate of infected CD4+ T-cells
           mu_v = 2.4,         #Natural death rate of viral particles
           mu_z = 0.04,        #Natural death rate of  CD8+ T-cells
           #
           N_RTI = N_RTI_min,  #Emtricitabine efficacy is equal to the result of N_RTI_min
           #
           k_v = 360,          #Average number of virions released per infected cell
           beta_z = 0.000005,  #CD8+ T-cell activation rate
           beta_v = 0.000024,  #Rate at which CD4+ T-cells become infected
           p_y = 0.02,         #Rate at which infected CD4+ T-cells are destroyed
           p_v = 0.02,         #Rate at which viruses and virions are destroyed
           lambda_x = 20,      #Rate at which CD4+ T-cells are created by the body
           lambda_z = 20       #Rate at which CD8+ T-cells are created by the body
           )
	   

#Maximum efficacy parameters where N_RTI = N_RTI_max
parms_max <- c(mu_x = 0.02,    #Natural death rate of uninfected CD4+ T-cells
           mu_y = 0.24,        #Natural death rate of infected CD4+ T-cells
           mu_v = 2.4,         #Natural death rate of viral particles
           mu_z = 0.04,        #Natural death rate of  CD8+ T-cells
           #
           N_RTI = N_RTI_max,  #Emtricitabine efficacy is equal to the result of N_RTI_max
           #
           k_v = 360,          #Average number of virions released per infected cell
           beta_z = 0.000005,  #CD8+ T-cell activation rate
           beta_v = 0.000024,  #Rate at which CD4+ T-cells become infected
           p_y = 0.02,         #Rate at which infected CD4+ T-cells are destroyed
           p_v = 0.02,         #Rate at which viruses and virions are destroyed
           lambda_x = 20,      #Rate at which CD4+ T-cells are created by the body
           lambda_z = 20       #Rate at which CD8+ T-cells are created by the body
           )


#Initial values
Pstart <- c(T = 1000,  #CD4+ T-cells
       	    I = 0,     #Infected CD4+ T-cells
            V = 0.001, #Viral particles
	    Z = 500,   #CD8+ T-cells
	    ZA = 0)    #Activated CD8+ T-cells

#Time steps to solve equation for (full collection of x-axis points)
times <- seq(from=0, to=365, by= 1)


#Now inoke the ODE function.
#It requires four things to run: the right hand side of the function
#(rti_func), the initial values (Pstart), time steps that need to be
#solved for (times), and equation constants (parms_min or parms_max)

#Min
diffeq_result_min <- ode(
    func=rti_func,
    y=Pstart,
    times=times,
    parms=parms_min)

#Max
diffeq_result_max <- ode(
    func=rti_func,
    y=Pstart,
    times=times,
    parms=parms_max)

#Take the solved equations (diffeq_result_min and diffeq_result_max) and put them into a data frame.
#ggplot2 requires such large data sets to be in a data frame for its
#plotting function
diffeq_result_min <- as.data.frame(diffeq_result_min)
diffeq_result_max <- as.data.frame(diffeq_result_max)

#The R gather() function, extracts the data needed by ggplot.
#It also helps to massage the data in a way that preps it to be plotted. 
gathered_result_min <- gather(diffeq_result_min, variable, value, -time)

#Create a separate column of type "min" for diffeq_result_min values
#This will help ggplot2 tell the difference between the outputs for the same
#time (x value), and enable it to make two different plots
gathered_result_min$type <- "min"

#Perform the same operations for max
gathered_result_max <- gather(diffeq_result_max, variable, value, -time)
gathered_result_max$type <- "max"

#Puts minimum and maximum efficacy values into one variable to be plotted
gathered_result_rti <- rbind(gathered_result_min, gathered_result_max)

#This provides each variable, and therefore each graph, with a descriptive title.
gathered_result_rti$variable[gathered_result_rti$variable == "T"] <- "CD4+ T-Cells"
gathered_result_rti$variable[gathered_result_rti$variable == "I"] <- "Infected CD4+ T-Cells"
gathered_result_rti$variable[gathered_result_rti$variable == "V"] <- "Viral Load"
gathered_result_rti$variable[gathered_result_rti$variable == "Z"] <- "CD8+ T-Cells"
gathered_result_rti$variable[gathered_result_rti$variable == "ZA"] <- "Activated CD8+ T-Cells"
#The default order of ggplot2 is alphabetic, the factor() function enables 
#you to change the order of the plots into something more logical
#or relevant to your purposes.
gathered_result_rti$variable <- factor(gathered_result_rti$variable, levels = c("CD4+ T-Cells", "Infected CD4+ T-Cells", "Viral Load",
                                                                        "CD8+ T-Cells", "Activated CD8+ T-Cells"))

gathered_result_wider <- gathered_result_rti %>%
    pivot_wider(
        id_cols = c("time", "variable"),
        names_from = type,
        values_from = value     
        )


#choicesVar rerpresents types of cells to be plotted as selected by the user
#choicesDT will represent which ART is selected by the user
#TODO: assign an appropriate method to establish which model to use (i.e. HIV, RTI, TODO: PI, II)
shinyRTI <- function(choicesVar){###, choicesDT){

#if(choicesDT == "Reverse Transcriptase Inhibitor (RTI)"){
#plot the data and save it in a new variable: plot_result.
plot_result <- ggplot(data = subset(gathered_result_rti, variable %in% choicesVar),
                      mapping=aes(x=time, y = value, color = type)) +
                      geom_line(linewidth = 1) + #Linewidth
                      theme_classic() + #ggplot has several different themes
                      #that can be added to a plot here we are using "classic"
                      scale_color_manual("Drug Efficacy", values = c(min = "red", max = "blue")) + #Make legend
                      theme(legend.position = "bottom") + #Locate legend at the bottom of the image
                      ylab("Concentration (mm3)") + #Label y-axis "Concentration (mm3)"
                      xlab("Time (days)") +  #Label x-axis "Time (days)"
                      facet_wrap(~variable, scales = "free_y") +  #plot each equation on its own graph
                      coord_cartesian(
                      xlim = NULL,
                      ylim = c(0, NA)) #Y-axis begins at 0
    return(plot_result)#}
}

