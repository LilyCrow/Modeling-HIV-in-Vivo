library(deSolve)
library(ggplot2)
library(tidyr)

############
#Saquinavir#
############

# hiv_func puts the equations, initial values, and parameters
#into one variable so that the ODE can solve it.
pi_func <- function(t, PI, params) {

    #the following two blocks of code explain to R what
    #what the variables inside the function are (i.e. makes
    #the connection between the in func. variables and the
    #external variables that hold their values.

    #Dependent variables (what you are solving for)
    T = PI["T"]    #CD4+ T-cells
    I = PI["I"]    #Infected CD4+ T-cells
    VI = PI["VI"]    #Viral load (number of viruses/ virions in blood)
    VNI = PI["VNI"]
    Z = PI["Z"]    #CD8+ T-cells
    ZA = PI["ZA"]  #Activated CD8+ T-cells


    #Independent variables (physical constants, etc.)
    mu_T = params["mu_T"]          #Natural death rate of uninfected CD4+ T-cells
    mu_I = params["mu_I"]          #Natural death rate of infected CD4+ T-cells
    mu_v = params["mu_v"]          #Natural death rate of viruses and virions
    mu_z = params["mu_z"]          #Natural death rate of CD8+ T-cells
    k_v = params["k_v"]            #Average number of virions released per infected cell
    N_PI = params["N_PI"]
    beta_z = params["beta_z"]      #CD8+ T-cell activation rate
    beta_v = params["beta_v"]      #Rate at which CD4+ T-cells become infected
    p_I = params["p_I"]            #Rate at which infected CD4+ T-cells are destroyed
    p_v = params["p_v"]            #Rate at which viruses and virions are destroyed
    lambda_T = params["lambda_T"]  #Rate at which CD4+ T-cells are created by the body
    lambda_z = params["lambda_z"]  #Rate at which CD8+ T-cells are created by the body


    # now that we have extracted the parameters, we can evaluate the
    # differential equation's right hand side
    dT = lambda_T - mu_T*T - beta_v*T*VI
    dI = beta_v*T*VI - mu_I*I - p_I*I*ZA
    dVI = (1 - N_PI)*k_v*mu_I*I - mu_v*VI - p_v*VI*ZA
    dVNI = N_PI*k_v*mu_I*I - mu_v*VNI - p_v*VNI*ZA
    dZ = lambda_z - mu_z*Z - beta_z*Z*VI
    dZA = beta_z*Z*I -mu_z*ZA

    #Put all 5 equations into a list so the ODE can read it  
    list(c(dT, dI, dVI, dVNI, dZ, dZA))
}

C_min <- 4.78
C_max <- 9.09
EC_50_min <- 40.9
EC_50_max <- 49.7
Hill_Coefficient <- 3.68

N_PI_min = log(1 + ((C_min/EC_50_min) ^ Hill_Coefficient))
N_PI_max = log(1 + ((C_max/EC_50_max) ^ Hill_Coefficient))
#N_PI_min = (C_min ^ Hill_Coefficient) / ((C_min ^ Hill_Coefficient) + (EC_50_min ^ Hill_Coefficient))
#N_PI_max = (C_max ^ Hill_Coefficient) / ((C_max ^ Hill_Coefficient) + (EC_50_max ^ Hill_Coefficient))

#Parameters
parms_min <- c(mu_T = 0.02,        #Natural death rate of uninfected CD4+ T-cells
           mu_I = 0.24,        #Natural death rate of infected CD4+ T-cells
           mu_v = 2.4,         #Natural death rate of viruses and virions
           mu_z = 0.04,        #Natural death rate of CD8+ T-cells
           k_v = 360,          #Average number of virions released per infected cell
           N_PI = N_PI_min,
           beta_z = 0.000005,  #CD8+ T-cell activation rate
           beta_v = 0.000024,  #Rate at which CD4+ T-cells become infected
           p_I = 0.02,         #Rate at which infected CD4+ T-cells are destroyed
           p_v = 0.02,         #Rate at which viruses and virions are destroyed
           lambda_T = 20,      #Rate at which CD4+ T-cells are created by the body
           lambda_z = 20       #Rate at which CD8+ T-cells are created by the body
           )

#Parameters
parms_max <- c(mu_T = 0.02,        #Natural death rate of uninfected CD4+ T-cells
           mu_I = 0.24,        #Natural death rate of infected CD4+ T-cells
           mu_v = 2.4,         #Natural death rate of viruses and virions
           mu_z = 0.04,        #Natural death rate of CD8+ T-cells
           k_v = 360,          #Average number of virions released per infected cell
           N_PI = N_PI_max,
           beta_z = 0.000005,  #CD8+ T-cell activation rate
           beta_v = 0.000024,  #Rate at which CD4+ T-cells become infected
           p_I = 0.02,         #Rate at which infected CD4+ T-cells are destroyed
           p_v = 0.02,         #Rate at which viruses and virions are destroyed
           lambda_T = 20,      #Rate at which CD4+ T-cells are created by the body
           lambda_z = 20       #Rate at which CD8+ T-cells are created by the body
           )
	   

#Initial values
Pstart <- c(T = 1000,   #CD4+ T-cells
       	    I = 0,      #Infected CD4+ T-cells
            VI = 0.001,  #Viral load 
            VNI = 0,
            Z = 500,    #CD8+ T-cells
	    ZA = 0)     #Activated CD8+ T-cells

#Time steps to solve equation for (full collection of x-axis points)
times <- seq(from=0, to=365, by= 0.01)



#Now inoke the ODE function.
#It requires four things to run: the right hand side of the function
#(hiv_func), the initial values (Pstart), time steps that need to be
#solved for (times), and equation constants (parms)
#Min
diffeq_result_min <- ode(
    func=pi_func,
    y=Pstart,
    times=times,
    parms=parms_min)

#Max
diffeq_result_max <- ode(
    func=pi_func,
    y=Pstart,
    times=times,
    parms=parms_max)

#Take the solved equations (diffeq_result) and put it into a data frame
#ggplot2 requires such large data sets to be in a data frame for its
#plotting function
diffeq_result_min <- as.data.frame(diffeq_result_min)
diffeq_result_max <- as.data.frame(diffeq_result_max)
#The R gather() function, extracts the data needed by ggplot.
                                        #It also helps to massage the data in a way that preps it to be plotted.
gathered_result_min <- gather(diffeq_result_min, variable, value, -time)
gathered_result_min$type <- "min"

gathered_result_max <- gather(diffeq_result_max, variable, value, -time)
gathered_result_max$type <- "max"

gathered_result <- rbind(gathered_result_min, gathered_result_max)

#This provides each variable, and therefore each graph, with a descriptive title.
gathered_result$variable[gathered_result$variable == "T"] <-  "CD4+ T-Cells"
gathered_result$variable[gathered_result$variable == "I"] <-  "HIV Infected T-Cells"
gathered_result$variable[gathered_result$variable == "VI"] <-  "Viral Load"
gathered_result$variable[gathered_result$variable == "VNI"] <-  "Non-infectious Virions"
gathered_result$variable[gathered_result$variable == "Z"] <-  "CD8+ T-Cells"
gathered_result$variable[gathered_result$variable == "ZA"] <- "Activated CD8+ T-Cells"

#The default order of ggplot2 is alphabetic, the factor() function enables 
#you to change the order of the plots into something more logical
#or relevant to your purposes. 
gathered_result$variable <- factor(gathered_result$variable, levels = c("CD4+ T-Cells", "HIV Infected T-Cells", "Viral Load",
                                                                        "Non-infectious Virions", "CD8+ T-Cells", "Activated CD8+ T-Cells"))

gathered_result_wider <- gathered_result %>%
    pivot_wider(
        id_cols = c("time", "variable"),
        names_from = type,
        values_from = value
        )



#Now plot the data and save it in a new variable: plot_result.
plot_result <- ggplot(data = gathered_result,
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


#Save image to file
fname_base <- "result_pi" #name of file
extensions = c("png", "pdf") #file extension
for (ext in extensions) {
    fname <- paste(fname_base, ".", ext, sep="")
    ggsave(fname, plot_result, height = 5, width = 9) #height, width of image
    print(paste("wrote output to file ", fname)) #If successful a message will print
    #with the name of the new file
    }
