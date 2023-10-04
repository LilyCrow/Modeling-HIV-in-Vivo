library(deSolve)
library(ggplot2)
library(tidyr)


# hiv_func puts the equations, initial values, and parameters
#into one variable so that the ODE can solve it.
hiv_func <- function(t, HIV, params) {

    #the following two blocks of code explain to R what
    #what the variables inside the function are (i.e. makes
    #the connection between the in func. variables and the
    #external variables that hold their values.

    #Dependent variables (what you are solving for)
    T = HIV["T"]    #CD4+ T-cells
    I = HIV["I"]    #Infected CD4+ T-cells
    V = HIV["V"]    #Viral load (number of viruses/ virions in blood)
    Z = HIV["Z"]    #CD8+ T-cells
    ZA = HIV["ZA"]  #Activated CD8+ T-cells


    #Independent variables (physical constants, etc.)
    mu_T = params["mu_T"]          #Natural death rate of uninfected CD4+ T-cells
    mu_I = params["mu_I"]          #Natural death rate of infected CD4+ T-cells
    mu_v = params["mu_v"]          #Natural death rate of viruses and virions
    mu_z = params["mu_z"]          #Natural death rate of CD8+ T-cells
    k_v = params["k_v"]            #Average number of virions released per infected cell
    beta_z = params["beta_z"]      #CD8+ T-cell activation rate
    beta_v = params["beta_v"]      #Rate at which CD4+ T-cells become infected
    p_I = params["p_I"]            #Rate at which infected CD4+ T-cells are destroyed
    p_v = params["p_v"]            #Rate at which viruses and virions are destroyed
    lambda_T = params["lambda_T"]  #Rate at which CD4+ T-cells are created by the body
    lambda_z = params["lambda_z"]  #Rate at which CD8+ T-cells are created by the body


    # now that we have extracted the parameters, we can evaluate the
    # differential equation's right hand side
    dT = lambda_T - mu_T*T - beta_v*T*V
    dI = beta_v*T*V - mu_I*I - p_I*I*ZA
    dV = k_v*mu_I*I - mu_v*V - p_v*V*ZA
    dZ = lambda_z - mu_z*Z - beta_z*Z*V
    dZA = beta_z*Z*I -mu_z*ZA

    #Put all 5 equations into a list so the ODE can read it  
    list(c(dT, dI, dV, dZ, dZA))
}

#Parameters
parms <- c(mu_T = 0.02,        #Natural death rate of uninfected CD4+ T-cells
           mu_I = 0.24,        #Natural death rate of infected CD4+ T-cells
           mu_v = 2.4,         #Natural death rate of viruses and virions
           mu_z = 0.04,        #Natural death rate of CD8+ T-cells
           k_v = 360,          #Average number of virions released per infected cell
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
            V = 0.001,  #Viral load 
	    Z = 500,    #CD8+ T-cells
	    ZA = 0)     #Activated CD8+ T-cells

#Time steps to solve equation for (full collection of x-axis points)
times <- seq(from=0, to=365, by= 0.01)



#Now inoke the ODE function.
#It requires four things to run: the right hand side of the function
#(hiv_func), the initial values (Pstart), time steps that need to be
#solved for (times), and equation constants (parms)
diffeq_result <- ode(
    func=hiv_func,
    y=Pstart,
    times=times,
    parms=parms)


#Take the solved equations (diffeq_result) and put it into a data frame
#ggplot2 requires such large data sets to be in a data frame for its
#plotting function
diffeq_result <- as.data.frame(diffeq_result)

#The R gather() function, extracts the data needed by ggplot.
#It also helps to massage the data in a way that preps it to be plotted. 
gathered_result <- gather(diffeq_result, variable, value, -time)


#This provides each variable, and therefore each graph, with a descriptive title.
gathered_result$variable[gathered_result$variable == "T"] <-  "CD4+ T-Cells"
gathered_result$variable[gathered_result$variable == "I"] <-  "Infected CD4+ T-Cells"
gathered_result$variable[gathered_result$variable == "V"] <-  "Viral Load"
gathered_result$variable[gathered_result$variable == "Z"] <-  "CD8+ T-Cells"
gathered_result$variable[gathered_result$variable == "ZA"] <- "Activated CD8+ T-Cells"

#The default order of ggplot2 is alphabetic, the factor() function enables 
#you to change the order of the plots into something more logical
#or relevant to your purposes. 
gathered_result$variable <- factor(gathered_result$variable, levels = c("CD4+ T-Cells", "Infected CD4+ T-Cells",
                                                                        "Viral Load","CD8+ T-Cells", "Activated CD8+ T-Cells"))




shinyPlot <- function(choicesVar){#), choicesDT){

   ### if(is.null(choicesDT)){
#Now plot the data and save it in a new variable: plot_result.
    plot_result <- ggplot(data = subset(gathered_result, variable %in% choicesVar),
                      mapping=aes(x=time, y=value, color = variable)) +
                      geom_line(linewidth = 1) + #linewidth
                      theme_classic() + #ggplot has several different themes
                      #that can be added to a plot, here we are using "classic"
                      theme(legend.position = "none") +
                      ylab("Concentration (cells/mm3)")+ #Y-axis title
                      xlab("Time (days)")+ #X-axis title
                      facet_wrap(~variable, scales = "free_y") + #plot each equation on its own graph
                      coord_cartesian(
                          xlim = NULL,
                          ylim = c(0, NA)) #begin each plot at y = 0
    return(plot_result)#}
}




