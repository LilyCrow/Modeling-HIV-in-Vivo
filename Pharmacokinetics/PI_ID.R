library(deSolve)
library(ggplot2)
library(tidyr)

####################################################################
###TODO: FIND OUT HOW TO USE CURRENT TIME FOR INTRACELLULAR DELAY###
####################################################################

###############
###Ritonavir###
###############

pi_id_func <- function(t, PI_ID, params) {

    C_b = PI_ID["C_b"]
    C_c = PI_ID["C_c"]
    epsilon_PI = PI_ID["epsilon_PI"]
    T = PI_ID["T"]
    I = PI_ID["I"]
    V_I = PI_ID["V_I"]
    V_NI = PI_ID["V_NI"]

    F = params["F"]  
    D = params["D"] 
    V_d = params["V_d"] 
    k_a = params["k_a"] 
    k_e = params["k_e"]  
    I_d = params["I_d"]   
    N_d = params["N_d"]   
    f_b = params["f_b"]  
    H = params["H"]   
    IC_50_PI = params["IC_50_PI"] 
    lambda = params["lambda"]
    d = params["d"]           
    k = params["k"]          
    tau = params["tau"]       
    delta = params["delta"] 
    N = params["N"]         
    c = params["c"]      
    m = params["m"]          
    

    if (t < tau){
        C_b = 0        #Plasma concentration
        C_c = 0        #Intracellular concentration
        epsilon_PI = 0 #Instantatneous efficacy
        dT = lambda - d*T - k*T*V_I      #CD4+ T-cells
        dI = k*T*V_I*exp(-m*tau)         #Infected CD4+ T-cells
        dV_I = N*delta*I*(1 - epsilon_PI) - c*V_I    #Infectious Virions
        dV_NI = N*delta*I*epsilon_PI - c*V_NI        #Non-infectious Virions
    }else if (t >= tau){
        lag_T = lagvalue(t - tau)[4]
        lag_V_I = lagvalue(t - tau)[6]

            C_b = ((F*D) / (V_d)) * ((k_a)/(k_e - k_a)) * ((exp(-k_e*t)) / ((exp(k_a * I_d)) - 1)) * (1 - (exp(k_e - (k_a * t)))) * (1 - (exp(N_d * k_a * I_d))) + ((exp(k_e * I_d)) - (exp(k_a * I_d))) * ((exp((N_d - 1) * k_e * I_d) - 1) / ((exp(k_e * I_d)) - 1)) - (exp((((N_d - 1) * k_e) + k_a) * I_d))
            C_c = (1 - f_b)*H*C_b
            epsilon_PI = (C_c) / (IC_50_PI + C_c)
            dT = lambda - d*T - k*T*V_I
            dI = k*lag_T*lag_V_I*exp(-m*tau) - delta*I
            dV_I = N*delta*I*(1 - epsilon_PI) - c*V_I
            dV_NI = N*delta*I*epsilon_PI - c*V_NI
}
    list(c(C_b, C_c, epsilon_PI, dT, dI, dV_I, dV_NI))
}

#Parameters
parms <- c(F = 1,                #Fraction of drug that is absorbed
           D = 600,              #Mass of the drug administered in a dose
           V_d = 28000,          #Volume of distribution
           k_a = 14.64,          #Rate with which the drug is absorbed
           k_e = 6.86,           #Rate with which the drug is elminated from the blood
           I_d = 0.5,            #Dosing interval
           N_d = 731,            # number of doses until time t
           f_b = 0.99,           #fractionof drug that can't be transported into the cells
           H = 0.052,            #Quantified effect of the cell membrane, partition coefficient
           IC_50_PI = 0.0000009, #Plasma concentration that reduces viral production rate by 50%
           lambda = 10000,       #Rate at which CD4+ T-cells are produced by the body
           d = 0.01,             #Natural death rate of uninfected CD4+ T-cells
           k = 0.000000024,      #Infection rate of CD4+ T-cells
           tau = 3,           #Intracellular delay
           delta = 0.74,         #Natural death rate of infected CD4+ T-cells
           N = 1249,             #Viral burst size (# of virions produced per infected CD4+ T-cell)
           c = 23,               #Rate at which virions are cleared by the body
           m = 0.01              #Rate at which infected cells die before producing virions
           #t = times
           )

#Initial values
Pstart <- c(C_b = 0,        #Plasma concentration
            C_c = 0,        #Intracellular concentration
            epsilon_PI = 0, #Instantatneous efficacy
            T = 1000000,       #CD4+ T-cells
       	    I = 0,          #Infected CD4+ T-cells
            V_I = 0.001,    #Infectious Virions
            V_NI = 0        #Non-infectious Virions
	    )    

#Time steps to solve equation for (full collection of x-axis points)
times <- seq(from=0, to= 10, by= 1)


#Invoke the ODE function.
diffeq_result <- dede(
    func=pi_id_func,
    y=Pstart,
    times=times,
    parms=parms)

#Take the solved equations and put them into a data frame.
#ggplot2 requires such large data sets to be in a data frame for its
#plotting function
diffeq_result <- as.data.frame(diffeq_result) 

#The R gather() function, extracts the data needed by ggplot.
#It also helps to massage the data in a way that preps it to be plotted. 
gathered_result <- gather(diffeq_result, variable, value, -time)


#This provides each variable, and therefore each graph, with a descriptive title.
gathered_result$variable[gathered_result$variable == "C_b"] <- "Plasma (Blood) Concentration"
gathered_result$variable[gathered_result$variable == "C_c"] <- "Intracellular Concentration"
gathered_result$variable[gathered_result$variable == "epsilon_PI"] <- "Instantaneous Efficacy"
gathered_result$variable[gathered_result$variable == "T"] <- "CD4+ T-Cells"
gathered_result$variable[gathered_result$variable == "I"] <- "Infected CD4+ T-Cells"
gathered_result$variable[gathered_result$variable == "V_I"] <- "Infectious Virions" 
gathered_result$variable[gathered_result$variable == "V_NI"] <- "Non-infectious Virions"

#The default order of ggplot2 is alphabetic, the factor() function enables 
#you to change the order of the plots into something more logical
#or relevant to your purposes.
gathered_result$variable <- factor(gathered_result$variable, levels = c("Plasma (Blood) Concentration",
                                                                        "Intracellular Concentration",
                                                                        "Instantaneous Efficacy",
                                                                        "CD4+ T-Cells", "Infected CD4+ T-Cells",
                                                                        "Infectious Virions", "Non-infectious Virions"))

#Now plot the data and save it in a new variable: plot_result.
plot_result <- ggplot(data = gathered_result, 
                      mapping=aes(x=time, y = value, color = variable)) +
                      geom_line(linewidth = 1) + #Linewidth
                      theme_classic() + #ggplot has several different themes
                      #that can be added to a plot here we are using "classic"
                      theme(legend.position = "none") + #Locate legend at the bottom of the image
                      ylab("Concentration (mm3)") + #Label y-axis "Concentration (mm3)"
                      xlab("Time (days)") +  #Label x-axis "Time (days)"
                      facet_wrap(~variable, scales = "free_y") +  #plot each equation on its own graph
                      coord_cartesian(
                      xlim = NULL,
                      ylim = c(0, NA)) #Y-axis begins at 0
		

#Save image to file
fname_base <- "result_pi_id"  #Name of file
extensions = c("png", "pdf")   #File extension


for (ext in extensions) {
    fname <- paste(fname_base, ".", ext, sep="")
    ggsave(fname, plot_result, height = 9, width = 9) #height, width of png/ pdf
    print(paste("wrote output to file ", fname)) 
    }
