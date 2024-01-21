library(deSolve)
library(ggplot2)
library(tidyr)

####################################################################
###TODO: FIND OUT HOW TO USE CURRENT TIME FOR INTRACELLULAR DELAY###
####################################################################

##################
###TENOFOVIR DF###
##################


# rti_func is the function on the right hand side of the
# differential equation. 
rti_id_func <- function(t, RTI_ID, params) {

    C_b = RTI_ID["C_b"]
    C_x = RTI_ID["C_x"]
    C_c = RTI_ID["C_c"]
    C_cp = RTI_ID["C_cp"]
    C_cpp = RTI_ID["C_cpp"]
    epsilon_RTI = RTI_ID["epsilon_RTI"]
    T = RTI_ID["T"]
    I = RTI_ID["I"]
    V = RTI_ID["V"]

    F = params["F"]                 #Fraction of the drug that is absorbed into the blood
    D = params["D"]                 #Mass of the drug administered in a dose
    V_d = params["V_d"]             #Volume of distribuition (pt weight = 70 kg)
    k_a = params["k_a"]             #Rate at which drug is absorbed into the blood
    k_e = params["k_e"]             #Rate at which drug is eliminated from the blood
    I_d = params["I_d"]             #Dosing Interval
    N_d = params["N_d"]             #Number of doses until time t
    f_b = params["f_b"]             #Fraction of the drug that cannot be transported into cells
    H = params["H"]                 #Effect of the cell membrane, partition coefficient
    k_acell = params["k_acell"]     #Rate constant for cellular absorption
    k_ecell = params["k_ecell"]     #Rate constant for cellular elimination
    k_1f = params["k_1f"]           #Rate constant of phosphorylation to monophosphorylated form (C_cp)
    k_1b = params["k_1b"]           #Rate constant of dephosphorylation (C_cp to C_c)
    k_2f = params["k_2f"]           #Rate constant of phosphorylation to diphosphorylated from of drug (C_cp to C_cpp)
    k_2b = params["k_2b"]           #Rate constant of dephosphorylation (C_cpp to C_cp)
    IC_50_RTI = params["IC_50_RTI"] #Plasma concentration that reduces viral production rate by 50%

    lambda = params["lambda"] #Rate at which CD4+ T-cells are produced by the body
    d = params["d"]           #Natural death rate of uninfected CD4+ T-cells
    k = params["k"]           #Infection rate of CD4+ T-cells
    tau = params["tau"]       #Intracellular delay
    delta = params["delta"]   #Natural death rate of infected CD4+ T-cells
    N = params["N"]           #Viral burst size (# of virions produced per infected CD4+ T-cell)
    c = params["c"]           #Rate at which virions are cleared by the body
    m = params["m"]           #Rate at which infected cells die before producing virions
    #t = params["t"]

    # now that we have extracted the parameters, we can evaluate the
    # differential equation's right hand side

    if (t < tau){
        dC_b = ((F*D) / (V_d)) * ((k_a)/(k_e - k_a)) * ((exp(-k_e*t)) / ((exp(k_a * I_d)) - 1)) * (1 - (exp(k_e - (k_a * t)))) * (1 - (exp(N_d * k_a * I_d))) + ((exp(k_e * I_d)) - (exp(k_a * I_d))) * ((exp((N_d - 1) * k_e * I_d) - 1) / ((exp(k_e * I_d)) - 1)) - (exp((((N_d - 1) * k_e) + k_a) * I_d))
        C_x = (1 - f_b)*H*C_b - C_c
        dC_c = k_acell*C_x - k_ecell*C_c - k_1f*C_c + k_1b*C_cp
        dC_cp = -k_ecell*C_cp + k_1f*C_c - k_1b*C_cp - k_2f*C_cp + k_2b*C_cpp
        dC_cpp = -k_ecell*C_cpp + k_2f*C_cp - k_2b*C_cpp
        epsilon_RTI = (C_cpp) / (IC_50_RTI + C_cpp)
    
        dT = lambda - d*T - (1 - epsilon_RTI)*k*T*V
        dI = (1 - epsilon_RTI)*k*T*V*exp(-m*tau)
        dV = - c*V
    }else if(t >= tau){
        lag_epsilon_RTI = lagvalue(t - tau)[6]
        lag_T = lagvalue(t - tau)[7]
        lag_V = lagvalue(t - tau)[9]
        
        dC_b = ((F*D) / (V_d)) * ((k_a)/(k_e - k_a)) * ((exp(-k_e*t)) / ((exp(k_a * I_d)) - 1)) * (1 - (exp(k_e - (k_a * t)))) * (1 - (exp(N_d * k_a * I_d))) + ((exp(k_e * I_d)) - (exp(k_a * I_d))) * ((exp((N_d - 1) * k_e * I_d) - 1) / ((exp(k_e * I_d)) - 1)) - (exp((((N_d - 1) * k_e) + k_a) * I_d))
        C_x = (1 - f_b)*H*C_b - C_c
        dC_c = k_acell*C_x - k_ecell*C_c - k_1f*C_c + k_1b*C_cp
        dC_cp = -k_ecell*C_cp + k_1f*C_c - k_1b*C_cp - k_2f*C_cp + k_2b*C_cpp
        dC_cpp = -k_ecell*C_cpp + k_2f*C_cp - k_2b*C_cpp
        epsilon_RTI = (C_cpp) / (IC_50_RTI + C_cpp)
    
        dT = lambda - d*T - (1 - epsilon_RTI)*k*T*V
        dI = (1 - lag_epsilon_RTI)*k*lag_T*lag_V*exp(-m*tau) - delta*I
        dV = N*delta*I - c*V
        }

    list(c(dC_b, C_x, dC_c, dC_cp, dC_cpp, epsilon_RTI, dT, dI, dV))
}

#Parameters
parms <- c(F = 0.39,       #Fraction of the drug that is absorbed into the blood
           D = 300,            #Mass of the drug administered in a dose
           V_d = 87,500,    #Volume of distribuition (pt weight = 70 kg)
           k_a = 14.64,    #Rate at which drug is absorbed into the blood
           k_e = 9.60,     #Rate at which drug is eliminated from the blood
           I_d = 1,        #Dosing Interval
           N_d = 366,      #Number of doses until time t
           f_b = 0.07,     #fractionof the drug that cannot be transported into cells
           H = 1800,       #Effect of the cell membrane, partition coefficient
           k_acell = 6.624,#Rate constant for cellular absorption
           k_ecell = 1.1,  #Rate constant for cellular elimination
           k_1f = 9.6,     #Rate constant of phosphorylation to monophosphorylated form (C_cp)
           k_1b = 30.3,    #Rate constant of dephosphorylation (C_cp to C_c)
           k_2f = 270.7,   #Rate constant of phophorylation to diphosphorylated from of drug (C_cp to C_cpp)
           k_2b = 95.5,    #Rate constant of dephosphorylation (C_cpp to C_cp)
           IC_50_RTI = 0.54, #Plasma concentration that reduces viral production rate by 50%
           
           lambda = 10000, #Rate at which CD4+ T-cells are produced by the body
           d = 0.01,        #Natural death rate of uninfected CD4+ T-cells
           k = 0.000000024, #Infection rate of CD4+ T-cells
           tau = 0.34,       #Intracellular delay
           delta = 1.06,       #Natural death rate of infected CD4+ T-cells
           N = 1110,       #Viral burst size (# of virions produced per infected CD4+ T-cell)
           c = 23,          #Rate at which virions are cleared by the body
           m = 0.01         #Rate at which infected cells die before producing virions
           #t = times
           )

#Initial values
Pstart <- c(C_b = 0,     #Concentration of the drug in the blood
            C_x = 0,     #Driving force for drug transport
            C_c = 0,     #Intracellular concentration of native form of the drug
            C_cp = 0,    #Intracellular concentration of monophosphorylated form of the drug
            C_cpp = 0,   #Intracellular concentration of diphosphorylated form of the drug
            epsilon_RTI = 0, #Instantaneous efficacy of the drug
            T = 1000000,    #CD4+ T-cells
       	    I = 0,       #Infected CD4+ T-cells
            V = 0.001    #Viral particles
	    )    

#Time steps to solve equation for (full collection of x-axis points)
times <- seq(from=0, to=365, by= 1)


#Now inoke the ODE function.
#It requires four things to run: the right hand side of the function
#(rti_func), the initial values (Pstart), time steps that need to be
#solved for (times), and equation constants (parms_min or parms_max)

diffeq_result <- dede(
    func=rti_id_func,
    y=Pstart,
    times=times,
    parms=parms)

#Take the solved equations (diffeq_result_min and diffeq_result_max) and put them into a data frame.
#ggplot2 requires such large data sets to be in a data frame for its
#plotting function
diffeq_result <- as.data.frame(diffeq_result)

#The R gather() function, extracts the data needed by ggplot.
#It also helps to massage the data in a way that preps it to be plotted. 
gathered_result <- gather(diffeq_result, variable, value, -time)


#This provides each variable, and therefore each graph, with a descriptive title.
gathered_result$variable[gathered_result$variable == "C_b"] <- "Blood Concentration"
gathered_result$variable[gathered_result$variable == "C_x"] <- "Driving factor across cell membrane"
gathered_result$variable[gathered_result$variable == "C_c"] <- "Native Concentration"
gathered_result$variable[gathered_result$variable == "C_cp"] <- "Monophosphorylated Concentration"
gathered_result$variable[gathered_result$variable == "C_cpp"] <- "Diphosphorylated Concentration"
gathered_result$variable[gathered_result$variable == "epsilon_RTI"] <- "Instantaneous Efficacy"
gathered_result$variable[gathered_result$variable == "T"] <- "CD4+ T-Cells"
gathered_result$variable[gathered_result$variable == "I"] <- "Infected CD4+ T-Cells"
gathered_result$variable[gathered_result$variable == "V"] <- "Viral Load"

#The default order of ggplot2 is alphabetic, the factor() function enables 
#you to change the order of the plots into something more logical
#or relevant to your purposes.
gathered_result$variable <- factor(gathered_result$variable, levels = c("Blood Concentration",
                                                                        "Driving factor across cell membrane",
                                                                        "Native Concentration",
                                                                        "Monophosphorylated Concentration",
                                                                        "Diphosphorylated Concentration",
                                                                        "Instantaneous Efficacy", 
                                                                        "CD4+ T-Cells",
                                                                        "Infected CD4+ T-Cells",
                                                                        "Viral Load"))

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
fname_base <- "result_rti_id" #Name of file
extensions = c("png", "pdf")   #File extension


for (ext in extensions) {
    fname <- paste(fname_base, ".", ext, sep="")
    ggsave(fname, plot_result, height = 9, width = 9) #height, width of png/ pdf
    print(paste("wrote output to file ", fname))      #If successful a message
    #will print out the name of the new image 
    }
