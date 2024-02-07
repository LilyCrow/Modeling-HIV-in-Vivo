library(deSolve)
library(ggplot2)
library(tidyr)

# rti_func is the function on the right hand side of the
# differential equation. 
id_func <- function(t, ID, params) {

    T = ID["T"]
    I = ID["I"]
    V = ID["V"]

    lambda = params["lambda"] #Rate at which CD4+ T-cells are produced by the body
    d = params["d"]           #Natural death rate of uninfected CD4+ T-cells
    k = params["k"]           #Infection rate of CD4+ T-cells
    tau = params["tau"]       #Intracellular delay
    delta = params["delta"]   #Natural death rate of infected CD4+ T-cells
    N = params["N"]           #Viral burst size (# of virions produced per infected CD4+ T-cell)
    c = params["c"]           #Rate at which virions are cleared by the body
    m = params["m"]           #Rate at which infected cells die before producing virions


    if (t < tau){
        #t < tau cells still get infected but don't yet produce virions (take time for cell cycle)
        dT = lambda - d*T - k*T*V
        dI = k*T*V*exp(-m*tau) - delta*I
        #dV = ifelse(lagvalue(t - 1)[3] - c*V < 0, 0, -c*V)
        dV = -c*V
       
    }else if (t >= tau){
        lag_T = lagvalue(t - tau)[1]
        lag_V = lagvalue(t - tau)[3]
        
        dT = lambda - d*T - k*T*V
        dI = k*lag_T*lag_V*exp(-m*tau) - delta*I
        dV = N*delta*I - c*V
    }

    list(c(dT, dI, dV))
}

#Parameters
parms <- c(lambda = 10000,  #Rate at which CD4+ T-cells are produced by the body
           d = 0.01,        #Natural death rate of uninfected CD4+ T-cells
           k = 0.000000024, #Infection rate of CD4+ T-cells
           tau = 1.5,       #Intracellular delay
           delta = 1,       #Natural death rate of infected CD4+ T-cells
           N = 2500,        #Viral burst size (# of virions produced per infected CD4+ T-cell)
           c = 23,          #Rate at which virions are cleared by the body
           m = 0.01         #Rate at which infected cells die before producing virions
           )

#Initial values
Pstart <- c(T = 1000000,    #CD4+ T-cells
       	    I = 0,       #Infected CD4+ T-cells
            V = 0.01    #Viral particles
	    )    

#Time steps to solve equation for (full collection of x-axis points)
times <- seq(from=0, to=365, by= 1)


#Invoke the ODE function.
#It requires four things to run: the right hand side of the function
#(rti_func), the initial values (Pstart), time steps that need to be
#solved for (times), and equation constants
diffeq_result <- dede(
    func=id_func,
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
gathered_result$variable[gathered_result$variable == "T"] <- "CD4+ T-Cells"
gathered_result$variable[gathered_result$variable == "I"] <- "Infected CD4+ T-Cells"
gathered_result$variable[gathered_result$variable == "V"] <- "Viral Load"

#The default order of ggplot2 is alphabetic, the factor() function enables 
#you to change the order of the plots into something more logical
#or relevant to your purposes.
gathered_result$variable <- factor(gathered_result$variable, levels = c("CD4+ T-Cells", "Infected CD4+ T-Cells",
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
fname_base <- "result_id" #Name of file
extensions = c("png", "pdf")   #File extension

for (ext in extensions) {
    fname <- paste(fname_base, ".", ext, sep="")
    ggsave(fname, plot_result, height = 5, width = 9) #height, width of png/ pdf
    print(paste("wrote output to file ", fname))      #If successful a message
    #will print out the name of the new image 
    }
