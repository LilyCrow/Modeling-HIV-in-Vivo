require(deSolve)
require(ggplot2)
require(gridExtra)
require(dplyr)

# parameters
v_HIV_parm1 = c(s =6 , p = 2, tmax = 54, dt = 1/7, k =.3 , delta = .49, N = 10.3e9, c = 3.07)

# initial values
v_HIV_init_val <- c(T = 0.999999, I = 1e-6, V = 0)

t <- seq(0, 200, by = 5)


#Model Equations
model.HIV <- function(initval, params, t, ...){
    hiv <- lsoda(y = initval, times = t,
                 func = function(time, state_vars, param){ 
                     with(as.list(c(state_vars, param)), {
                         
                         dT <- s + p*T-dt*T-k*V*T #equation for proportion of susceptible population
                         dI <- k*V*T - delta*I #equation for proportion of infected population
                         dV <- N*delta*I - c*V #equation for proportion of recovered population
                         
                         return(list(c(dT, dI, dV)))
                     })
                 },
                 parms = params, ...)

    print(paste(hiv))
    return(as.data.frame(hiv))

}

df_hiv_model1 <- model.HIV(v_HIV_init_val, v_HIV_parm1, 0:150)
plot(df_hiv_model1, main = "Modeling HIV in Vivo")#, xlab = "Time", ylab = "Percent of Population" )
#par(new=TRUE)                                  

        
        

                         
                         
