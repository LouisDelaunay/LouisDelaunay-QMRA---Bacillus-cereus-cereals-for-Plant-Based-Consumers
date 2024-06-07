#June 2024 

#######################
#---- Input files ----#
#######################
    # R_input_data.xlsx

########################
#---- Output files ----#
########################
    # Cooling_time.RData
    # T_frigde.RData

    # gamma_fridge_Bc3.RData
    # gamma_fridge_Bc4.RData

    # Qty_cat16_day.RData
    # dish_size.RData
    # Qty_cat16_tot.RData

    # N0_tot_C.RData
    # prev_Bernoulli_C.RData
    # Qty_Bc_t0_C.RData
    # N0_tot_Bc3_C.RData
    # N0_tot_Bc4_C.RData
    # N_ht_Bc3_C.RData
    # N_ht_Bc4_C.RData
    # N_cooling_time_Bc3_C.RData
    # N_24h_nofridge_Bc3_C.RData
    # N_cooling_time_Bc4_C.RData
    # N_24h_nofridge_Bc4_C.RData
    # N_24h_Bc3_C.RData
    # N_48h_Bc3_C.RData
    # N_72h_Bc3_C.RData
    # N_96h_Bc3_C.RData
    # N_120h_Bc3_C.RData
    # N_24h_Bc4_C.RData
    # N_48h_Bc4_C.RData
    # N_72h_Bc4_C.RData
    # N_96h_Bc4_C.RData
    # N_120h_Bc4_C.RData

    # Model_cereal_full_200000.RData

######################
#---- References ----#
######################
    # Albaridi, 2022 - Risk of Bacillus cereus contamination in cooked rice - Cienc Technol Aliment. 2022 Mar 18;42:e108221 - DOI: 10.1590/fst.108221
    # Anses, 2021 - Fiche de description de danger biologique transmissible par les aliments - Bacillus cereus; Saisine n°2016-SA-0075; 
    # Daelman et al., 2013 - Development of a time-to-detect growth model for heat-treated Bacillus cereus spores -  Int J Food Microbiol. 2013 Aug 1;165(3):231-40 - DOI: 10.1016/j.ijfoodmicro.2013.04.018
    # Ellouze et al., 2021 - Modeling Bacillus cereus Growth and Cereulide Formation in Cereal-, Dairy-, Meat-, Vegetable-Based Food and Culture Medium -  Front Microbiol. 2021 Feb 17:12:639546 - DOI: 10.3389/fmicb.2021.639546
    # Laurent, Arino, Rosso, 1999 - A quantitative approach for studying the effect of heat treatment conditions on resistance and recovery of Bacillus cereus spores -  Int J Food Microbiol. 1999 May 1;48(2):149-57 - DOI: 10.1016/S0168-1605(99)00039-2
    # Roccato et al., 2017 - Analysis of domestic refrigerator temperatures and home storage time distributions for shelf-life studies and food safety risk assessment - Food Res Int. 2017 Jun:96:171-181 - DOI: 10.1016/j.foodres.2017.02.017


# Cleaning the environment
rm(list=ls(all=TRUE))
gc()

library(mc2d)
library(fitdistrplus)
library(readxl)
library(nlstools)

setwd("C:/Users/XXXX")
getwd

ndvar(200000) # number of iterations in variability dimension
ndunc(1000) # number of iterations in uncertainty dimension

#-----------------------FIXED PARAMETERS-----------------------#

# Amount of B. cereus associated to a risk of food poisoning
mean_portion_size <- 182.96
cutoff            <- 10^5 #Anses, 2021: food poisoning usually observed with foods containing at least 10^5 CFU.g-1

# Heat treatment (cooking) and cooling (ambient) temperatures
Tht       <- 95
Tcooling  <- 19

# B. cereus group III
Tmin_Bc3    <- 7.3
Topt_Bc3    <- 39.06
Tmax_Bc3    <- 47.76
muopt_Bc3_C <- 1.41 # µopt for cereals (Ellouze et al., 2021 - see table 3) 
muopt_Bc3_P <- 2.05 # µopt for tubers and lentils - Combase ID: Car_45 ; pH 6.5 ; temperature 30; G= 0.346 µ = 2.05
z_Bc3       <- 6.2
D90_Bc3     <- 101.9
Dht_Bc3 = D90_Bc3 * 10^((90-Tht)/z_Bc3)

# B. cereus group IV
Tmin_Bc4    <- 7.44
Topt_Bc4    <- 38.71
Tmax_Bc4    <- 47.44
muopt_Bc4_C <- 1.41 # µopt for cereals (Ellouze et al., 2021 - see table 3) 
muopt_Bc4_P <- 2.05 # µopt for tubers and lentils - Combase ID: Car_45 ; pH 6.5 ; temperature 30; G= 0.346 µ = 2.05
z_Bc4       <- 7.128
D100_Bc4    <- 6.08
Dht_Bc4 = D100_Bc4 * 10^((100-Tht)/z_Bc4)


#--------------------------------------------------------------#


###########################
#   Heat treatment time   #
###########################

# Heat treatment time - 8 to 12 minutes for cereal-based dishes
min_tht_C <- 8
max_tht_C <- 12
tht_C     <- mcstoc(runif, type="V", min=min_tht_C, max=max_tht_C)

# Heating treatment time - 25 to 30 minutes for (tuber- and) pulse-based dishes #FRENCH - on pourrait supprimer cette information qui n'est pas utilisée par la suite
min_tht_P <- 25
max_tht_P <- 30
tht_P     <- mcstoc(runif, type="V", min=min_tht_P, max=max_tht_P)


#####################
#   Cooling time    #
#####################

Data_cooling_time         <- read_excel("R_input_data.xlsx", sheet="Cooling_time") # Data from INCA3
descdist(as.numeric(Data_cooling_time$Cooling_time), discrete=FALSE, method="unbiased")
Cooling_time_fit          <- fitdist(Data_cooling_time$Cooling_time,"exp")

# Summary(Cooling_time)
Boot_cooling_time  <- bootdist(Cooling_time_fit, bootmethod="nonpara", niter=ndunc())

# Boot_cooling_time$CI
rate_R        <- mcdata(Boot_cooling_time$estim$rate, type="U")
Cooling_time  <- mcstoc(rexp, type="VU", rate=rate_R)

# Cooling time, truncating the distribution at 24h (more realistic)
Cooling_time  <- ifelse(Cooling_time > 24, 24, Cooling_time) 




###########################
#   FRIDGE TEMPERATURE    #
###########################

# Data from Roccato et al., 2017  (data from INCA3 deemed not representative enough - no enough fridge at +10°C)
Data_temp_fridge  <- read_excel("R_input_data.xlsx", sheet="Roccato_temp_fridge")

# Min and max mean temperature - based on data from Northern countries
min_mean_temp <- quantile(subset(Data_temp_fridge, zone=='N')$mean, 0.25)
max_mean_temp <- quantile(subset(Data_temp_fridge, zone=='N')$mean, 0.75)

# Min and max sd temperature - based on data from Northern countries 
min_sd_temp   <- quantile(subset(Data_temp_fridge, zone=='N')$sd, 0.25)
max_sd_temp   <- quantile(subset(Data_temp_fridge, zone=='N')$sd, 0.75)

#Temperature distributions
T_ref_mean    <- mcstoc(runif, type = "U", min=min_mean_temp, max=max_mean_temp) 
T_ref_sd      <- mcstoc(runif, type = "U", min=min_sd_temp, max=max_sd_temp)
T_frigde      <- mcstoc(rnorm, type = "VU", mean=T_ref_mean, sd=T_ref_sd)

# Truncating the distribution so that temperatures do not fall below 0°C or exceed 12°C
T_fridge_min  <- 0
T_fridge_max  <- 15
T_frigde      <- ifelse(T_frigde < T_fridge_min, T_fridge_min, T_frigde)
T_frigde      <- ifelse(T_frigde > T_fridge_max, T_fridge_max, T_frigde)




#######################
#   Gamma in fridge   #
#######################
gamma_fridge_Bc3 <- ifelse(T_frigde<Tmin_Bc3, 0, ifelse(T_frigde>Tmax_Bc3, 0, ((T_frigde-Tmin_Bc3)^2*(T_frigde-Tmax_Bc3)/((Topt_Bc3-Tmin_Bc3)*((Topt_Bc3-Tmin_Bc3)*(T_frigde-Topt_Bc3)-(Topt_Bc3-Tmax_Bc3)*(Topt_Bc3+Tmin_Bc3-2*T_frigde))))))
gamma_fridge_Bc4 <- ifelse(T_frigde<Tmin_Bc4, 0, ifelse(T_frigde>Tmax_Bc4, 0, ((T_frigde-Tmin_Bc4)^2*(T_frigde-Tmax_Bc4)/((Topt_Bc4-Tmin_Bc4)*((Topt_Bc4-Tmin_Bc4)*(T_frigde-Topt_Bc4)-(Topt_Bc4-Tmax_Bc4)*(Topt_Bc4+Tmin_Bc4-2*T_frigde))))))


#----------- Model 1: Cereal-based dishes (pasta, rice, etc.) B. cereus -----------#

##############----------------------------------------------------

# Portion size of corresponding food category (16: Composite meal with no meat or fish and a majority of starchy foods) 
Qty_cat16_input       <- read_excel("R_input_data.xlsx", sheet="Quantity")
descdist(as.numeric(Qty_cat16_input$Quantity_cat16), discrete=FALSE, method="unbiased")
Qty_cat16_fit         <- fitdist(Qty_cat16_input$Quantity_cat16, "lnorm")
Boot_Qty_cat16_fit    <- bootdist(Qty_cat16_fit, bootmethod="nonpara", niter=ndunc())
mean_log_Qty_cat16    <- mcdata(Boot_Qty_cat16_fit$estim$meanlog, type="U")
sd_log_Qty_cat16      <- mcdata(Boot_Qty_cat16_fit$estim$sdlog, type="U")



Qty_cat16             <- mcstoc(rlnorm,type="VU", meanlog=mean_log_Qty_cat16, sdlog=sd_log_Qty_cat16)
Qty_cat16_max         <- ifelse(Qty_cat16 >3000, 3000, Qty_cat16) #Truncating data so that max 3 kgs consumed in 3 days
Qty_cat16_day         <- Qty_cat16_max/3 #INCA3 consumptions correspond to 3 24h-recalls - dividing by 3 gives the amount consumed per day
Nb_portions_per_dish  <- 6 
dish_size             <- Qty_cat16_day*Nb_portions_per_dish



#------------------------------- Results and heat inactivation treatment #FRENCH - pas sure que ce soit le meilleur titre mais je laisse pour le moment

#################################################################################################
##                                                                                             ##
##------------------- MODULE 1: Initial contamination before heat treatment -------------------##
##                                                                                             ##
#################################################################################################

# Data from Albaridi, 2022 
N0          <- 20  # Number of uncontaminated samples
N_inf_2log  <- 19  # Number of samples containing < 2 log CFU.g-1 of B. cereus
N_sup_2log  <- 1   # Number of samples containing 2 to 3 log CFU.g-1 of B. cereus
Ntot        <- N0 + N_inf_2log + N_sup_2log  # Total number of collected samples
# No sample was containing more than 3 log CFU.g-1

# Generating probabilities of B. cereus 
  # Dishes contaminated with less than 2 log CFU.g-1
N0_Bc_inf_2log_C  <- runif(100, min=1, max=100*dish_size) #quantity of B. cereus per dish
pct_Bc_inf_2log_C <- (N_inf_2log/(N_inf_2log + N_sup_2log))
P_Bc_inf_2log_C   <- rep(pct_Bc_inf_2log_C, 100)

  # Dishes contaminated with 2 to 3 log CFU.g-1
N0_Bc_sup_2log_C  <- runif(100, min=100*dish_size, max=1000*dish_size) #quantity of B. cereus per dish
pct_Bc_sup_2log_C <- (N_sup_2log/(N_inf_2log + N_sup_2log))
P_Bc_sup_2log_C   <- rep(pct_Bc_sup_2log_C, 100)

N0_tot_C          <- mcstoc(rempiricalD, type="V", values=c(N0_Bc_inf_2log_C,N0_Bc_sup_2log_C), prob=c(P_Bc_inf_2log_C,P_Bc_sup_2log_C)) #quantity of B. cereus (all groups combined)



# Calculating prevalence and initial quantities of group III & group IV B. cereus before heat treatment
## Prevalence and initial quantity of B. cereus
prev_t0_C         <- mcstoc(rbeta, type="U", shape1=N_inf_2log + N_sup_2log+1, shape2=Ntot-(N_inf_2log + N_sup_2log)+1, ncp=0) #prevalence 


##################################################################################################
##                                                                                              ##
##------------------- MODULE 2: Reduction of contamination by heat treatment -------------------##
##                                                                                              ##
##################################################################################################

# Survival of bacterial spores from group III B. cereus after heat treatment
p_surv_Bc3_C  <- 10^(-tht_C/Dht_Bc3) #survival probability of bacterial spores from group III B. cereus

# Survival of bacterial spores from group IV B. cereus after heat treatment
p_surv_Bc4_C  <- 10^(-tht_C/Dht_Bc4) #survival probability of bacterial spores from group IV B. cereus


#################################################################################################
##                                                                                             ##
##------------------- MODULE 3: Recontamination during cooling and storage  -------------------##
##                                                                                             ##
#################################################################################################

# Data from Laurent, Arino, Rosso, 1999 - see Figure 3, (b): 
lag_min_cooling <- 6
lag_max_cooling <- 8

# Lag time during cooling at ambient temperature
lag_cooling <- mcstoc(runif, type="U", min=lag_min_cooling, max=lag_max_cooling)

# Calculating gamma and quantity of bacteria after cooling time for group III B. cereus
gamma_cooling_Bc3_C   <- ifelse(Tcooling<Tmin_Bc3, 0, ifelse(Tcooling>Tmax_Bc3, 0, ((Tcooling-Tmin_Bc3)^2*(Tcooling-Tmax_Bc3)/((Topt_Bc3-Tmin_Bc3)*((Topt_Bc3-Tmin_Bc3)*(Tcooling-Topt_Bc3)-(Topt_Bc3-Tmax_Bc3)*(Topt_Bc3+Tmin_Bc3-2*Tcooling))))))

# Calculating gamma and quantity of bacteria after cooling time for group IV B. cereus
gamma_cooling_Bc4_C   <- ifelse(Tcooling<Tmin_Bc4, 0, ifelse(Tcooling>Tmax_Bc4, 0, ((Tcooling-Tmin_Bc4)^2*(Tcooling-Tmax_Bc4)/((Topt_Bc4-Tmin_Bc4)*((Topt_Bc4-Tmin_Bc4)*(Tcooling-Topt_Bc4)-(Topt_Bc4-Tmax_Bc4)*(Topt_Bc4+Tmin_Bc4-2*Tcooling))))))


# Data from Daelman et al., 2013 - see Figure 4: 
lag_min_fridge <- 12
lag_max_fridge <- 72

# Lag time during storage in the fridge
lag_storage <- mcstoc(runif, type="U", min=lag_min_fridge, max=lag_max_fridge)

#---- Growth during storage in the fridge
T_lagtime_gt1week <- 10 # Lag time was assumed to be more than one week if the refrigerator temperature was low enough. 


#------------------------------- Data for sensitivity analysis

#######################################################################
#   Percentage of Individuals keeping dishes at ambient temperature   #
#######################################################################
No_Indiv_wo_fridge  <- 21
No_Indiv_w_fridge   <- 161
Tot_No_Indiv        <- No_Indiv_wo_fridge + No_Indiv_w_fridge
Pct_wo_fridge       <- mcstoc(rpert, type="U", min=0, mode=0.01, max=No_Indiv_wo_fridge/Tot_No_Indiv) #mode set at 1% of the population
Pct_w_fridge        <- 1 - Pct_wo_fridge


#--- Risks per portion

#Data from a internal survey
data_storage_leftovers <- read_excel("R_input_data.xlsx", sheet="Storage_leftovers")

#Number of individuals storing leftovers from cereal-based dished for: 
storage_cereal      <- subset(data_storage_leftovers, dish_type=="cereal")

#No leftovers
nb_noleftovers_cons_C<-subset(storage_cereal, storage_in_days=="0")$number_indiv
# 1 day (i.e 24h) (max)
nb_24h_cons_C       <- subset(storage_cereal, storage_in_days=="1")$number_indiv
# 2 days (i.e. 48h) (max)
nb_48h_cons_C       <- subset(storage_cereal, storage_in_days=="2")$number_indiv 
# 3 days (i.e. 72h) (max)
nb_72h_cons_C       <- subset(storage_cereal, storage_in_days=="3")$number_indiv
# 4 days (i.e. 96h) (max) - option "more than 3 days" considered to be 4 or 5 days, with the same probability
nb_96h_cons_C       <- subset(storage_cereal, storage_in_days==">3")$number_indiv/2  
# 5 days (i.e. 120h) (max) - option "more than 3 days" considered to be 4 or 5 days, with the same probability
nb_120h_cons_C      <- subset(storage_cereal, storage_in_days==">3")$number_indiv/2

Sum_nb_cons_C       <- nb_noleftovers_cons_C + nb_24h_cons_C + nb_48h_cons_C + nb_72h_cons_C + nb_96h_cons_C + nb_120h_cons_C


prev_no_leftovers 	<- mcstoc(rbeta, type="U", shape1=nb_noleftovers_cons_C + 1, shape2=Sum_nb_cons_C-nb_noleftovers_cons_C + 1, ncp=0) 
prev_24h_cons_C     <- mcstoc(rbeta, type="U", shape1=nb_24h_cons_C+1, shape2=Sum_nb_cons_C-nb_24h_cons_C+1, ncp=0)  
prev_48h_cons_C     <- mcstoc(rbeta, type="U", shape1=nb_48h_cons_C+1, shape2=Sum_nb_cons_C-nb_48h_cons_C+1, ncp=0)  
prev_72h_cons_C     <- mcstoc(rbeta, type="U", shape1=nb_72h_cons_C+1, shape2=Sum_nb_cons_C-nb_72h_cons_C+1, ncp=0)  
prev_96h_cons_C     <- mcstoc(rbeta, type="U", shape1=nb_96h_cons_C+1, shape2=Sum_nb_cons_C-nb_96h_cons_C+1, ncp=0)
prev_120h_cons_C    <- mcstoc(rbeta, type="U", shape1=nb_120h_cons_C+1, shape2=Sum_nb_cons_C-nb_120h_cons_C+1, ncp=0)


Resultat_C<- mcmodel({


##############----------------------------------------------------

#################################################################################################
##                                                                                             ##
##------------------- MODULE 1: Initial contamination before heat treatment -------------------##
##                                                                                             ##
#################################################################################################



# Calculating prevalence and initial quantities of group III & group IV B. cereus before heat treatment
## Prevalence and initial quantity of B. cereus

N0_tot_C          <- mcstoc(rempiricalD, type="V", values=c(N0_Bc_inf_2log_C,N0_Bc_sup_2log_C), prob=c(P_Bc_inf_2log_C,P_Bc_sup_2log_C)) #quantity of B. cereus (all groups combined)
prev_Bernoulli_C  <- mcstoc(rbern, type="VU", prev_t0_C)

Qty_Bc_t0_C       <- N0_tot_C*prev_Bernoulli_C 


# Generating proportions of group III and group IV B. cereus
Pct_Bc3_C <- mcstoc(runif, type="U", min=0.2, max=0.8)
Pct_Bc4_C <- 1-Pct_Bc3_C

# Calculating corresponding amounts of bacteria
N0_tot_Bc3_C <- Qty_Bc_t0_C*Pct_Bc3_C
N0_tot_Bc4_C <- Qty_Bc_t0_C*Pct_Bc4_C 


##################################################################################################
##                                                                                              ##
##------------------- MODULE 2: Reduction of contamination by heat treatment -------------------##
##                                                                                              ##
##################################################################################################

# Survival of bacterial spores from group III B. cereus after heat treatment
N_surv_Bc3_C  <- N0_tot_Bc3_C*p_surv_Bc3_C
N_ht_Bc3_C    <- mcstoc(rpois, type="VU", lambda=N_surv_Bc3_C)

# Survival of bacterial spores from group IV B. cereus after heat treatment
N_surv_Bc4_C  <- N0_tot_Bc4_C*p_surv_Bc4_C
N_ht_Bc4_C    <- mcstoc(rpois, type="VU", lambda=N_surv_Bc4_C)


#################################################################################################
##                                                                                             ##
##------------------- MODULE 3: Recontamination during cooling and storage  -------------------##
##                                                                                             ##
#################################################################################################

N_cooling_time_Bc3_C  <- ifelse(Cooling_time<lag_cooling, N_ht_Bc3_C, N_ht_Bc3_C*exp(muopt_Bc3_C*gamma_cooling_Bc3_C*(Cooling_time-lag_cooling))) #Quantity of group III B. cereus after cooling
N_cooling_time_Bc4_C  <- ifelse(Cooling_time<lag_cooling, N_ht_Bc4_C, N_ht_Bc4_C*exp(muopt_Bc4_C*gamma_cooling_Bc4_C*(Cooling_time-lag_cooling))) #Quantity of group IV B. cereus after cooling 


#--- Growth during storage at ambient temperature (i.e. no fridge)
###------ 24h growth for group III B. cereus
N_24h_nofridge_Bc3_C <- N_ht_Bc3_C*exp(muopt_Bc3_C*gamma_cooling_Bc3_C*(24-lag_cooling))

###------ 24h growth for group IV B. cereus
N_24h_nofridge_Bc4_C <- N_ht_Bc4_C*exp(muopt_Bc4_C*gamma_cooling_Bc4_C*(24-lag_cooling)) 


# Lag time during storage in the fridge
lag_storage <- mcstoc(runif, type="U", min=lag_min_fridge, max=lag_max_fridge)

#---- Growth during storage in the fridge
T_lagtime_gt1week <- 10 # Lag time was assumed to be more than one week if the refrigerator temperature was low enough. 
##----- Group III B. cereus
###------ 24h growth
N_24h_Bc3_C   <- ifelse(T_frigde<T_lagtime_gt1week, N_cooling_time_Bc3_C, ifelse(24<lag_storage, N_cooling_time_Bc3_C, N_cooling_time_Bc3_C*exp(muopt_Bc3_C*gamma_fridge_Bc3*(24-lag_storage)))) 
###------ 48h growth
N_48h_Bc3_C   <- ifelse(T_frigde<T_lagtime_gt1week, N_cooling_time_Bc3_C, ifelse(48<lag_storage, N_cooling_time_Bc3_C, N_cooling_time_Bc3_C*exp(muopt_Bc3_C*gamma_fridge_Bc3*(48-lag_storage))))
###------ 72h growth
N_72h_Bc3_C   <- ifelse(T_frigde<T_lagtime_gt1week, N_cooling_time_Bc3_C, ifelse(72<lag_storage, N_cooling_time_Bc3_C, N_cooling_time_Bc3_C*exp(muopt_Bc3_C*gamma_fridge_Bc3*(72-lag_storage)))) 
###------ 96h growth
N_96h_Bc3_C   <- ifelse(T_frigde<T_lagtime_gt1week, N_cooling_time_Bc3_C, ifelse(96<lag_storage, N_cooling_time_Bc3_C, N_cooling_time_Bc3_C*exp(muopt_Bc3_C*gamma_fridge_Bc3*(96-lag_storage)))) 
###------ 120h growth
N_120h_Bc3_C  <- ifelse(T_frigde<T_lagtime_gt1week, N_cooling_time_Bc3_C, ifelse(120<lag_storage, N_cooling_time_Bc3_C, N_cooling_time_Bc3_C*exp(muopt_Bc3_C*gamma_fridge_Bc3*(120-lag_storage)))) 


##----- Group IV B. cereus
###------ 24h growth
N_24h_Bc4_C   <- ifelse(T_frigde<T_lagtime_gt1week, N_cooling_time_Bc4_C, ifelse(24<lag_storage, N_cooling_time_Bc4_C, N_cooling_time_Bc4_C*exp(muopt_Bc4_C*gamma_fridge_Bc4*(24-lag_storage)))) 
###------ 48h growth
N_48h_Bc4_C   <- ifelse(T_frigde<T_lagtime_gt1week, N_cooling_time_Bc4_C, ifelse(48<lag_storage, N_cooling_time_Bc4_C, N_cooling_time_Bc4_C*exp(muopt_Bc4_C*gamma_fridge_Bc4*(48-lag_storage)))) 
###------ 72h growth
N_72h_Bc4_C   <- ifelse(T_frigde<T_lagtime_gt1week, N_cooling_time_Bc4_C, ifelse(72<lag_storage, N_cooling_time_Bc4_C, N_cooling_time_Bc4_C*exp(muopt_Bc4_C*gamma_fridge_Bc4*(72-lag_storage)))) 
###------ 96h growth
N_96h_Bc4_C   <- ifelse(T_frigde<T_lagtime_gt1week, N_cooling_time_Bc4_C, ifelse(96<lag_storage, N_cooling_time_Bc4_C, N_cooling_time_Bc4_C*exp(muopt_Bc4_C*gamma_fridge_Bc4*(96-lag_storage)))) 
###------ 120h growth
N_120h_Bc4_C  <- ifelse(T_frigde<T_lagtime_gt1week, N_cooling_time_Bc4_C, ifelse(120<lag_storage, N_cooling_time_Bc4_C, N_cooling_time_Bc4_C*exp(muopt_Bc4_C*gamma_fridge_Bc4*(120-lag_storage)))) 

##----- Concentration of B. cereus per g
C0_tot_Bc_C<-(N0_tot_Bc3_C+N0_tot_Bc4_C)/Nb_portions_per_dish/mean_portion_size
C_ht_Bc_C<-(N_ht_Bc3_C+N_ht_Bc4_C)/Nb_portions_per_dish/mean_portion_size
C_cooling_time_Bc_C<-(N_cooling_time_Bc3_C + N_cooling_time_Bc4_C)/Nb_portions_per_dish/mean_portion_size
C_24h_Bc_C<-(N_24h_Bc3_C+N_24h_Bc4_C)/Nb_portions_per_dish/mean_portion_size
C_48h_Bc_C<-(N_48h_Bc3_C+N_48h_Bc4_C)/Nb_portions_per_dish/mean_portion_size
C_72h_Bc_C<-(N_72h_Bc3_C+N_72h_Bc4_C)/Nb_portions_per_dish/mean_portion_size
C_96h_Bc_C<-(N_96h_Bc3_C+N_96h_Bc4_C)/Nb_portions_per_dish/mean_portion_size
C_120h_Bc_C<-(N_120h_Bc3_C+N_120h_Bc4_C)/Nb_portions_per_dish/mean_portion_size
C_24h_nofridge_Bc_C<-(N_24h_nofridge_Bc3_C+N_24h_nofridge_Bc4_C)/Nb_portions_per_dish/mean_portion_size


#For figure 1
concentration_C_5days        <- C_120h_Bc_C
log_concentration_0_C_5days  <- ifelse(concentration_C_5days<=0.0001, 0, log10(concentration_C_5days))
log_concentration_08_C_5days <- ifelse(log_concentration_0_C_5days >8, 8, log_concentration_0_C_5days)


#---- Calculating risks
risk_hot_C          <- C_ht_Bc_C > cutoff 

risk_24h_C          <- C_24h_Bc_C > cutoff 
risk_48h_C          <- C_48h_Bc_C > cutoff 
risk_72h_C          <- C_72h_Bc_C	> cutoff
risk_96h_C          <- C_96h_Bc_C > cutoff 
risk_120h_C         <- C_120h_Bc_C > cutoff
risk_24h_nofridge_C <- C_24h_nofridge_Bc_C > cutoff

risk_portion_fridge_C     <- (1-Pct_wo_fridge)*(prev_24h_cons_C*risk_24h_C + prev_48h_cons_C*risk_48h_C + prev_72h_cons_C*risk_72h_C + prev_96h_cons_C*risk_96h_C + prev_120h_cons_C*risk_120h_C)



mc(prev_Bernoulli_C, T_frigde, Cooling_time, Qty_cat16_day,
  N0_tot_C,
  N0_tot_Bc_C, N0_tot_Bc3_C, N0_tot_Bc4_C, 
  N_ht_Bc_C, N_ht_Bc3_C, N_ht_Bc4_C, 
  N_cooling_time_Bc_C, N_cooling_time_Bc3_C, N_cooling_time_Bc4_C, 
  N_24h_Bc_C, N_24h_Bc3_C, N_24h_Bc4_C, 
  N_48h_Bc_C, N_48h_Bc3_C, N_48h_Bc4_C,
  N_72h_Bc_C,N_72h_Bc3_C, N_72h_Bc4_C,
  N_96h_Bc_C, N_96h_Bc3_C, N_96h_Bc4_C,
  N_120h_Bc_C, N_120h_Bc3_C, N_120h_Bc4_C,
  N_24h_nofridge_Bc_C, N_24h_nofridge_Bc3_C, N_24h_nofridge_Bc4_C,
  risk_hot_C, risk_24h_C, risk_48h_C, risk_72h_C, risk_96h_C, risk_120h_C, risk_portion_no_fridge_C, risk_portion_hot_C, risk_portion_fridge_C,risk_24h_nofridge_C,
  log_concentration_08_C_5days)

})

rm(Qty_cat16,Qty_cat16_max)

Resu_C        <- evalmcmod(Resultat_C, nsv=ndvar(), nsu=ndunc())
summary(Resu_C, probs=c(0.5, 0.95, 0.99), lim=c(0.025,0.975))

#------------------------------- Sensitivity analysis (only for variability) -------------------------------#

# Plotting results in log CFU.g-1
plot(Resu_C$log_concentration_08_C_5days, ylab="Cumulative probability density - 5-day chilled storage", xlab ="log CFU.g-1", ymin=0, xlim=c(0,6), xmax=6) #main="Bacterial concentration in cereal-based dishes (CFU.g-1)", adj = 0.5, line = 6
abline(v=5, col="red", lty=2, lwd=3)

# Sensitivity Analysis
Ana_sens_C120 <- tornado(Resu_C, output="log_concentration_08_C_5days", use = "pairwise.complete.obs", method = c("spearman"), lim=c(0.025, 0.975))
plot(Ana_sens_C120)
Ana_sens_C120

save.image("C:/Users/XX.Rdata")