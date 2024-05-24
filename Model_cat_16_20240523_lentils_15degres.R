
#REMARQUE GENERALE - il faudrait arriver à faire un seul code pour lentilles et céréales - en principe, faisable en créant une fonction


#May 2024 

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

    # N0_tot_L.RData
    # prev_Bernoulli_L.RData
    # Qty_Bc_t0_L.RData
    # N0_tot_Bc3_L.RData
    # N0_tot_Bc4_L.RData
    # N_ht_Bc3_L.RData
    # N_ht_Bc4_L.RData
    # N_cooling_time_Bc3_L.RData
    # N_24h_nofridge_Bc3_L.RData
    # N_cooling_time_Bc4_L.RData
    # N_24h_nofridge_Bc4_L.RData
    # N_24h_Bc3_L.RData
    # N_48h_Bc3_L.RData
    # N_72h_Bc3_L.RData
    # N_96h_Bc3_L.RData
    # N_120h_Bc3_L.RData
    # N_24h_Bc4_L.RData
    # N_48h_Bc4_L.RData
    # N_72h_Bc4_L.RData
    # N_96h_Bc4_L.RData
    # N_120h_Bc4_L.RData

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

setwd("C:/Users/SECALIM20221/Desktop/recherche/Projets/Thes_Pauline_Mombert/Publi/Major_revision_V2")
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
muopt_Bc3_L <- 2.05 # µopt for tubers and lentils - Combase ID: Car_45 ; pH 6.5 ; temperature 30; G= 0.346 µ = 2.05 #FRENCH - on pourrait supprimer cette information qui n'est pas utilisée par la suite
z_Bc3       <- 6.2
D90_Bc3     <- 101.9
Dht_Bc3 = D90_Bc3 * 10^((90-Tht)/z_Bc3)

# B. cereus group IV
Tmin_Bc4    <- 7.44
Topt_Bc4    <- 38.71
Tmax_Bc4    <- 47.44
muopt_Bc4_C <- 1.41 # µopt for cereals (Ellouze et al., 2021 - see table 3) 
muopt_Bc4_L <- 2.05 # µopt for tubers and lentils - Combase ID: Car_45 ; pH 6.5 ; temperature 30; G= 0.346 µ = 2.05 #FRENCH - on pourrait supprimer cette information qui n'est pas utilisée par la suite
z_Bc4       <- 7.128
D100_Bc4    <- 6.08
Dht_Bc4 = D100_Bc4 * 10^((100-Tht)/z_Bc4)


#--------------------------------------------------------------#


###########################
#   Heat treatment time   #
###########################

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

# Truncating the distribution so that temperatures do not fall below 0°C or exceed 15°C
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
Nb_portions_per_dish  <- 6 #FRENCH - il serait bien de justifier cette valeur
dish_size             <- Qty_cat16_day*Nb_portions_per_dish

#------------------------------- Results and heat inactivation treatment #FRENCH - pas sure que ce soit le meilleur titre mais je laisse pour le moment

#################################################################################################
##                                                                                             ##
##------------------- MODULE 1: Initial contamination before heat treatment -------------------##
##                                                                                             ##
#################################################################################################

## Lentil-based dishes => Data from Blakey and Priest, 1980 
N0          <- (8-6) +(4-2) + (3-1) + (3-2) + 3 # Number of uncontaminated samples
N_inf_2log  <- 1 + 1                            # Number of samples containing <= 2 log CFU.g-1 of B. cereus
N_sup_2log  <- (6-1) + 1 + 2 + 1 + 1 + 2 + 1    # Number of samples containing 2 to 3 log CFU.g-1 of B. cereus
Cmax        <- 45000 # CFU.g-1
Ntot        <- N0 + N_inf_2log + N_sup_2log

# Generating probabilities of B. cereus 
  # Dishes contaminated with less than 2 log CFU.g-1
N0_Bc_inf_2log  <- runif(100, min=1, max=100*dish_size) #quantity of B. cereus per dish
pct_Bc_inf_2log <- (N_inf_2log/(N_inf_2log + N_sup_2log))/100
P_Bc_inf_2log   <- rep(pct_Bc_inf_2log, 100)

  # Dishes contaminated with 2 to 3 log CFU.g-1
N0_Bc_sup_2log  <- runif(100, min=100*dish_size, max=Cmax*dish_size) #quantity of B. cereus per dish
pct_Bc_sup_2log <- (N_sup_2log/(N_inf_2log + N_sup_2log))/100
P_Bc_sup_2log   <- rep(pct_Bc_sup_2log, 100)

N0_tot_L          <- mcstoc(rempiricalD, type="V", values=c(N0_Bc_inf_2log,N0_Bc_sup_2log), prob=c(P_Bc_inf_2log,P_Bc_sup_2log)) #quantity of B. cereus (all groups combined)

# Calculating prevalence and initial quantities of group III & group IV B. cereus before heat treatment
## Prevalence and initial quantity of B. cereus
prev_t0_L         <- mcstoc(rbeta, type="U", shape1=N_inf_2log + N_sup_2log+1, shape2=Ntot-(N_inf_2log + N_sup_2log)+1, ncp=0) #prevalence 


##################################################################################################
##                                                                                              ##
##------------------- MODULE 2: Reduction of contamination by heat treatment -------------------##
##                                                                                              ##
##################################################################################################

# Survival of bacterial spores from group III B. cereus after heat treatment
p_surv_Bc3_L  <- 10^(-Tht/Dht_Bc3) #survival probability of bacterial spores from group III B. cereus

# Survival of bacterial spores from group IV B. cereus after heat treatment
p_surv_Bc4_L  <- 10^(-Tht/Dht_Bc4) #survival probability of bacterial spores from group IV B. cereus


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
gamma_cooling_Bc3_L   <- ifelse(Tcooling<Tmin_Bc3, 0, ifelse(Tcooling>Tmax_Bc3, 0, ((Tcooling-Tmin_Bc3)^2*(Tcooling-Tmax_Bc3)/((Topt_Bc3-Tmin_Bc3)*((Topt_Bc3-Tmin_Bc3)*(Tcooling-Topt_Bc3)-(Topt_Bc3-Tmax_Bc3)*(Topt_Bc3+Tmin_Bc3-2*Tcooling))))))

# Calculating gamma and quantity of bacteria after cooling time for group IV B. cereus
gamma_cooling_Bc4_L   <- ifelse(Tcooling<Tmin_Bc4, 0, ifelse(Tcooling>Tmax_Bc4, 0, ((Tcooling-Tmin_Bc4)^2*(Tcooling-Tmax_Bc4)/((Topt_Bc4-Tmin_Bc4)*((Topt_Bc4-Tmin_Bc4)*(Tcooling-Topt_Bc4)-(Topt_Bc4-Tmax_Bc4)*(Topt_Bc4+Tmin_Bc4-2*Tcooling))))))


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
storage_lentil      <- subset(data_storage_leftovers, dish_type=="lentil")

# FRENCH : AJOUTER DANS TABLEUR EXCEL 1 NO LEFT OVER
#No leftovers
nb_noleftovers_cons_L<-subset(storage_lentil, storage_in_days=="0")$number_indiv
# 1 day (i.e 24h) (max)
nb_24h_cons_L       <- subset(storage_lentil, storage_in_days=="1")$number_indiv
# 2 days (i.e. 48h) (max)
nb_48h_cons_L       <- subset(storage_lentil, storage_in_days=="2")$number_indiv 
# 3 days (i.e. 72h) (max)
nb_72h_cons_L       <- subset(storage_lentil, storage_in_days=="3")$number_indiv
# 4 days (i.e. 96h) (max) - option "more than 3 days" considered to be 4 or 5 days, with the same probability
nb_96h_cons_L       <- subset(storage_lentil, storage_in_days==">3")$number_indiv/2  
# 5 days (i.e. 120h) (max) - option "more than 3 days" considered to be 4 or 5 days, with the same probability
nb_120h_cons_L      <- subset(storage_lentil, storage_in_days==">3")$number_indiv/2

Sum_nb_cons_L       <- nb_noleftovers_cons_L + nb_24h_cons_L + nb_48h_cons_L + nb_72h_cons_L + nb_96h_cons_L + nb_120h_cons_L


prev_no_leftovers 	<- mcstoc(rbeta, type="U", shape1=nb_noleftovers_cons_L + 1, shape2=Sum_nb_cons_L-nb_noleftovers_cons_L + 1, ncp=0) 
prev_24h_cons_L     <- mcstoc(rbeta, type="U", shape1=nb_24h_cons_L+1, shape2=Sum_nb_cons_L-nb_24h_cons_L+1, ncp=0)  
prev_48h_cons_L     <- mcstoc(rbeta, type="U", shape1=nb_48h_cons_L+1, shape2=Sum_nb_cons_L-nb_48h_cons_L+1, ncp=0)  
prev_72h_cons_L     <- mcstoc(rbeta, type="U", shape1=nb_72h_cons_L+1, shape2=Sum_nb_cons_L-nb_72h_cons_L+1, ncp=0)  
prev_96h_cons_L     <- mcstoc(rbeta, type="U", shape1=nb_96h_cons_L+1, shape2=Sum_nb_cons_L-nb_96h_cons_L+1, ncp=0)
prev_120h_cons_L    <- mcstoc(rbeta, type="U", shape1=nb_120h_cons_L+1, shape2=Sum_nb_cons_L-nb_120h_cons_L+1, ncp=0)


Resultat_L<- mcmodel({


##############----------------------------------------------------

#################################################################################################
##                                                                                             ##
##------------------- MODULE 1: Initial contamination before heat treatment -------------------##
##                                                                                             ##
#################################################################################################



# Calculating prevalence and initial quantities of group III & group IV B. cereus before heat treatment
## Prevalence and initial quantity of B. cereus

N0_tot_L          <- mcstoc(rempiricalD, type="V", values=c(N0_Bc_inf_2log,N0_Bc_sup_2log), prob=c(P_Bc_inf_2log,P_Bc_sup_2log)) #quantity of B. cereus (all groups combined)
prev_Bernoulli_L  <- mcstoc(rbern, type="VU", prev_t0_L)

Qty_Bc_t0_L       <- N0_tot_L*prev_Bernoulli_L 


# Generating proportions of group III and group IV B. cereus
Pct_Bc3_L <- mcstoc(runif, type="U", min=0.2, max=0.8)
Pct_Bc4_L <- 1-Pct_Bc3_L

# Calculating corresponding amounts of bacteria
N0_tot_Bc3_L <- Qty_Bc_t0_L*Pct_Bc3_L
N0_tot_Bc4_L <- Qty_Bc_t0_L*Pct_Bc4_L 


##################################################################################################
##                                                                                              ##
##------------------- MODULE 2: Reduction of contamination by heat treatment -------------------##
##                                                                                              ##
##################################################################################################

# Survival of bacterial spores from group III B. cereus after heat treatment
N_surv_Bc3_L  <- N0_tot_Bc3_L*p_surv_Bc3_L
N_ht_Bc3_L    <- mcstoc(rpois, type="VU", lambda=N_surv_Bc3_L)

# Survival of bacterial spores from group IV B. cereus after heat treatment
N_surv_Bc4_L  <- N0_tot_Bc4_L*p_surv_Bc4_L
N_ht_Bc4_L    <- mcstoc(rpois, type="VU", lambda=N_surv_Bc4_L)


#################################################################################################
##                                                                                             ##
##------------------- MODULE 3: Recontamination during cooling and storage  -------------------##
##                                                                                             ##
#################################################################################################

N_cooling_time_Bc3_L  <- ifelse(Cooling_time<lag_cooling, N_ht_Bc3_L, N_ht_Bc3_L*exp(muopt_Bc3_L*gamma_cooling_Bc3_L*(Cooling_time-lag_cooling))) #Quantity of group III B. cereus after cooling
N_cooling_time_Bc4_L  <- ifelse(Cooling_time<lag_cooling, N_ht_Bc4_L, N_ht_Bc4_L*exp(muopt_Bc4_L*gamma_cooling_Bc4_L*(Cooling_time-lag_cooling))) #Quantity of group IV B. cereus after cooling 


#--- Growth during storage at ambient temperature (i.e. no fridge)
###------ 24h growth for group III B. cereus
N_24h_nofridge_Bc3_L <- N_ht_Bc3_L*exp(muopt_Bc3_L*gamma_cooling_Bc3_L*(24-lag_cooling))

###------ 24h growth for group IV B. cereus
N_24h_nofridge_Bc4_L <- N_ht_Bc4_L*exp(muopt_Bc4_L*gamma_cooling_Bc4_L*(24-lag_cooling))


# Lag time during storage in the fridge
lag_storage <- mcstoc(runif, type="U", min=lag_min_fridge, max=lag_max_fridge)

#---- Growth during storage in the fridge
T_lagtime_gt1week <- 10 # Lag time was assumed to be more than one week if the refrigerator temperature was low enough. 
##----- Group III B. cereus
###------ 24h growth
N_24h_Bc3_L   <- ifelse(T_frigde<T_lagtime_gt1week, N_cooling_time_Bc3_L, ifelse(24<lag_storage, N_cooling_time_Bc3_L, N_cooling_time_Bc3_L*exp(muopt_Bc3_L*gamma_fridge_Bc3*(24-lag_storage)))) #FRENCH - j'ai changé la formule en supprimant le *1/6 pour raisonner à l'échelle du plat jusqu'au bout pour tout le monde
###------ 48h growth
N_48h_Bc3_L   <- ifelse(T_frigde<T_lagtime_gt1week, N_cooling_time_Bc3_L, ifelse(48<lag_storage, N_cooling_time_Bc3_L, N_cooling_time_Bc3_L*exp(muopt_Bc3_L*gamma_fridge_Bc3*(48-lag_storage)))) #FRENCH - j'ai changé la formule en supprimant le *1/6 pour raisonner à l'échelle du plat jusqu'au bout pour tout le monde
###------ 72h growth
N_72h_Bc3_L   <- ifelse(T_frigde<T_lagtime_gt1week, N_cooling_time_Bc3_L, ifelse(72<lag_storage, N_cooling_time_Bc3_L, N_cooling_time_Bc3_L*exp(muopt_Bc3_L*gamma_fridge_Bc3*(72-lag_storage)))) #FRENCH - j'ai changé la formule en supprimant le *1/6 pour raisonner à l'échelle du plat jusqu'au bout pour tout le monde
###------ 96h growth
N_96h_Bc3_L   <- ifelse(T_frigde<T_lagtime_gt1week, N_cooling_time_Bc3_L, ifelse(96<lag_storage, N_cooling_time_Bc3_L, N_cooling_time_Bc3_L*exp(muopt_Bc3_L*gamma_fridge_Bc3*(96-lag_storage)))) #FRENCH - j'ai changé la formule en supprimant le *1/6 pour raisonner à l'échelle du plat jusqu'au bout pour tout le monde
###------ 120h growth
N_120h_Bc3_L  <- ifelse(T_frigde<T_lagtime_gt1week, N_cooling_time_Bc3_L, ifelse(120<lag_storage, N_cooling_time_Bc3_L, N_cooling_time_Bc3_L*exp(muopt_Bc3_L*gamma_fridge_Bc3*(120-lag_storage)))) #FRENCH - j'ai changé la formule en supprimant le *1/6 pour raisonner à l'échelle du plat jusqu'au bout pour tout le monde


##----- Group IV B. cereus
###------ 24h growth
N_24h_Bc4_L   <- ifelse(T_frigde<T_lagtime_gt1week, N_cooling_time_Bc4_L, ifelse(24<lag_storage, N_cooling_time_Bc4_L, N_cooling_time_Bc4_L*exp(muopt_Bc4_L*gamma_fridge_Bc4*(24-lag_storage)))) #FRENCH - j'ai changé la formule en supprimant le *1/6 pour raisonner à l'échelle du plat jusqu'au bout pour tout le monde
###------ 48h growth
N_48h_Bc4_L   <- ifelse(T_frigde<T_lagtime_gt1week, N_cooling_time_Bc4_L, ifelse(48<lag_storage, N_cooling_time_Bc4_L, N_cooling_time_Bc4_L*exp(muopt_Bc4_L*gamma_fridge_Bc4*(48-lag_storage)))) #FRENCH - j'ai changé la formule en supprimant le *1/6 pour raisonner à l'échelle du plat jusqu'au bout pour tout le monde
###------ 72h growth
N_72h_Bc4_L   <- ifelse(T_frigde<T_lagtime_gt1week, N_cooling_time_Bc4_L, ifelse(72<lag_storage, N_cooling_time_Bc4_L, N_cooling_time_Bc4_L*exp(muopt_Bc4_L*gamma_fridge_Bc4*(72-lag_storage)))) #FRENCH - j'ai changé la formule en supprimant le *1/6 pour raisonner à l'échelle du plat jusqu'au bout pour tout le monde
###------ 96h growth
N_96h_Bc4_L   <- ifelse(T_frigde<T_lagtime_gt1week, N_cooling_time_Bc4_L, ifelse(96<lag_storage, N_cooling_time_Bc4_L, N_cooling_time_Bc4_L*exp(muopt_Bc4_L*gamma_fridge_Bc4*(96-lag_storage)))) #FRENCH - j'ai changé la formule en supprimant le *1/6 pour raisonner à l'échelle du plat jusqu'au bout pour tout le monde
###------ 120h growth
N_120h_Bc4_L  <- ifelse(T_frigde<T_lagtime_gt1week, N_cooling_time_Bc4_L, ifelse(120<lag_storage, N_cooling_time_Bc4_L, N_cooling_time_Bc4_L*exp(muopt_Bc4_L*gamma_fridge_Bc4*(120-lag_storage)))) #FRENCH - j'ai changé la formule en supprimant le *1/6 pour raisonner à l'échelle du plat jusqu'au bout pour tout le monde


#####Addition des (N_24h_nofridge_Bc3_L + N_24h_nofridge_Bc4_L)/Nb_portions_per_dish/mean_portion_size



#####Addition des (N_24h_nofridge_Bc3_L + N_24h_nofridge_Bc4_L)/Nb_portions_per_dish/mean_portion_size
C0_tot_Bc_L<-(N0_tot_Bc3_L+N0_tot_Bc4_L)/Nb_portions_per_dish/mean_portion_size
C_ht_Bc_L<-(N_ht_Bc3_L+N_ht_Bc4_L)/Nb_portions_per_dish/mean_portion_size
C_cooling_time_Bc_L<-(N_cooling_time_Bc3_L + N_cooling_time_Bc4_L)/Nb_portions_per_dish/mean_portion_size
C_24h_Bc_L<-(N_24h_Bc3_L+N_24h_Bc4_L)/Nb_portions_per_dish/mean_portion_size
C_48h_Bc_L<-(N_48h_Bc3_L+N_48h_Bc4_L)/Nb_portions_per_dish/mean_portion_size
C_72h_Bc_L<-(N_72h_Bc3_L+N_72h_Bc4_L)/Nb_portions_per_dish/mean_portion_size
C_96h_Bc_L<-(N_96h_Bc3_L+N_96h_Bc4_L)/Nb_portions_per_dish/mean_portion_size
C_120h_Bc_L<-(N_120h_Bc3_L+N_120h_Bc4_L)/Nb_portions_per_dish/mean_portion_size
C_24h_nofridge_Bc_L<-(N_24h_nofridge_Bc3_L+N_24h_nofridge_Bc4_L)/Nb_portions_per_dish/mean_portion_size


#For figure XXX
concentration_L_5days        <- C_120h_Bc_L
log_concentration_0_L_5days  <- ifelse(concentration_L_5days<=0.0001, 0, log10(concentration_L_5days))
log_concentration_08_L_5days <- ifelse(log_concentration_0_L_5days >8, 8, log_concentration_0_L_5days)


#---- Calculating risks
risk_hot_L          <- C_ht_Bc_L > cutoff #FRENCH - j'ai changé la formule pour récupérer la quantité de B. cereus par portion (et pas par plat)

risk_24h_L          <- C_24h_Bc_L > cutoff #FRENCH - j'ai changé la formule en divisant par le nbre de portions pour récupérer la quantité de B. cereus par portion (et pas par plat)
risk_48h_L          <- C_48h_Bc_L > cutoff #FRENCH - j'ai changé la formule en divisant par le nbre de portions pour récupérer la quantité de B. cereus par portion (et pas par plat)
risk_72h_L          <- C_72h_Bc_L	> cutoff #FRENCH - j'ai changé la formule en divisant par le nbre de portions pour récupérer la quantité de B. cereus par portion (et pas par plat)
risk_96h_L          <- C_96h_Bc_L > cutoff #FRENCH - j'ai changé la formule en divisant par le nbre de portions pour récupérer la quantité de B. cereus par portion (et pas par plat)
risk_120h_L         <- C_120h_Bc_L > cutoff #FRENCH - j'ai changé la formule en divisant par le nbre de portions pour récupérer la quantité de B. cereus par portion (et pas par plat)

risk_24h_nofridge_L <- C_24h_nofridge_Bc_L > cutoff #FRENCH - j'ai changé la formule pour récupérer la quantité de B. cereus par portion (et pas par plat)

risk_portion_fridge_L     <- (1-Pct_wo_fridge)*(prev_24h_cons_L*risk_24h_L + prev_48h_cons_L*risk_48h_L + prev_72h_cons_L*risk_72h_L + prev_96h_cons_L*risk_96h_L + prev_120h_cons_L*risk_120h_L)


mc(prev_Bernoulli_L, T_frigde, Cooling_time, Qty_cat16_day,
  N0_tot_L, C0_tot_Bc_L, C_ht_Bc_L, C_cooling_time_Bc_L, C_24h_Bc_L, C_48h_Bc_L, C_72h_Bc_L, C_96h_Bc_L, C_120h_Bc_L, C_24h_nofridge_Bc_L,
  risk_hot_L, risk_24h_L, risk_48h_L, risk_72h_L, risk_96h_L, risk_120h_L,risk_24h_nofridge_L,risk_portion_fridge_L,
  log_concentration_08_L_5days) #lag_storage, lag_cooling




#mc(prev_Bernoulli_L, T_frigde, Cooling_time, Qty_cat16_day,
#  N0_tot_L,
#  N0_tot_Bc_L, N0_tot_Bc3_L, N0_tot_Bc4_L, 
#  N_ht_Bc_L, N_ht_Bc3_L, N_ht_Bc4_L, 
#  N_cooling_time_Bc_L, N_cooling_time_Bc3_L, N_cooling_time_Bc4_L, 
#  N_24h_Bc_L, N_24h_Bc3_L, N_24h_Bc4_L, 
#  N_48h_Bc_L, N_48h_Bc3_L, N_48h_Bc4_L,
#  N_72h_Bc_L,N_72h_Bc3_L, N_72h_Bc4_L,
#  N_96h_Bc_L, N_96h_Bc3_L, N_96h_Bc4_L,
#  N_120h_Bc_L, N_120h_Bc3_L, N_120h_Bc4_L,
#  N_24h_nofridge_Bc_L, N_24h_nofridge_Bc3_L, N_24h_nofridge_Bc4_L,
#  risk_hot_L, risk_24h_L, risk_48h_L, risk_72h_L, risk_96h_L, risk_120h_L, risk_portion_no_fridge_L, risk_portion_hot_L, risk_portion_fridge_L,risk_24h_nofridge_L,
#  log_concentration_08_L_5days) #lag_storage, lag_cooling

})

rm(Qty_cat16,Qty_cat16_max)

Resu_L        <- evalmcmod(Resultat_L, nsv=ndvar(), nsu=ndunc())
#Resu.final_L  <- summary(Resu_L, probs=c(0.5, 0.95, 0.99), lim=c(0.025,0.975))

#------------------------------- Sensitivity analysis (only for variability) -------------------------------#

summary(Resu_L, probs=c(0.5, 0.95, 0.99), lim=c(0.025,0.975))


Resu.final_L$N0_tot_L
Resu.final_L$N_ht_Bc_L
Resu.final_L$N_cooling_time_Bc_L
Resu.final_L$N_24h_Bc_L
Resu.final_L$N_48h_Bc_L
Resu.final_L$N_72h_Bc_L
Resu.final_L$N_96h_Bc_L
Resu.final_L$N_120h_Bc_L
Resu.final_L$N_24h_nofridge_Bc_L

Resu.final_L$risk_hot_L
Resu.final_L$risk_24h_L
Resu.final_L$risk_48h_L
Resu.final_L$risk_72h_L
Resu.final_L$risk_96h_L
Resu.final_L$risk_120h_L
Resu.final_L$risk_24h_nofridge_L
Resu.final_L$risk_portion_no_fridge_L
Resu.final_L$risk_portion_hot_L
Resu.final_L$risk_portion_fridge_L


# Plotting results in log CFU.g-1

plot(Resu_L$log_concentration_08_L_5days, ylab="Cumulative probability density - 5-day chilled storage", xlab ="log CFU.g-1", ymin=0, xlim=c(0,6)) #main="Bacterial concentration in cereal-based dishes (CFU.g-1)", adj = 0.5, line = 6
abline(v=5, col="red", lty=2, lwd=3)

Resu.final_L$log_concentration_08_L_5days


plot(Resu_L$N_120h_Bc_L)

Ana_sens_L120 <- tornado(Resu_L, output="log_concentration_08_L_5days", use = "pairwise.complete.obs", method = c("spearman"), lim=c(0.025, 0.975))
plot(Ana_sens_L120)
Ana_sens_L120

save.image("C:/Users/SECALIM20221/Desktop/recherche/Projets/Thes_Pauline_Mombert/Publi/Major_revision_V2/Results/Model lentils full 200000 20240523.RData")
