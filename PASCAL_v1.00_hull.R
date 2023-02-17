###################
#Seeding simulation
###################
ST <- Sys.time()

#Configuration settings
#----------------------
RUNID <- 30
set.seed(RUNID)
IDCode <- "60N_BasicRun_@67K"

#Variable limits
#---------------
N <- 1000000
MinDepth <- 1
MaxDepth <- 500
STL <- min(which(FA[1, 1:8760] > 0)) - (10 * 24) #10d prior to the pps onset
STU <- min(which(FA[1, 1:8760] > 0)) + (90 * 24) #30d after the pps onset
SDL <- 1
SDU <- 50
DVL <- 0
DVU <- roundup(max(VIrad) - 100, 100)
SDPL <- 0.4
SDPU <- 0.7
SAPL <- 0.7
SAPU <- 1.0
GAPL <- 0.5
GAPU <- 0.7
DVPMU <- 10
DVPML <- 2
STimePool <- seq(STL, STU, 24)
SDepthPool <- seq(SDL, SDU, 1)
SDPPool <- seq(SDPL, SDPU, 0.01)
SAPPool <- seq(SAPL, SAPU, 0.01)
GAPPool <- seq(GAPL, GAPU, 0.01)
DVPPool <- seq(DVL, DVU, 100)
ODPPool <- seq(50, MaxDepth, 50)
DVPMPool <- seq(DVPML, DVPMU, 1)

#Optimized parameters
#--------------------
STime <- sample(STimePool, N, replace = TRUE)
SDepth <- sample(SDepthPool, N, replace = TRUE)
DVP <- sample(DVPPool, N, replace = TRUE)
SDP <- sample(SDPPool, N, replace = TRUE)
SAP <- sample(SAPPool, N, replace = TRUE)
GAP <- sample(GAPPool, N, replace = TRUE)
ODP <- sample(ODPPool, N, replace = TRUE)
DVPM <- sample(DVPMPool, N, replace = TRUE)

#Factors
#-------
EDF <- 717 #Egg Development Factor (modified coefficient)
RAF <- 0.01
RAFScalar <- 15 #Scaling factor for risk adjustment
MRF <- 0.10 #0.0716277660995 #Max. 365d overwintering capability
WCF <- 2.50 #Carbon to total weight conversion assuming 40% C composition
GEF <- 0.903 #Growth factor assuming 40% C composition
FCF <- 30 #Food conversion factor (Chl-a to C)
GLF <- 0 #Growth limitation factor
SCF <- 1.0 #Swimming cost adjustment factor
MFW <- WCrc[13, 3] #Maximum adult female mass (ugC)
UEM <- 0.55 #Unit egg mass (ugC)
THZ <- 15000 #Termination horizon, ping 15000
RAModel <- lm(c(RAF / RAFScalar, RAF)~c(WCrc[1, 3], MFW))
OFfsetTreshold <- 720 #Max.deviation from life cycle coherance allowed (hrs)
OFfsetScalar <- 0.01 #Penalty for failing coherance within the treshold
UKey <- "uci2xd9oqbasg8b84ft7be5zx1jthp"
AToken <- "abovwrofnhqc4ij7nxjb5ambsia7c6"

#I/O parameters
#--------------
IterationID <- 0
DIni <- rep(NA, N)
DTer <- rep(NA, N)
Fecundity <- rep(NA, N)
DStg <- rep(NA, N)
LCTer <- rep(NA, N)
IFitness <- rep(NA, N)
FEp <- rep(NA, N)
NBREAK <- 15000

#Byte code compilation
#---------------------
require(pushoverr)
require(compiler)
Simulator_c <- cmpfun(Simulator)

#Simulation
#----------
PB <- txtProgressBar(min = 0, max = N, style = 3)
BreakPoints <- c(1, seq(1000, N, by = 1000))

for(i in 1:N){
      IterationID <- i
      Simulator_c()
      
      if(any(i == BreakPoints)){
            Sys.sleep(0.1)
            setTxtProgressBar(PB, i)
      }else{
            #proceed
      }
}

#Fitness weighting
#-----------------
SOffsetWModel_1 <- lm(c(OFfsetScalar, (OFfsetScalar/OFfsetTreshold)) ~ c(1, OFfsetTreshold))

if(any(IFitness > 0) == TRUE){
      FPIndex <- which(IFitness > 0)
      FP_IFitness <- IFitness[FPIndex]
      FP_STime <- STime[FPIndex]
      FP_FEp <- FEp[FPIndex]
      
      SOffset <- FP_STime - (FP_FEp - 8760)
      
      for(i in 1:length(SOffset)){
            if(SOffset[i] >= 0){
                  FP_IFitness[i] <- FP_IFitness[i]
            }else{
                  if(SOffset[i] >= -OFfsetTreshold){
                        OFS <- abs(SOffset[i])
                        FWeight <- (OFS * SOffsetWModel_1$coefficients[2]) + SOffsetWModel_1$coefficients[1]
                        FP_IFitness[i] <- FP_IFitness[i] * FWeight
                  }else{
                        FP_IFitness[i] <- 0
                  }
            }
      }
      
      IFitness[FPIndex] <- FP_IFitness
      rm(FPIndex, FP_IFitness, FP_FEp, FP_STime, SOffset, OFS, FWeight, SOffsetWModel_1)
}else{
      #No fitness weighting
}

SOffsetWModel_2 <- lm(c(OFfsetScalar, (OFfsetScalar/OFfsetTreshold)) ~ c(1, OFfsetTreshold))

if(any(IFitness > 0) == TRUE){
      FPIndex <- which(IFitness > 0)
      FP_IFitness <- IFitness[FPIndex]
      FP_STime <- STime[FPIndex]
      FP_LCTer <- LCTer[FPIndex]
      
      SOffset <- FP_STime - (FP_LCTer - 8760)
      
      for(i in 1:length(SOffset)){
            if(SOffset[i] <= 0){
                  FP_IFitness[i] <- FP_IFitness[i]
            }else{
                  if(SOffset[i] <= OFfsetTreshold){
                        OFS <- abs(SOffset[i])
                        FWeight <- (OFS * SOffsetWModel_2$coefficients[2]) + SOffsetWModel_2$coefficients[1]
                        FP_IFitness[i] <- FP_IFitness[i] * FWeight
                  }else{
                        FP_IFitness[i] <- 0
                  }
            }
      }
      
      IFitness[FPIndex] <- FP_IFitness
      rm(FPIndex, FP_IFitness, FP_STime, SOffset, OFS, FWeight, SOffsetWModel_2)
      
}else{
      #No fitness weighting
}

#File assembly and writing
#-------------------------
ID <- 1:N
ParameterSummary <- cbind(ID, STime, SDepth, DVP, SDP, SAP, GAP, ODP, DVPM, DIni, DTer, DStg, LCTer, FEp, Fecundity, IFitness)
write.table(ParameterSummary, file = "SeedOutput.txt", sep = " ", row.names = FALSE, col.names = TRUE)
rm(ID, STime, SDepth, DVP, SDP, SAP, GAP, ODP, DVPM, DIni, DTer, DStg, LCTer, FEp, Fecundity, IFitness)
gc()

ET <- Sys.time()
print(ET - ST)

#Completion push message to Nexus-5
#----------------------------------
msg <- paste0("Seeding completed: ", IDCode)
pushover(message = msg, user = UKey, app = AToken, sound = "pushover")