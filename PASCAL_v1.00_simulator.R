Simulator <- function(){
      
      #Object retrival - individual-specific parameters & other things
      IND <- IterationID                  
      StartPing <- STime[IND]             
      StartDepth <- SDepth[IND]           
      DV <- DVP[IND]                      
      DVM <- DVPM[IND]                    
      SD <- SDP[IND]
      SA <- SAP[IND]
      GA <- GAP[IND]
      OD <- ODP[IND]
      
      RAMSlo <- RAModel$coefficients[2]
      RAMInt <- RAModel$coefficients[1]
      LifetimeMaxCW <- WCrc[1, 3]
      LifeState <- as.factor("ALIVE")
      DiapauseStage <- NA
      CessasionStage <- NA
      DiapauseStartTime <- NA
      DiapauseEndTime <- NA
      DiapauseStartSW <- NA
      FEPState <- 0
      FEPTime <- NA
      
      #ParOutput <- NULL
      
      ###############
      #Egg-Nauplius I
      ###############
      S <- 0
      DTTrack <- NULL
      CurrentCW <- .subset2(WCrc, 3)[(S + 1)]
      CurrentSW <- 0
      CurrentRW <- 0
      Suvivorship <- 1
      CurrentEggProduction <- 0
      CurrentFitness <- 0
      CurrentPing <- StartPing - 1
      CurrentDepth <- StartDepth
      
      repeat{
            CurrentPing <- CurrentPing + 1
            CurrentTemp <- .subset2(Temperature, CurrentPing)[CurrentDepth]
            CurrentPProb <- .subset2(PProb, CurrentPing)[CurrentDepth]
            CurrentDT <- round(((EDF * (CurrentTemp + 9.11)^-2.05) * 24), 0)
            CurrentPRisk <- CurrentPProb * ((CurrentCW * RAMSlo) + RAMInt)
            Suvivorship <- Suvivorship * (1 - CurrentPRisk)
            DTTrack <- append(DTTrack, CurrentDT)
            ElapsedPing <- CurrentPing - StartPing
            
            if(Suvivorship < 0.0001){
                  LifeState <- as.factor("DEAD")
                  CessasionStage <- S
                  break
            }else if(ElapsedPing >= mean(DTTrack)){
                  S <- 1
                  break
            }else{
                  #Proceed iteration
            }
      }
      
      #############################
      #DVM Pre-processing directive
      #############################
      
      #DVPM Adjust - one time only
      DVModel_PreC4 <- lm(c(DV, DV / DVM) ~ c(WCrc[1, 3], MFW))
      DVModel_PostC4 <- lm(c(DV, DV / DVM) ~ c(WCrc[1, 3], MFW))
      
      #############################
      #Non-feeding nauplii (NI-NII)
      #############################
      
      if(LifeState == "DEAD"){
            #Do nothing and proceed
      }else{
            repeat{
                  CurrentPing <- CurrentPing + 1
                  CurrentTWMG <- (CurrentCW * WCF) / 1000
                  MaxDistance <- round((3.7036 * CurrentCW^0.5853), 0)
                  MigrationCeiling <- CurrentDepth - MaxDistance
                  MigrationFloor <- CurrentDepth + MaxDistance
                  
                  MigrationCeiling <- ifelse(MigrationCeiling < MinDepth, yes = MinDepth, no = MigrationCeiling)
                  MigrationFloor <- ifelse(MigrationFloor > MaxDepth, yes = MaxDepth, no = MigrationFloor)
                  SearchRange <- MigrationCeiling:MigrationFloor
                  
                  DVAdj <- (CurrentCW * DVModel_PreC4$coefficients[2]) + DVModel_PreC4$coefficients[1]
                  CurrentVI <- .subset2(VIrad, CurrentPing)[CurrentDepth]
                  
                  if(CurrentVI <= DVAdj){
                        #1.Growth potential estimation
                        TempRange <- .subset2(Temperature, CurrentPing)[SearchRange]
                        GmaxRange <- GEF * exp(0.110 * TempRange)
                        GPotentialRange <- GmaxRange * CurrentTWMG
                        
                        #2.Decision making
                        GloTempRange <- .subset2(Temperature, CurrentPing)[1:100]               #Values of upper 100m only used
                        GloGmaxRange <- GEF * exp(0.110 * GloTempRange)
                        GloGPotentialRange <- GloGmaxRange * CurrentTWMG
                        GPotentialMaxRange <- which(GPotentialRange == max(GPotentialRange))    #Set of index locations towards MinDepth -> MaxDepth
                        
                        if(max(GPotentialRange) < max(GloGPotentialRange)){
                              GPotentialMaxPos <- min(GPotentialMaxRange)                       #Single index location
                              FutureDepth <- SearchRange[GPotentialMaxPos]                     #Depth corresponding to the single index location above
                        }else{
                              GPotentialMaxPos <- max(GPotentialMaxRange)                       #Single index location
                              FutureDepth <- SearchRange[GPotentialMaxPos]                      #Depth corresponding to the single index location above
                        }
                        
                        #3.Final estimations
                        MovingTime <- abs(FutureDepth - CurrentDepth) / MaxDistance
                        k <- 0.375 * exp(0.0546 * TempRange[GPotentialMaxPos])
                        m <- 0.858 * exp(-0.008 * TempRange[GPotentialMaxPos])
                        CurrentSWCost <- k * CurrentTWMG^m * MovingTime * SCF
                        
                        NetGPotential <- max(GPotentialRange) - CurrentSWCost
                        CurrentGPotential <- CurrentCW + NetGPotential
                        
                        if(CurrentGPotential >= LifetimeMaxCW){
                              LifetimeMaxCW <- CurrentGPotential
                              CurrentSRisk <- 0
                        }else{
                              CWLTreshold <- LifetimeMaxCW / 2
                              
                              if(CurrentGPotential <= CWLTreshold){
                                    CurrentSRisk <- 1
                                    CurrentGPotential <- CWLTreshold
                              }else{
                                    CWLoss <- (50 / (LifetimeMaxCW - CWLTreshold)) * (LifetimeMaxCW - CurrentGPotential)
                                    CurrentSRisk <- 0.02 * CWLoss
                              }
                        }
                        
                        CurrentCW <- CurrentGPotential
                        CurrentDepth <- FutureDepth
                        CurrentPProb <- .subset2(PProb, CurrentPing)[CurrentDepth]
                        CurrentPRisk <- CurrentPProb * ((CurrentCW * RAMSlo) + RAMInt)
                        Suvivorship <- Suvivorship * (1 - (CurrentPRisk + CurrentSRisk))
                        
                  }else{
                        #1.Decision making
                        VIRadRange <- .subset2(VIrad, CurrentPing)[SearchRange]
                        PRiskMinRange <-  which(VIRadRange <= DVAdj)
                        
                        if(length(PRiskMinRange) == 0){
                              PRiskMinPos <- which.max(SearchRange)
                              FutureDepth <- SearchRange[PRiskMinPos]
                        }else if(length(PRiskMinRange) == 1){
                              PRiskMinPos <- PRiskMinRange
                              FutureDepth <- SearchRange[PRiskMinPos]
                        }else{
                              PRiskMinPos <- PRiskMinRange[1]
                              FutureDepth <- SearchRange[PRiskMinPos]
                        }
                        
                        #2.Estimations
                        CurrentTemp <- .subset2(Temperature, CurrentPing)[FutureDepth]
                        CurrentGmax <- GEF * exp(0.110 * CurrentTemp)
                        
                        MovingTime <- abs(FutureDepth - CurrentDepth) / MaxDistance
                        k <- 0.375 * exp(0.0546 * CurrentTemp)
                        m <- 0.858 * exp(-0.008 * CurrentTemp)
                        CurrentSWCost <- k * CurrentTWMG^m * MovingTime * SCF
                        
                        NetGPotential <- (CurrentGmax * CurrentTWMG) - CurrentSWCost
                        CurrentGPotential <- CurrentCW + NetGPotential
                        
                        if(CurrentGPotential >= LifetimeMaxCW){
                              LifetimeMaxCW <- CurrentGPotential
                              CurrentSRisk <- 0
                        }else{
                              CWLTreshold <- LifetimeMaxCW / 2
                              
                              if(CurrentGPotential <= CWLTreshold){
                                    CurrentSRisk <- 1
                                    CurrentGPotential <- CWLTreshold
                              }else{
                                    CWLoss <- (50 / (LifetimeMaxCW - CWLTreshold)) * (LifetimeMaxCW - CurrentGPotential)
                                    CurrentSRisk <- 0.02 * CWLoss
                              }
                        }
                        
                        CurrentCW <- CurrentGPotential
                        CurrentDepth <- FutureDepth
                        CurrentPProb <- .subset2(PProb, CurrentPing)[CurrentDepth]
                        CurrentPRisk <- CurrentPProb * ((CurrentCW * RAMSlo) + RAMInt)
                        Suvivorship <- Suvivorship * (1 - (CurrentPRisk + CurrentSRisk))
                        
                  }
                  
                  if(Suvivorship < 0.0001){
                        LifeState <- as.factor("DEAD")
                        CessasionStage <- S
                        break
                  }else if(CurrentCW >= WCrc[3, 3] & CurrentCW < WCrc[4, 3]){
                        S <- 2
                        #Proceed iteration
                  }else if(CurrentCW >= WCrc[4, 3]){
                        S <- S + 1
                        break
                  }else{
                        #Proceed iteration
                  }
            }
      }
      
      ##########################################
      #Feeding stages without a store (NII-CIII)
      ##########################################
      
      if(LifeState == "DEAD"){
            #Do nothing and proceed
      }else{
            repeat{
                  CurrentPing <- CurrentPing + 1
                  CurrentTWMG <- (CurrentCW * WCF) / 1000
                  MaxDistance <- round((3.7036 * CurrentCW^0.5853), 0)
                  MigrationCeiling <- CurrentDepth - MaxDistance
                  MigrationFloor <- CurrentDepth + MaxDistance
                  
                  MigrationCeiling <- ifelse(MigrationCeiling < MinDepth, yes = MinDepth, no = MigrationCeiling)
                  MigrationFloor <- ifelse(MigrationFloor > MaxDepth, yes = MaxDepth, no = MigrationFloor)
                  SearchRange <- MigrationCeiling:MigrationFloor
                  
                  DVAdj <- (CurrentCW * DVModel_PreC4$coefficients[2]) + DVModel_PreC4$coefficients[1]
                  CurrentVI <- .subset2(VIrad, CurrentPing)[CurrentDepth]
                  
                  if(CurrentVI <= DVAdj){
                        #1.Growth potential estimation
                        TempRange <- .subset2(Temperature, CurrentPing)[SearchRange]
                        FARange <- .subset2(FA, CurrentPing)[SearchRange]
                        FACRange <- FARange * (FCF / 1000)
                        GmaxRange <- GEF * exp(0.110 * TempRange)
                        a <- 0.7
                        b <- 1.777 * exp(0.234 * TempRange)
                        n <- 0.681 * exp(0.0199 * TempRange)
                        k <- 0.375 * exp(0.0546 * TempRange)
                        m <- 0.858 * exp(-0.008 * TempRange)
                        FA_CRange <- (k * CurrentTWMG^(m - n)) / (a * b)
                        FA_SRange <- FA_CRange + (GmaxRange / (a * b * CurrentTWMG^(n - 1)))
                        GPotentialRange <- rep(NA, length(FACRange))
                        
                        for(w in 1:length(FACRange)){
                              if(FACRange[w] >= FA_SRange[w]){
                                    GPotentialRange[w] <- GmaxRange[w] * CurrentTWMG
                              }else{
                                    GPotentialRange[w] <- (a * b[w] * CurrentTWMG^n[w] * FACRange[w]) - (k[w] * CurrentTWMG^m[w])
                              }
                        }
                        
                        #2.Decision making
                        GloTempRange <- .subset2(Temperature, CurrentPing)[1:100]
                        GloGmaxRange <- GEF * exp(0.110 * GloTempRange)
                        GloFARange <- .subset2(FA, CurrentPing)[1:100]
                        GloFACRange <- GloFARange * (FCF / 1000)
                        a <- 0.7
                        Glo_b <- 1.777 * exp(0.234 * GloTempRange)
                        Glo_n <- 0.681 * exp(0.0199 * GloTempRange)
                        Glo_k <- 0.375 * exp(0.0546 * GloTempRange)
                        Glo_m <- 0.858 * exp(-0.008 * GloTempRange)
                        GloFA_CRange <- (Glo_k * CurrentTWMG^(Glo_m - Glo_n)) / (a * Glo_b)
                        GloFA_SRange <- GloFA_CRange + (GloGmaxRange / (a * Glo_b * CurrentTWMG^(Glo_n - 1)))
                        
                        for(x in 1:length(GloFACRange)){
                              if(GloFACRange[x] >= GloFA_SRange[x]){
                                    GloGPotentialRange[x] <- GloGmaxRange[x] * CurrentTWMG
                              }else{
                                    GloGPotentialRange[x] <- (a * Glo_b[x] * CurrentTWMG^Glo_n[x] * GloFACRange[x]) - (Glo_k[x] * CurrentTWMG^Glo_m[x])
                              }
                        }
                        
                        GPotentialMaxRange <- which(GPotentialRange == max(GPotentialRange))
                        
                        if(max(GPotentialRange) < max(GloGPotentialRange)){
                              GPotentialMaxPos <- min(GPotentialMaxRange)                       #Single index location
                              FutureDepth <- SearchRange[GPotentialMaxPos]                     #Depth corresponding to the single index location above
                        }else{
                              GPotentialMaxPos <- max(GPotentialMaxRange)                       #Single index location
                              FutureDepth <- SearchRange[GPotentialMaxPos]                      #Depth corresponding to the single index location above
                        }
                        
                        #3.Estimation
                        MovingTime <- abs(FutureDepth - CurrentDepth) / MaxDistance
                        k <- 0.375 * exp(0.0546 * TempRange[GPotentialMaxPos])
                        m <- 0.858 * exp(-0.008 * TempRange[GPotentialMaxPos])
                        CurrentSWCost <- k * CurrentTWMG^m * MovingTime * SCF
                        
                        GrossGPotential <- max(GPotentialRange)
                        NetGPotential <-  GrossGPotential - CurrentSWCost
                        CurrentGPotential <- CurrentCW + NetGPotential
                        
                        if(CurrentGPotential >= LifetimeMaxCW){
                              LifetimeMaxCW <- CurrentGPotential
                              CurrentSRisk <- 0
                        }else{
                              CWLTreshold <- LifetimeMaxCW / 2
                              
                              if(CurrentGPotential <= CWLTreshold){
                                    CurrentSRisk <- 1
                                    CurrentGPotential <- CWLTreshold
                              }else{
                                    CWLoss <- (50 / (LifetimeMaxCW - CWLTreshold)) * (LifetimeMaxCW - CurrentGPotential)
                                    CurrentSRisk <- 0.02 * CWLoss
                              }
                        }
                        
                        CurrentCW <- CurrentGPotential
                        CurrentDepth <- FutureDepth
                        CurrentPProb <- .subset2(PProb, CurrentPing)[CurrentDepth]
                        CurrentPRisk <- CurrentPProb * ((CurrentCW * RAMSlo) + RAMInt)
                        Suvivorship <- Suvivorship * (1 - (CurrentPRisk + CurrentSRisk))
                        
                  }else{
                        #1.Decision making
                        VIRadRange <- .subset2(VIrad, CurrentPing)[SearchRange]
                        PRiskMinRange <-  which(VIRadRange <= DVAdj)
                        
                        if(length(PRiskMinRange) == 0){
                              PRiskMinPos <- which.max(SearchRange)
                              FutureDepth <- SearchRange[PRiskMinPos]
                        }else if(length(PRiskMinRange) == 1){
                              PRiskMinPos <- PRiskMinRange
                              FutureDepth <- SearchRange[PRiskMinPos]
                        }else{
                              PRiskMinPos <- PRiskMinRange[1]
                              FutureDepth <- SearchRange[PRiskMinPos]
                        }
                        
                        #2.Estimations
                        CurrentTemp <- .subset2(Temperature, CurrentPing)[FutureDepth]
                        CurrentGmax <- GEF * exp(0.110 * CurrentTemp)
                        CurrentFA <- .subset2(FA, CurrentPing)[FutureDepth]
                        CurrentFAC <- CurrentFA * (FCF / 1000)
                        a <- 0.7
                        b <- 1.777 * exp(0.234 * CurrentTemp)
                        n <- 0.681 * exp(0.0199 * CurrentTemp)
                        k <- 0.375 * exp(0.0546 * CurrentTemp)
                        m <- 0.858 * exp(-0.008 * CurrentTemp)
                        CurrentFA_C <- (k * CurrentTWMG^(m - n)) / (a * b)
                        CurrentFA_S <- CurrentFA_C + (CurrentGmax / (a * b * CurrentTWMG^(n - 1)))
                        
                        if(CurrentFAC >= CurrentFA_S){
                              GrossGPotential <- CurrentGmax * CurrentTWMG
                        }else{
                              GrossGPotential <- (a * b * CurrentTWMG^n * CurrentFAC) - (k * CurrentTWMG^m)
                        }
                        
                        MovingTime <- abs(FutureDepth - CurrentDepth) / MaxDistance
                        k <- 0.375 * exp(0.0546 * CurrentTemp)
                        m <- 0.858 * exp(-0.008 * CurrentTemp)
                        CurrentSWCost <- k * CurrentTWMG^m * MovingTime * SCF
                        
                        NetGPotential <- GrossGPotential - CurrentSWCost
                        CurrentGPotential <- CurrentCW + NetGPotential
                        
                        if(CurrentGPotential >= LifetimeMaxCW){
                              LifetimeMaxCW <- CurrentGPotential
                              CurrentSRisk <- 0
                        }else{
                              CWLTreshold <- LifetimeMaxCW / 2
                              
                              if(CurrentGPotential <= CWLTreshold){
                                    CurrentSRisk <- 1
                                    CurrentGPotential <- CWLTreshold
                              }else{
                                    CWLoss <- (50 / (LifetimeMaxCW - CWLTreshold)) * (LifetimeMaxCW - CurrentGPotential)
                                    CurrentSRisk <- 0.02 * CWLoss
                              }
                        }
                        
                        CurrentCW <- CurrentGPotential
                        CurrentDepth <- FutureDepth
                        CurrentPProb <- .subset2(PProb, CurrentPing)[CurrentDepth]
                        CurrentPRisk <- CurrentPProb * ((CurrentCW * RAMSlo) + RAMInt)
                        Suvivorship <- Suvivorship * (1 - (CurrentPRisk + CurrentSRisk))
                  }
                  
                  if(Suvivorship < 0.0001){
                        LifeState <- as.factor("DEAD")
                        CessasionStage <- S
                        break
                  }else if(CurrentCW >= WCrc[5, 3] & CurrentCW < WCrc[6, 3]){
                        S <- 4
                  }else if(CurrentCW >= WCrc[6, 3] & CurrentCW < WCrc[7, 3]){
                        S <- 5
                  }else if(CurrentCW >= WCrc[7, 3] & CurrentCW < WCrc[8, 3]){
                        S <- 6
                  }else if(CurrentCW >= WCrc[8, 3] & CurrentCW < WCrc[9, 3]){
                        S <- 7
                  }else if(CurrentCW >= WCrc[9, 3] & CurrentCW < WCrc[10, 3]){
                        S <- 8
                  }else if(CurrentCW >= WCrc[10, 3] & CurrentCW < WCrc[11, 3]){
                        S <- 9
                  }else if(CurrentCW >= WCrc[11, 3]){
                        S <- 10
                        break
                  }else{
                        #Proceed iteration
                  }
            }
      }
      
      #####################################
      #Feeding stages with a store (CIV-CV)
      #####################################
      
      if(LifeState == "DEAD"){
            #Do nothing and proceed
      }else{
            repeat{
                  CurrentPing <- CurrentPing + 1
                  CurrentTWMG <- ((CurrentCW + CurrentSW) * WCF) / 1000
                  MaxDistance <- round((3.7036 * CurrentCW^0.5853), 0)
                  MigrationCeiling <- CurrentDepth - MaxDistance
                  MigrationFloor <- CurrentDepth + MaxDistance
                  
                  MigrationCeiling <- ifelse(MigrationCeiling < MinDepth, yes = MinDepth, no = MigrationCeiling)
                  MigrationFloor <- ifelse(MigrationFloor > MaxDepth, yes = MaxDepth, no = MigrationFloor)
                  SearchRange <- MigrationCeiling:MigrationFloor
                  
                  DVAdj <- (CurrentCW * DVModel_PostC4$coefficients[2]) + DVModel_PostC4$coefficients[1]
                  CurrentVI <- .subset2(VIrad, CurrentPing)[CurrentDepth]
                  
                  if(CurrentVI <= DVAdj){
                        #1.Growth potential estimation
                        TempRange <- .subset2(Temperature, CurrentPing)[SearchRange]
                        FARange <- .subset2(FA, CurrentPing)[SearchRange]
                        FACRange <- FARange * (FCF / 1000)
                        GmaxRange <- GEF * exp(0.110 * TempRange)
                        a <- 0.7
                        b <- 1.777 * exp(0.234 * TempRange)
                        n <- 0.681 * exp(0.0199 * TempRange)
                        k <- 0.375 * exp(0.0546 * TempRange)
                        m <- 0.858 * exp(-0.008 * TempRange)
                        FA_CRange <- (k * CurrentTWMG^(m - n)) / (a * b)
                        FA_SRange <- FA_CRange + (GmaxRange / (a * b * CurrentTWMG^(n - 1)))
                        GPotentialRange <- rep(NA, length(FACRange))
                        
                        for(w in 1:length(FACRange)){
                              if(FACRange[w] >= FA_SRange[w]){
                                    GPotentialRange[w] <- GmaxRange[w] * CurrentTWMG
                              }else{
                                    GPotentialRange[w] <- (a * b[w] * CurrentTWMG^n[w] * FACRange[w]) - (k[w] * CurrentTWMG^m[w])
                              }
                        }
                        
                        #2.Decision making
                        GloTempRange <- .subset2(Temperature, CurrentPing)[1:100]
                        GloGmaxRange <- GEF * exp(0.110 * GloTempRange)
                        GloFARange <- .subset2(FA, CurrentPing)[1:100]
                        GloFACRange <- GloFARange * (FCF / 1000)
                        a <- 0.7
                        Glo_b <- 1.777 * exp(0.234 * GloTempRange)
                        Glo_n <- 0.681 * exp(0.0199 * GloTempRange)
                        Glo_k <- 0.375 * exp(0.0546 * GloTempRange)
                        Glo_m <- 0.858 * exp(-0.008 * GloTempRange)
                        GloFA_CRange <- (Glo_k * CurrentTWMG^(Glo_m - Glo_n)) / (a * Glo_b)
                        GloFA_SRange <- GloFA_CRange + (GloGmaxRange / (a * Glo_b * CurrentTWMG^(Glo_n - 1)))
                        
                        for(x in 1:length(GloFACRange)){
                              if(GloFACRange[x] >= GloFA_SRange[x]){
                                    GloGPotentialRange[x] <- GloGmaxRange[x] * CurrentTWMG
                              }else{
                                    GloGPotentialRange[x] <- (a * Glo_b[x] * CurrentTWMG^Glo_n[x] * GloFACRange[x]) - (Glo_k[x] * CurrentTWMG^Glo_m[x])
                              }
                        }
                        
                        GPotentialMaxRange <- which(GPotentialRange == max(GPotentialRange))
                        
                        if(max(GPotentialRange) < max(GloGPotentialRange)){
                              GPotentialMaxPos <- min(GPotentialMaxRange)                       #Single index location
                              FutureDepth <- SearchRange[GPotentialMaxPos]                      #Depth corresponding to the single index location above
                        }else{
                              GPotentialMaxPos <- max(GPotentialMaxRange)                       #Single index location
                              FutureDepth <- SearchRange[GPotentialMaxPos]                      #Depth corresponding to the single index location above
                        }
                        
                        #2.Estimation
                        MovingTime <- abs(FutureDepth - CurrentDepth) / MaxDistance
                        k <- 0.375 * exp(0.0546 * TempRange[GPotentialMaxPos])
                        m <- 0.858 * exp(-0.008 * TempRange[GPotentialMaxPos])
                        CurrentSWCost <- k * CurrentTWMG^m * MovingTime * SCF
                        
                        GrossGPotential <- max(GPotentialRange)
                        NetGPotential <-  GrossGPotential - CurrentSWCost
                        
                        if(NetGPotential >= 0){
                              CurrentSPotential <- CurrentSW + (NetGPotential * GA)             #Allocation to store
                              NetGPotential <- NetGPotential * (1 - GA)
                        }else{
                              if(CurrentSW >= abs(NetGPotential)){
                                    CurrentSPotential <- CurrentSW - abs(NetGPotential)         #Loss compensated 100% from store
                                    NetGPotential <- 0
                              }else{
                                    CurrentSPotential <- 0
                                    NetGPotential <- NetGPotential + CurrentSW                  #Loss compensated to max.extent if possible
                              }
                        }
                        
                        CurrentGPotential <- CurrentCW + NetGPotential
                        
                        if(CurrentGPotential >= LifetimeMaxCW){
                              LifetimeMaxCW <- CurrentGPotential
                              CurrentSRisk <- 0
                        }else{
                              CWLTreshold <- LifetimeMaxCW / 2
                              
                              if(CurrentGPotential <= CWLTreshold){
                                    CurrentSRisk <- 1
                                    CurrentGPotential <- CWLTreshold
                              }else{
                                    CWLoss <- (50 / (LifetimeMaxCW - CWLTreshold)) * (LifetimeMaxCW - CurrentGPotential)
                                    CurrentSRisk <- 0.02 * CWLoss
                              }
                        }
                        
                        CurrentCW <- CurrentGPotential
                        CurrentSW <- CurrentSPotential
                        CurrentDepth <- FutureDepth
                        CurrentPProb <- .subset2(PProb, CurrentPing)[CurrentDepth]
                        CurrentPRisk <- CurrentPProb * ((CurrentCW * RAMSlo) + RAMInt)
                        Suvivorship <- Suvivorship * (1 - (CurrentPRisk + CurrentSRisk))
                        
                  }else{
                        #1.Decision making
                        VIRadRange <- .subset2(VIrad, CurrentPing)[SearchRange]
                        PRiskMinRange <-  which(VIRadRange <= DVAdj)
                        
                        if(length(PRiskMinRange) == 0){
                              PRiskMinPos <- which.max(SearchRange)
                              FutureDepth <- SearchRange[PRiskMinPos]
                        }else if(length(PRiskMinRange) == 1){
                              PRiskMinPos <- PRiskMinRange
                              FutureDepth <- SearchRange[PRiskMinPos]
                        }else{
                              PRiskMinPos <- PRiskMinRange[1]
                              FutureDepth <- SearchRange[PRiskMinPos]
                        }
                        
                        #2.Estimations
                        CurrentTemp <- .subset2(Temperature, CurrentPing)[FutureDepth]
                        CurrentGmax <- GEF * exp(0.110 * CurrentTemp)
                        CurrentFA <- .subset2(FA, CurrentPing)[FutureDepth]
                        CurrentFAC <- CurrentFA * (FCF / 1000)
                        a <- 0.7
                        b <- 1.777 * exp(0.234 * CurrentTemp)
                        n <- 0.681 * exp(0.0199 * CurrentTemp)
                        k <- 0.375 * exp(0.0546 * CurrentTemp)
                        m <- 0.858 * exp(-0.008 * CurrentTemp)
                        CurrentFA_C <- (k * CurrentTWMG^(m - n)) / (a * b)
                        CurrentFA_S <- CurrentFA_C + (CurrentGmax / (a * b * CurrentTWMG^(n - 1)))
                        
                        if(CurrentFAC >= CurrentFA_S){
                              GrossGPotential <- CurrentGmax * CurrentTWMG
                        }else{
                              GrossGPotential <- (a * b * CurrentTWMG^n * CurrentFAC) - (k * CurrentTWMG^m)
                        }
                        
                        MovingTime <- abs(FutureDepth - CurrentDepth) / MaxDistance
                        k <- 0.375 * exp(0.0546 * CurrentTemp)
                        m <- 0.858 * exp(-0.008 * CurrentTemp)
                        CurrentSWCost <- k * CurrentTWMG^m * MovingTime * SCF
                        
                        NetGPotential <- GrossGPotential - CurrentSWCost
                        
                        if(NetGPotential >= 0){
                              CurrentSPotential <- CurrentSW + (NetGPotential * GA)             #Allocation to store
                              NetGPotential <- NetGPotential * (1 - GA)
                        }else{
                              if(CurrentSW >= abs(NetGPotential)){
                                    CurrentSPotential <- CurrentSW - abs(NetGPotential)         #Loss compensated 100% from store
                                    NetGPotential <- 0
                              }else{
                                    CurrentSPotential <- 0
                                    NetGPotential <- NetGPotential + CurrentSW                  #Loss compensated to max.extent if possible
                              }
                        }
                        
                        CurrentGPotential <- CurrentCW + NetGPotential
                        
                        if(CurrentGPotential >= LifetimeMaxCW){
                              LifetimeMaxCW <- CurrentGPotential
                              CurrentSRisk <- 0
                        }else{
                              CWLTreshold <- LifetimeMaxCW / 2
                              
                              if(CurrentGPotential <= CWLTreshold){
                                    CurrentSRisk <- 1
                                    CurrentGPotential <- CWLTreshold
                              }else{
                                    CWLoss <- (50 / (LifetimeMaxCW - CWLTreshold)) * (LifetimeMaxCW - CurrentGPotential)
                                    CurrentSRisk <- 0.02 * CWLoss
                              }
                        }
                        
                        CurrentCW <- CurrentGPotential
                        CurrentSW <- CurrentSPotential
                        CurrentDepth <- FutureDepth
                        CurrentPProb <- .subset2(PProb, CurrentPing)[CurrentDepth]
                        CurrentPRisk <- CurrentPProb * ((CurrentCW * RAMSlo) + RAMInt)
                        Suvivorship <- Suvivorship * (1 - (CurrentPRisk + CurrentSRisk))
                  }
                  
                  if(Suvivorship < 0.0001){
                        LifeState <- as.factor("DEAD")
                        CessasionStage <- S
                        break
                  }else if((CurrentSW / (CurrentCW * WCF)) >= SD){
                        DiapauseStage <- S
                        break
                  }else if(CurrentCW >= WCrc[12, 3] & CurrentCW < WCrc[13, 3]){
                        S <- 11
                  }else if(CurrentCW >= WCrc[13, 3]){
                        S <- 12
                        break
                  }else{
                        #Proceed iteration
                  }
            }
      }
      
      ####################
      #Diapause simulation
      ####################
      
      if(LifeState == "DEAD"){
            #Do nothing and proceed
      }else{
            if(is.na(DiapauseStage) == TRUE){
                  #No diapause - direct development to CVI-F
            }else{
                  ###############
                  #Diapause entry
                  ###############
                  
                  repeat{
                        CurrentPing <- CurrentPing + 1
                        CurrentTWMG <- ((CurrentCW + CurrentSW) * WCF) / 1000
                        MaxDistance <- round((3.7036 * CurrentCW^0.5853), 0)
                        MigrationCeiling <- CurrentDepth - MaxDistance
                        MigrationFloor <- CurrentDepth + MaxDistance
                        
                        MigrationCeiling <- ifelse(MigrationCeiling < MinDepth, yes = MinDepth, no = MigrationCeiling)
                        MigrationFloor <- ifelse(MigrationFloor > MaxDepth, yes = MaxDepth, no = MigrationFloor)
                        SearchRange <- MigrationCeiling:MigrationFloor
                        
                        #1.Decision making
                        CurrentOD <- which(SearchRange == OD)
                        
                        if(length(CurrentOD) == 1){
                              FutureDepth <- OD
                        }else{
                              if(max(SearchRange) > OD){
                                    FutureDepth <- SearchRange[1]
                              }else{
                                    FutureDepth <- SearchRange[length(SearchRange)]
                              }
                        }
                        
                        #2.Estimation
                        CurrentTemp <- .subset2(Temperature, CurrentPing)[FutureDepth]
                        CurrentGmax <- GEF * exp(0.110 * CurrentTemp)
                        CurrentFA <- .subset2(FA, CurrentPing)[FutureDepth]
                        CurrentFAC <- CurrentFA * (FCF / 1000)
                        a <- 0.7
                        b <- 1.777 * exp(0.234 * CurrentTemp)
                        n <- 0.681 * exp(0.0199 * CurrentTemp)
                        k <- 0.375 * exp(0.0546 * CurrentTemp)
                        m <- 0.858 * exp(-0.008 * CurrentTemp)
                        CurrentFA_C <- (k * CurrentTWMG^(m - n)) / (a * b)
                        CurrentFA_S <- CurrentFA_C + (CurrentGmax / (a * b * CurrentTWMG^(n - 1)))
                        
                        if(CurrentFAC >= CurrentFA_S){
                              GrossGPotential <- CurrentGmax * CurrentTWMG
                        }else{
                              GrossGPotential <- (a * b * CurrentTWMG^n * CurrentFAC) - (k * CurrentTWMG^m)
                        }
                        
                        MovingTime <- abs(FutureDepth - CurrentDepth) / MaxDistance
                        k <- 0.375 * exp(0.0546 * CurrentTemp)
                        m <- 0.858 * exp(-0.008 * CurrentTemp)
                        CurrentSWCost <- k * CurrentTWMG^m * MovingTime * SCF
                        
                        NetGPotential <- GrossGPotential - CurrentSWCost
                        
                        if(NetGPotential >= 0){
                              CurrentSPotential <- CurrentSW + (NetGPotential * GA)             #Allocation to store
                              NetGPotential <- NetGPotential * (1 - GA)
                        }else{
                              if(CurrentSW >= abs(NetGPotential)){
                                    CurrentSPotential <- CurrentSW - abs(NetGPotential)         #Loss compensated 100% from store
                                    NetGPotential <- 0
                              }else{
                                    CurrentSPotential <- 0
                                    NetGPotential <- NetGPotential + CurrentSW                  #Loss compensated to max.extent if possible
                              }
                        }
                        
                        CurrentGPotential <- CurrentCW + NetGPotential
                        
                        if(CurrentGPotential >= LifetimeMaxCW){
                              LifetimeMaxCW <- CurrentGPotential
                              CurrentSRisk <- 0
                        }else{
                              CWLTreshold <- LifetimeMaxCW / 2
                              
                              if(CurrentGPotential <= CWLTreshold){
                                    CurrentSRisk <- 1
                                    CurrentGPotential <- CWLTreshold
                              }else{
                                    CWLoss <- (50 / (LifetimeMaxCW - CWLTreshold)) * (LifetimeMaxCW - CurrentGPotential)
                                    CurrentSRisk <- 0.02 * CWLoss
                              }
                        }
                        
                        CurrentCW <- CurrentGPotential
                        CurrentSW <- CurrentSPotential
                        CurrentDepth <- FutureDepth
                        CurrentPProb <- .subset2(PProb, CurrentPing)[CurrentDepth]
                        CurrentPRisk <- CurrentPProb * ((CurrentCW * RAMSlo) + RAMInt)
                        Suvivorship <- Suvivorship * (1 - (CurrentPRisk + CurrentSRisk))
                        
                        if(Suvivorship < 0.0001){
                              LifeState <- as.factor("DEAD")
                              CessasionStage <- S
                              break
                        }else if(CurrentDepth == OD){
                              DiapauseStartTime <- CurrentPing
                              DiapauseStartSW <- CurrentSW
                              break
                        }else{
                              #Proceed iteration
                        }
                  }
                  
                  #####################
                  #Diapause progression
                  #####################
                  
                  if(LifeState == "DEAD"){
                        #Do nothing and proceed
                  }else{
                        CurrentTemp <- .subset2(ODParameters, 2)[ODParameters$Depths == OD]
                        CurrentPProb <- .subset2(ODParameters, 3)[ODParameters$Depths == OD]
                        
                        k <- 0.375 * exp(0.0546 * CurrentTemp)
                        m <- 0.858 * exp(-0.008 * CurrentTemp)
                        CurrentTWMG <- ((CurrentCW + CurrentSW) * WCF) / 1000
                        MCost_ceiling <- k * CurrentTWMG^m * MRF
                        ExpendedSW <- CurrentSW - (CurrentSW * (1 - SA))
                        EndSW <- CurrentSW - ExpendedSW
                        EndTWMG <- ((CurrentCW + EndSW) * WCF) / 1000
                        MCost_floor <- k * EndTWMG^m * MRF
                        MeanMCost <- (MCost_ceiling + MCost_floor) / 2
                        DiapauseHrs <- round(ExpendedSW / MeanMCost, 0)
                        round(ExpendedSW / MCost_floor, 0)
                        CurrentSW <- EndSW
                        
                        CurrentPRisk <- CurrentPProb * ((CurrentCW * RAMSlo) + RAMInt)
                        DiapauseSuvivorship <- (1 - CurrentPRisk)^DiapauseHrs
                        Suvivorship <- Suvivorship * DiapauseSuvivorship
                        
                        if(Suvivorship < 0.0001){
                              LifeState <- as.factor("DEAD")
                              CessasionStage <- S
                        }else{
                              #Do nothing
                        }
                        
                        DiapauseEndTime <- CurrentPing + DiapauseHrs
                        CurrentPing <- DiapauseEndTime
                  }
                  
                  ##############
                  #Diapause exit
                  ##############
                  
                  if(LifeState == "DEAD"){
                        #Do nothing and proceed
                  }else{
                        repeat{
                              CurrentPing <- CurrentPing + 1
                              CurrentTWMG <- ((CurrentCW + CurrentSW) * WCF) / 1000
                              MaxDistance <- round((3.7036 * CurrentCW^0.5853), 0)
                              MigrationCeiling <- CurrentDepth - MaxDistance
                              MigrationFloor <- CurrentDepth + MaxDistance
                              
                              MigrationCeiling <- ifelse(MigrationCeiling < MinDepth, yes = MinDepth, no = MigrationCeiling)
                              MigrationFloor <- ifelse(MigrationFloor > MaxDepth, yes = MaxDepth, no = MigrationFloor)
                              SearchRange <- MigrationCeiling:MigrationFloor
                              
                              DVAdj <- (CurrentCW * DVModel_PostC4$coefficients[2]) + DVModel_PostC4$coefficients[1]
                              CurrentVI <- .subset2(VIrad, CurrentPing)[CurrentDepth]
                              
                              if(CurrentVI <= DVAdj){
                                    #1.Growth potential estimation
                                    TempRange <- .subset2(Temperature, CurrentPing)[SearchRange]
                                    FARange <- .subset2(FA, CurrentPing)[SearchRange]
                                    FACRange <- FARange * (FCF / 1000)
                                    GmaxRange <- GEF * exp(0.110 * TempRange)
                                    a <- 0.7
                                    b <- 1.777 * exp(0.234 * TempRange)
                                    n <- 0.681 * exp(0.0199 * TempRange)
                                    k <- 0.375 * exp(0.0546 * TempRange)
                                    m <- 0.858 * exp(-0.008 * TempRange)
                                    FA_CRange <- (k * CurrentTWMG^(m - n)) / (a * b)
                                    FA_SRange <- FA_CRange + (GmaxRange / (a * b * CurrentTWMG^(n - 1)))
                                    GPotentialRange <- rep(NA, length(FACRange))
                                    
                                    for(w in 1:length(FACRange)){
                                          if(FACRange[w] >= FA_SRange[w]){
                                                GPotentialRange[w] <- GmaxRange[w] * CurrentTWMG
                                          }else{
                                                GPotentialRange[w] <- (a * b[w] * CurrentTWMG^n[w] * FACRange[w]) - (k[w] * CurrentTWMG^m[w])
                                          }
                                    }
                                    
                                    #2.Decision making
                                    GloTempRange <- .subset2(Temperature, CurrentPing)[1:100]
                                    GloGmaxRange <- GEF * exp(0.110 * GloTempRange)
                                    GloFARange <- .subset2(FA, CurrentPing)[1:100]
                                    GloFACRange <- GloFARange * (FCF / 1000)
                                    a <- 0.7
                                    Glo_b <- 1.777 * exp(0.234 * GloTempRange)
                                    Glo_n <- 0.681 * exp(0.0199 * GloTempRange)
                                    Glo_k <- 0.375 * exp(0.0546 * GloTempRange)
                                    Glo_m <- 0.858 * exp(-0.008 * GloTempRange)
                                    GloFA_CRange <- (Glo_k * CurrentTWMG^(Glo_m - Glo_n)) / (a * Glo_b)
                                    GloFA_SRange <- GloFA_CRange + (GloGmaxRange / (a * Glo_b * CurrentTWMG^(Glo_n - 1)))
                                    
                                    for(x in 1:length(GloFACRange)){
                                          if(GloFACRange[x] >= GloFA_SRange[x]){
                                                GloGPotentialRange[x] <- GloGmaxRange[x] * CurrentTWMG
                                          }else{
                                                GloGPotentialRange[x] <- (a * Glo_b[x] * CurrentTWMG^Glo_n[x] * GloFACRange[x]) - (Glo_k[x] * CurrentTWMG^Glo_m[x])
                                          }
                                    }
                                    
                                    GPotentialMaxRange <- which(GPotentialRange == max(GPotentialRange))
                                    
                                    if(max(GPotentialRange) < max(GloGPotentialRange)){
                                          GPotentialMaxPos <- min(GPotentialMaxRange)                       #Single index location
                                          FutureDepth <- SearchRange[GPotentialMaxPos]                      #Depth corresponding to the single index location above
                                    }else{
                                          GPotentialMaxPos <- max(GPotentialMaxRange)                       #Single index location
                                          FutureDepth <- SearchRange[GPotentialMaxPos]                      #Depth corresponding to the single index location above
                                    }
                                    
                                    #2.Estimation
                                    MovingTime <- abs(FutureDepth - CurrentDepth) / MaxDistance
                                    k <- 0.375 * exp(0.0546 * TempRange[GPotentialMaxPos])
                                    m <- 0.858 * exp(-0.008 * TempRange[GPotentialMaxPos])
                                    CurrentSWCost <- k * CurrentTWMG^m * MovingTime * SCF
                                    
                                    GrossGPotential <- max(GPotentialRange)
                                    NetGPotential <-  GrossGPotential - CurrentSWCost
                                    
                                    if(NetGPotential >= 0){
                                          CurrentSPotential <- CurrentSW + 0                                #No allocation to store after diapause (1-yr life cycle)
                                          NetGPotential <- NetGPotential                                    #All gains to somatic growth
                                    }else{
                                          MaxStoGAllocation <- GEF * exp(0.110 * .subset2(Temperature, CurrentPing)[FutureDepth]) * CurrentTWMG
                                          
                                          if(CurrentSW >= MaxStoGAllocation){
                                                CurrentSPotential <- CurrentSW - MaxStoGAllocation
                                                NetGPotential <- MaxStoGAllocation + NetGPotential
                                          }else if(CurrentSW < MaxStoGAllocation & CurrentSW >= abs(NetGPotential)){
                                                CurrentSPotential <- CurrentSW - abs(NetGPotential)
                                                NetGPotential <- 0
                                          }else{
                                                CurrentSPotential <- 0
                                                NetGPotential <- NetGPotential + CurrentSW
                                          }
                                    }
                                    
                                    CurrentGPotential <- CurrentCW + NetGPotential
                                    
                                    if(CurrentGPotential >= LifetimeMaxCW){
                                          LifetimeMaxCW <- CurrentGPotential
                                          CurrentSRisk <- 0
                                    }else{
                                          CWLTreshold <- LifetimeMaxCW / 2
                                          
                                          if(CurrentGPotential <= CWLTreshold){
                                                CurrentSRisk <- 1
                                                CurrentGPotential <- CWLTreshold
                                          }else{
                                                CWLoss <- (50 / (LifetimeMaxCW - CWLTreshold)) * (LifetimeMaxCW - CurrentGPotential)
                                                CurrentSRisk <- 0.02 * CWLoss
                                          }
                                    }
                                    
                                    CurrentCW <- CurrentGPotential
                                    CurrentSW <- CurrentSPotential
                                    CurrentDepth <- FutureDepth
                                    CurrentPProb <- .subset2(PProb, CurrentPing)[CurrentDepth]
                                    CurrentPRisk <- CurrentPProb * ((CurrentCW * RAMSlo) + RAMInt)
                                    Suvivorship <- Suvivorship * (1 - (CurrentPRisk + CurrentSRisk))
                                    
                              }else{
                                    #Minimize risk
                                    #1.Decision making
                                    VIRadRange <- .subset2(VIrad, CurrentPing)[SearchRange]
                                    PRiskMinRange <-  which(VIRadRange <= DVAdj)
                                    
                                    if(length(PRiskMinRange) == 0){
                                          PRiskMinPos <- which.max(SearchRange)
                                          FutureDepth <- SearchRange[PRiskMinPos]
                                    }else if(length(PRiskMinRange) == 1){
                                          PRiskMinPos <- PRiskMinRange
                                          FutureDepth <- SearchRange[PRiskMinPos]
                                    }else{
                                          PRiskMinPos <- PRiskMinRange[1]
                                          FutureDepth <- SearchRange[PRiskMinPos]
                                    }
                                    
                                    #2.Estimations
                                    CurrentTemp <- .subset2(Temperature, CurrentPing)[FutureDepth]
                                    CurrentGmax <- GEF * exp(0.110 * CurrentTemp)
                                    CurrentFA <- .subset2(FA, CurrentPing)[FutureDepth]
                                    CurrentFAC <- CurrentFA * (FCF / 1000)
                                    a <- 0.7
                                    b <- 1.777 * exp(0.234 * CurrentTemp)
                                    n <- 0.681 * exp(0.0199 * CurrentTemp)
                                    k <- 0.375 * exp(0.0546 * CurrentTemp)
                                    m <- 0.858 * exp(-0.008 * CurrentTemp)
                                    CurrentFA_C <- (k * CurrentTWMG^(m - n)) / (a * b)
                                    CurrentFA_S <- CurrentFA_C + (CurrentGmax / (a * b * CurrentTWMG^(n - 1)))
                                    
                                    if(CurrentFAC >= CurrentFA_S){
                                          GrossGPotential <- CurrentGmax * CurrentTWMG
                                    }else{
                                          GrossGPotential <- (a * b * CurrentTWMG^n * CurrentFAC) - (k * CurrentTWMG^m)
                                    }
                                    
                                    MovingTime <- abs(FutureDepth - CurrentDepth) / MaxDistance
                                    k <- 0.375 * exp(0.0546 * CurrentTemp)
                                    m <- 0.858 * exp(-0.008 * CurrentTemp)
                                    CurrentSWCost <- k * CurrentTWMG^m * MovingTime * SCF
                                    
                                    NetGPotential <- GrossGPotential - CurrentSWCost
                                    
                                    if(NetGPotential >= 0){
                                          CurrentSPotential <- CurrentSW + 0                                #No allocation to store after diapause (1-yr life cycle)
                                          NetGPotential <- NetGPotential                                    #All gains to somatic growth
                                    }else{
                                          MaxStoGAllocation <- GEF * exp(0.110 * .subset2(Temperature, CurrentPing)[FutureDepth]) * CurrentTWMG
                                          
                                          if(CurrentSW >= MaxStoGAllocation){
                                                CurrentSPotential <- CurrentSW - MaxStoGAllocation
                                                NetGPotential <- MaxStoGAllocation + NetGPotential
                                          }else if(CurrentSW < MaxStoGAllocation & CurrentSW >= abs(NetGPotential)){
                                                CurrentSPotential <- CurrentSW - abs(NetGPotential)
                                                NetGPotential <- 0
                                          }else{
                                                CurrentSPotential <- 0
                                                NetGPotential <- NetGPotential + CurrentSW
                                          }
                                    }
                                    
                                    CurrentGPotential <- CurrentCW + NetGPotential
                                    
                                    if(CurrentGPotential >= LifetimeMaxCW){
                                          LifetimeMaxCW <- CurrentGPotential
                                          CurrentSRisk <- 0
                                    }else{
                                          CWLTreshold <- LifetimeMaxCW / 2
                                          
                                          if(CurrentGPotential <= CWLTreshold){
                                                CurrentSRisk <- 1
                                                CurrentGPotential <- CWLTreshold
                                          }else{
                                                CWLoss <- (50 / (LifetimeMaxCW - CWLTreshold)) * (LifetimeMaxCW - CurrentGPotential)
                                                CurrentSRisk <- 0.02 * CWLoss
                                          }
                                    }
                                    
                                    CurrentCW <- CurrentGPotential
                                    CurrentSW <- CurrentSPotential
                                    CurrentDepth <- FutureDepth
                                    CurrentPProb <- .subset2(PProb, CurrentPing)[CurrentDepth]
                                    CurrentPRisk <- CurrentPProb * ((CurrentCW * RAMSlo) + RAMInt)
                                    Suvivorship <- Suvivorship * (1 - (CurrentPRisk + CurrentSRisk))
                              }
                              
                              if(Suvivorship < 0.0001){
                                    LifeState <- as.factor("DEAD")
                                    CessasionStage <- S
                                    break
                              }else if(CurrentCW >= WCrc[12, 3] & CurrentCW < WCrc[13, 3]){
                                    S <- 11
                              }else if(CurrentCW >= WCrc[13, 3]){
                                    S <- 12
                                    break
                              }else{
                                    #Proceed iteration
                              }
                        }
                  }
                  
                  ########################
                  #Female and reproduction
                  ########################
                  
                  if(LifeState == "DEAD"){
                        #Do nothing and proceed
                  }else{
                        if(DiapauseEndTime <= 8760){
                              #Zero fitness Case-II: Don't process
                        }else{
                              repeat{
                                    CurrentPing <- CurrentPing + 1
                                    CurrentTWMG <- ((CurrentCW + CurrentSW) * WCF) / 1000
                                    MaxDistance <- round((3.7036 * CurrentCW^0.5853), 0)
                                    MigrationCeiling <- CurrentDepth - MaxDistance
                                    MigrationFloor <- CurrentDepth + MaxDistance
                                    
                                    MigrationCeiling <- ifelse(MigrationCeiling < MinDepth, yes = MinDepth, no = MigrationCeiling)
                                    MigrationFloor <- ifelse(MigrationFloor > MaxDepth, yes = MaxDepth, no = MigrationFloor)
                                    SearchRange <- MigrationCeiling:MigrationFloor
                                    
                                    DVAdj <- (CurrentCW * DVModel_PostC4$coefficients[2]) + DVModel_PostC4$coefficients[1]
                                    CurrentVI <- .subset2(VIrad, CurrentPing)[CurrentDepth]
                                    
                                    if(CurrentVI <= DVAdj){
                                          #1.Growth potential estimation
                                          TempRange <- .subset2(Temperature, CurrentPing)[SearchRange]
                                          FARange <- .subset2(FA, CurrentPing)[SearchRange]
                                          FACRange <- FARange * (FCF / 1000)
                                          GmaxRange <- GEF * exp(0.110 * TempRange)
                                          a <- 0.7
                                          b <- 1.777 * exp(0.234 * TempRange)
                                          n <- 0.681 * exp(0.0199 * TempRange)
                                          k <- 0.375 * exp(0.0546 * TempRange)
                                          m <- 0.858 * exp(-0.008 * TempRange)
                                          FA_CRange <- (k * CurrentTWMG^(m - n)) / (a * b)
                                          FA_SRange <- FA_CRange + (GmaxRange / (a * b * CurrentTWMG^(n - 1)))
                                          GPotentialRange <- rep(NA, length(FACRange))
                                          
                                          for(w in 1:length(FACRange)){
                                                if(FACRange[w] >= FA_SRange[w]){
                                                      GPotentialRange[w] <- GmaxRange[w] * CurrentTWMG
                                                }else{
                                                      GPotentialRange[w] <- (a * b[w] * CurrentTWMG^n[w] * FACRange[w]) - (k[w] * CurrentTWMG^m[w])
                                                }
                                          }
                                          
                                          #2.Decision making
                                          GloTempRange <- .subset2(Temperature, CurrentPing)[1:100]
                                          GloGmaxRange <- GEF * exp(0.110 * GloTempRange)
                                          GloFARange <- .subset2(FA, CurrentPing)[1:100]
                                          GloFACRange <- GloFARange * (FCF / 1000)
                                          a <- 0.7
                                          Glo_b <- 1.777 * exp(0.234 * GloTempRange)
                                          Glo_n <- 0.681 * exp(0.0199 * GloTempRange)
                                          Glo_k <- 0.375 * exp(0.0546 * GloTempRange)
                                          Glo_m <- 0.858 * exp(-0.008 * GloTempRange)
                                          GloFA_CRange <- (Glo_k * CurrentTWMG^(Glo_m - Glo_n)) / (a * Glo_b)
                                          GloFA_SRange <- GloFA_CRange + (GloGmaxRange / (a * Glo_b * CurrentTWMG^(Glo_n - 1)))
                                          
                                          for(x in 1:length(GloFACRange)){
                                                if(GloFACRange[x] >= GloFA_SRange[x]){
                                                      GloGPotentialRange[x] <- GloGmaxRange[x] * CurrentTWMG
                                                }else{
                                                      GloGPotentialRange[x] <- (a * Glo_b[x] * CurrentTWMG^Glo_n[x] * GloFACRange[x]) - (Glo_k[x] * CurrentTWMG^Glo_m[x])
                                                }
                                          }
                                          
                                          GPotentialMaxRange <- which(GPotentialRange == max(GPotentialRange))
                                          
                                          if(max(GPotentialRange) < max(GloGPotentialRange)){
                                                GPotentialMaxPos <- min(GPotentialMaxRange)                       #Single index location
                                                FutureDepth <- SearchRange[GPotentialMaxPos]                      #Depth corresponding to the single index location above
                                          }else{
                                                GPotentialMaxPos <- max(GPotentialMaxRange)                       #Single index location
                                                FutureDepth <- SearchRange[GPotentialMaxPos]                      #Depth corresponding to the single index location above
                                          }
                                          
                                          #2.Estimation
                                          MovingTime <- abs(FutureDepth - CurrentDepth) / MaxDistance
                                          k <- 0.375 * exp(0.0546 * TempRange[GPotentialMaxPos])
                                          m <- 0.858 * exp(-0.008 * TempRange[GPotentialMaxPos])
                                          CurrentSWCost <- k * CurrentTWMG^m * MovingTime * SCF
                                          
                                          GrossGPotential <- max(GPotentialRange)
                                          NetGPotential <-  GrossGPotential - CurrentSWCost
                                          
                                          if(NetGPotential >= 0){
                                                MaxStoRAllocation <- GEF * exp(0.110 * .subset2(Temperature, CurrentPing)[FutureDepth]) * CurrentTWMG
                                                
                                                if(CurrentSW >= MaxStoRAllocation){
                                                      CurrentSPotential <- CurrentSW - MaxStoRAllocation
                                                      CurrentRPotential <- NetGPotential + MaxStoRAllocation
                                                      NetGPotential <- 0
                                                }else{
                                                      CurrentSPotential <- 0
                                                      CurrentRPotential <- CurrentSW + NetGPotential
                                                      NetGPotential <- 0
                                                }
                                          }else{
                                                MaxStoRAllocation <- GEF * exp(0.110 * .subset2(Temperature, CurrentPing)[FutureDepth]) * CurrentTWMG
                                                
                                                if(CurrentSW >= MaxStoRAllocation){
                                                      CurrentSPotential <- CurrentSW - MaxStoRAllocation
                                                      CurrentRPotential <- MaxStoRAllocation - abs(NetGPotential)
                                                      NetGPotential <- 0
                                                }else if(CurrentSW < MaxStoRAllocation & CurrentSW >= abs(NetGPotential)){
                                                      CurrentSPotential <- CurrentSW - abs(NetGPotential)
                                                      CurrentRPotential <- 0
                                                      NetGPotential <- 0
                                                }else{
                                                      CurrentSPotential <- 0
                                                      CurrentRPotential <- 0
                                                      NetGPotential <- NetGPotential + CurrentSW
                                                }
                                          }
                                          
                                          CurrentGPotential <- CurrentCW + NetGPotential
                                          NEggs <- round(CurrentRPotential / UEM, 0)
                                          CurrentRPotential <- CurrentRPotential - (NEggs * UEM)
                                          
                                          if(CurrentGPotential >= LifetimeMaxCW){
                                                LifetimeMaxCW <- CurrentGPotential
                                                CurrentSRisk <- 0
                                          }else{
                                                CWLTreshold <- LifetimeMaxCW / 2
                                                
                                                if(CurrentGPotential <= CWLTreshold){
                                                      CurrentSRisk <- 1
                                                      CurrentGPotential <- CWLTreshold
                                                }else{
                                                      CWLoss <- (50 / (LifetimeMaxCW - CWLTreshold)) * (LifetimeMaxCW - CurrentGPotential)
                                                      CurrentSRisk <- 0.02 * CWLoss
                                                }
                                          }
                                          
                                          CurrentCW <- CurrentGPotential
                                          CurrentSW <- CurrentSPotential
                                          CurrentRW <- CurrentRPotential
                                          CurrentDepth <- FutureDepth
                                          CurrentPProb <- .subset2(PProb, CurrentPing)[CurrentDepth]
                                          CurrentPRisk <- CurrentPProb * ((CurrentCW * RAMSlo) + RAMInt)
                                          CurrentEggProduction <- CurrentEggProduction + NEggs
                                          
                                          if(CurrentEggProduction > 0 & FEPState == 0){
                                                FEPState <- 1
                                                FEPTime <- CurrentPing
                                          }else{
                                                #Do nothing
                                          }
                                          
                                          Suvivorship <- Suvivorship * (1 - (CurrentPRisk + CurrentSRisk))
                                          CurrentFitness <- CurrentFitness + (NEggs * Suvivorship)
                                          
                                    }else{
                                          #1.Decision making
                                          VIRadRange <- .subset2(VIrad, CurrentPing)[SearchRange]
                                          PRiskMinRange <-  which(VIRadRange <= DVAdj)
                                          
                                          if(length(PRiskMinRange) == 0){
                                                PRiskMinPos <- which.max(SearchRange)
                                                FutureDepth <- SearchRange[PRiskMinPos]
                                          }else if(length(PRiskMinRange) == 1){
                                                PRiskMinPos <- PRiskMinRange
                                                FutureDepth <- SearchRange[PRiskMinPos]
                                          }else{
                                                PRiskMinPos <- PRiskMinRange[1]
                                                FutureDepth <- SearchRange[PRiskMinPos]
                                          }
                                          
                                          #2.Estimations
                                          CurrentTemp <- .subset2(Temperature, CurrentPing)[FutureDepth]
                                          CurrentGmax <- GEF * exp(0.110 * CurrentTemp)
                                          CurrentFA <- .subset2(FA, CurrentPing)[FutureDepth]
                                          CurrentFAC <- CurrentFA * (FCF / 1000)
                                          a <- 0.7
                                          b <- 1.777 * exp(0.234 * CurrentTemp)
                                          n <- 0.681 * exp(0.0199 * CurrentTemp)
                                          k <- 0.375 * exp(0.0546 * CurrentTemp)
                                          m <- 0.858 * exp(-0.008 * CurrentTemp)
                                          CurrentFA_C <- (k * CurrentTWMG^(m - n)) / (a * b)
                                          CurrentFA_S <- CurrentFA_C + (CurrentGmax / (a * b * CurrentTWMG^(n - 1)))
                                          
                                          if(CurrentFAC >= CurrentFA_S){
                                                GrossGPotential <- CurrentGmax * CurrentTWMG
                                          }else{
                                                GrossGPotential <- (a * b * CurrentTWMG^n * CurrentFAC) - (k * CurrentTWMG^m)
                                          }
                                          
                                          MovingTime <- abs(FutureDepth - CurrentDepth) / MaxDistance
                                          k <- 0.375 * exp(0.0546 * CurrentTemp)
                                          m <- 0.858 * exp(-0.008 * CurrentTemp)
                                          CurrentSWCost <- k * CurrentTWMG^m * MovingTime * SCF
                                          
                                          NetGPotential <- GrossGPotential - CurrentSWCost
                                          
                                          if(NetGPotential >= 0){
                                                MaxStoRAllocation <- GEF * exp(0.110 * .subset2(Temperature, CurrentPing)[FutureDepth]) * CurrentTWMG
                                                
                                                if(CurrentSW >= MaxStoRAllocation){
                                                      CurrentSPotential <- CurrentSW - MaxStoRAllocation
                                                      CurrentRPotential <- NetGPotential + MaxStoRAllocation
                                                      NetGPotential <- 0
                                                }else{
                                                      CurrentSPotential <- 0
                                                      CurrentRPotential <- CurrentSW + NetGPotential
                                                      NetGPotential <- 0
                                                }
                                          }else{
                                                MaxStoRAllocation <- GEF * exp(0.110 * .subset2(Temperature, CurrentPing)[FutureDepth]) * CurrentTWMG
                                                
                                                if(CurrentSW >= MaxStoRAllocation){
                                                      CurrentSPotential <- CurrentSW - MaxStoRAllocation
                                                      CurrentRPotential <- MaxStoRAllocation - abs(NetGPotential)
                                                      NetGPotential <- 0
                                                }else if(CurrentSW < MaxStoRAllocation & CurrentSW >= abs(NetGPotential)){
                                                      CurrentSPotential <- CurrentSW - abs(NetGPotential)
                                                      CurrentRPotential <- 0
                                                      NetGPotential <- 0
                                                }else{
                                                      CurrentSPotential <- 0
                                                      CurrentRPotential <- 0
                                                      NetGPotential <- NetGPotential + CurrentSW
                                                }
                                          }
                                          
                                          CurrentGPotential <- CurrentCW + NetGPotential
                                          NEggs <- round(CurrentRPotential / UEM, 0)
                                          CurrentRPotential <- CurrentRPotential - (NEggs * UEM)
                                          
                                          if(CurrentGPotential >= LifetimeMaxCW){
                                                LifetimeMaxCW <- CurrentGPotential
                                                CurrentSRisk <- 0
                                          }else{
                                                CWLTreshold <- LifetimeMaxCW / 2
                                                
                                                if(CurrentGPotential <= CWLTreshold){
                                                      CurrentSRisk <- 1
                                                      CurrentGPotential <- CWLTreshold
                                                }else{
                                                      CWLoss <- (50 / (LifetimeMaxCW - CWLTreshold)) * (LifetimeMaxCW - CurrentGPotential)
                                                      CurrentSRisk <- 0.02 * CWLoss
                                                }
                                          }
                                          
                                          CurrentCW <- CurrentGPotential
                                          CurrentSW <- CurrentSPotential
                                          CurrentRW <- CurrentRPotential
                                          CurrentDepth <- FutureDepth
                                          CurrentPProb <- .subset2(PProb, CurrentPing)[CurrentDepth]
                                          CurrentPRisk <- CurrentPProb * ((CurrentCW * RAMSlo) + RAMInt)
                                          CurrentEggProduction <- CurrentEggProduction + NEggs
                                          
                                          if(CurrentEggProduction > 0 & FEPState == 0){
                                                FEPState <- 1
                                                FEPTime <- CurrentPing
                                          }else{
                                                #Do nothing
                                          }
                                          
                                          Suvivorship <- Suvivorship * (1 - (CurrentPRisk + CurrentSRisk))
                                          CurrentFitness <- CurrentFitness + (NEggs * Suvivorship)
                                    }
                                    
                                    if(Suvivorship < 0.0001){
                                          LifeState <- as.factor("DEAD")
                                          CessasionStage <- S
                                          break
                                    }else if(CurrentEggProduction >= NBREAK){
                                          CessasionStage <- S
                                          LifeState <- as.factor("ALIVE")
                                          break
                                    }else{
                                          #Proceed iteration
                                    }
                              }
                        }
                  }
            }
      }
      
      ########
      #Outputs
      ########
      DIni[IND] <<- DiapauseStartTime
      DTer[IND] <<- DiapauseEndTime
      DStg[IND] <<- DiapauseStage
      Fecundity[IND] <<- CurrentEggProduction 
      IFitness[IND] <<- CurrentFitness
      LCTer[IND] <<- CurrentPing
      FEp[IND] <<- FEPTime 
      
}

