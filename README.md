# Investigating Item-level and Item-parcelling Models in Structural Equation Modelling
## The R-syntax for the simulation of item-level model
### _Normal_

    library("lavaan") # Required for SEM analysis
    library("matrixcalc") #Required for is.positive.definite function
    library("psych")  #required for reliability function


    data.f=function(n)
  	{
	KSI1=rnorm(n)
	X1Err=rnorm(n)
	X2Err=rnorm(n)
	X3Err=rnorm(n)
	X4Err=rnorm(n)
	X5Err=rnorm(n) 
	X6Err=rnorm(n)
	X7Err=rnorm(n)
	X8Err=rnorm(n)
	X9Err=rnorm(n)
	X10Err=rnorm(n)
	X11Err=rnorm(n)
	X12Err=rnorm(n)
	X1 = round(.70*KSI1 + .714*X1Err + 5.2)
	X2 = round(.80*KSI1 + .060*X2Err + 5.2)
	X3 = round(.90*KSI1 + .436*X3Err + 5.2)
	
	X4 = round(.70*KSI1 + .714*X4Err + 5.2)
	X5 = round(.80*KSI1 + .060*X5Err + 5.2)
	X6 = round(.90*KSI1 + .436*X6Err + 5.2)
	
	X7 = round(.70*KSI1 + .714*X7Err + 5.2)
	X8 = round(.80*KSI1 + .060*X8Err + 5.2)
	X9 = round(.90*KSI1 + .436*X9Err + 5.2)
	
	X10 = round(.70*KSI1 + .714*X10Err + 5.2)
	X11 = round(.80*KSI1 + .060*X11Err + 5.2)
	X12 = round(.90*KSI1 + .436*X12Err + 5.2)

	Eta1Err=rnorm(n)
	Y1Err=rnorm(n)
	Y2Err=rnorm(n)
	Y3Err=rnorm(n)

	ETA1=.585*KSI1 + .811*Eta1Err
	Y1=round(.70*ETA1 + .714*Y1Err + 5.2)
	Y2=round(.80*ETA1 + .060*Y2Err + 5.2)
	Y3=round(.90*ETA1 + .436*Y3Err + 5.2)

	data.norm = data.frame(X1, X2, X3, X4, X5, X6, X7, X8, X9, X10, X11, X12, Y1, Y2, Y3)
	
	return(data.norm)
	}
    -------------------------------------------------------------------------
    n=100 #sample size
    nrep=5000 # set number of replications
    rep<-1
    g<-100000 # set maximum iteration

    beta.norm <- matrix(NA, nrow = g, ncol = 1) #create space to store parameter estimates values
    FIT.norm <- matrix(NA, nrow = g, ncol = 7) #create space to store model fit values
    alphaX <- matrix(NA, nrow = g, ncol = 1) #create space to store reliability values
    alphaY <-  matrix(NA, nrow = g, ncol = 1) #create space to store reliability values

    set.seed(12345)

	for (i in 1:g)
		{

			if(rep <= nrep)
			{
			data.norm = data.f(n)
		
			#Fit the item-level model (Normal)
			model.norm <- '
			KSI1 =~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + X11 + X12
			ETA1 =~ Y1 + Y2 + Y3

			#REGRESSION
			ETA1 ~ KSI1
			'	
			myCov <- cov(data.norm)
			round(myCov, 2)
			# check whether any correlations are perfect (i.e., collinearity)
			myCor <- cov2cor(myCov)
			noDiag <- myCor
			diag(noDiag) <- 0
			corr <- any(noDiag == 1)
			# if not, check for multicollinearity (i.e., is one variable a linear combination of 2+ variables?)
			mult <- det(myCov) < 0
			## or
			eig <- any(eigen(myCov)$values < 0)
			#check for positive-definite matrix
			pd <- !isTRUE(is.positive.definite(myCov, tol=1e-8))

				if(corr=="TRUE" || mult=="TRUE" || eig=="TRUE" || pd=="TRUE" )
				{
				rep <- rep
				alphaX[i] = 0
				alphaY[i] = 0
				alphaXY[i] = 0
				beta.norm[i, 1:1] = 0
				FIT.norm[i,1:7] = cbind(0,0,0,0,0,0,0)
				colnames(FIT.norm) <- c("chisq","gfi", "agfi", "rmsea", "nfi", "nnfi", "cfi")
				} else {
			
				fit.norm <- cfa(model=model.norm, data = data.norm,estimator="ML")
				BN <- parameterEstimates (fit.norm, ci = FALSE, standardized = TRUE)[16:16, 4]
				BN <- as.numeric(BN)

					if(is.na(BN))
					{
					rep <- rep
					alphaX[i] = 0
					alphaY[i] = 0
					alphaXY[i] = 0
					beta.norm[i, 1:1] = 0
					FIT.norm[i,1:7] = cbind(0,0,0,0,0,0,0)
					colnames(FIT.norm) <- c("chisq","gfi", "agfi", "rmsea", "nfi", "nnfi", "cfi")
					
					} else {
					rep <- rep + 1
    #store the results of parameter estimates 
    beta.norm[i, 1:1] = parameterEstimates (fit.norm, ci = FALSE, standardized = TRUE)[16:16, 4]
    #store the results of model fit 
    FIT.norm[i,1:7] = fitMeasures(fit.norm, c("chisq","GFI", "AGFI", "RMSEA", "NFI", "NNFI","CFI"))
    #store the results of reliability 
    alphaX[i] = alpha(data.norm[1:12])$total[1]
    alphaY[i] = alpha(data.norm[13:15])$total[1]
					}
				}
		} else { 
		break
		}
	}
     ------------------------------------------------------
    #Output

    #Get rid of the non-convergence output
    alphaX <- alphaX[!is.na(alphaX)]
    alphaY <- alphaY[!is.na(alphaY)]
    alphaX <- alphaX[ alphaX != 0 ] 
    alphaY <- alphaY[ alphaY != 0 ]

    alpha_X <- round(colMeans(do.call("rbind",alphaX)),3)
    alpha_Y <- round(colMeans(do.call("rbind",alphaY)),3)
    cbind(alpha_X,alpha_Y)

    count(beta.norm!= 0)
    beta.norm <- beta.norm[!is.na(beta.norm)]
    beta.norm <- beta.norm[ beta.norm!= 0 ] 
    beta=mean(beta.norm) # to extract the means of beta for the 5000 simulations
    sd=sd(t(beta.norm))    #std. deviation 
    se=sd/sqrt(nrep)                 #std. of error
    MSE=round(sum((beta.norm-0.585)^2/nrep),6)
    round(cbind(beta,sd,se,MSE),3)    #Calculate mean in three decimal point

    FIT.norm <- FIT.norm[complete.cases(FIT.norm), ]
    FIT.norm <- FIT.norm[!apply(FIT.norm==0, 1, all),]
    round(colMeans(FIT.norm),3)    #Calculate mean in three decimal point
    
 #### _Non-Normal_
    #Generate item-level data (Moderate Non-normal)
    #sk= 1.00
    #k=1.5
    b<- 0.9530769
    c<- 0.1631943
    d<- 0.0065974
    a<- -c
    --------------------------------------------------------------------
    #Generate item-level data (Severe Non-normal)
    #sk= 1.75
    #k=3.75
    b<- 0.9296605
    c<- 0.3994967
    d<- -0.036467
    a<- -c
    --------------------------------------------------------------------
    data.f=function(n){
   	KSI1=rnorm(n)
  	X1Err=rnorm(n)
	X2Err=rnorm(n)
  	X3Err=rnorm(n)
    X4Err=rnorm(n)
	X5Err=rnorm(n) 
    X6Err=rnorm(n)
	X7Err=rnorm(n)
    X8Err=rnorm(n)
 	X9Err=rnorm(n)
  	X10Err=rnorm(n)
  	X11Err=rnorm(n)
	X12Err=rnorm(n)
	X1 = .70*KSI1 + .714*X1Err 
	X2 = .80*KSI1 + .060*X2Err 
	X3 = .90*KSI1 + .436*X3Err 
	X4 = .70*KSI1 + .714*X4Err 
  	X5 = .80*KSI1 + .060*X5Err 
  	X6 = .90*KSI1 + .436*X6Err 
	X7 = .70*KSI1 + .714*X7Err 
	X8 = .80*KSI1 + .060*X8Err 
  	X9 = .90*KSI1 + .436*X9Err 
  	X10 = .70*KSI1 + .714*X10Err 
	X11 = .80*KSI1 + .060*X11Err 
	X12 = .90*KSI1 + .436*X12Err 

  	NX1 = round(a + b*X1 + c*X1^2 + d*X1^3 + 5.2)
  	NX2 = round(a + b*X2 + c*X2^2 + d*X2^3 + 5.2)
  	NX3 = round(a + b*X3 + c*X3^2 + d*X3^3 + 5.2)
	NX4 = round(a + b*X4 + c*X4^2 + d*X4^3 + 5.2)
  	NX5 = round(a + b*X5 + c*X5^2 + d*X5^3 + 5.2)
	NX6 = round(a + b*X6 + c*X6^2 + d*X6^3 + 5.2)
  	NX7 = round(a + b*X7 + c*X7^2 + d*X7^3 + 5.2)
  	NX8 = round(a + b*X8 + c*X8^2 + d*X8^3 + 5.2)
	NX9 = round(a + b*X9 + c*X9^2 + d*X9^3 + 5.2)
  	NX10 = round(a + b*X10 + c*X10^2 + d*X10^3 + 5.2)
    NX11 = round(a + b*X11 + c*X11^2 + d*X11^3 + 5.2)
    NX12 = round(a + b*X12 + c*X12^2 + d*X12^3 + 5.2)

	  Eta1Err=rnorm(n)
  	  Y1Err=rnorm(n)
	  Y2Err=rnorm(n)
	  Y3Err=rnorm(n)

	  ETA1=.585*KSI1 + .811*Eta1Err
	  Y1=round(.70*ETA1 + .714*Y1Err + 5.2)
	  Y2=round(.80*ETA1 + .060*Y2Err + 5.2)
	  Y3=round(.90*ETA1 + .436*Y3Err + 5.2)

	  data.norm = data.frame(NX1, NX2, NX3, NX4, NX5, NX6, NX7, NX8, NX9, NX10, NX11, NX12, Y1, Y2, Y3)
	
  	return(data.norm)
    }}

    #Fit the item-level model (Non-normal)
    model.nonnormal <- '
      KSI1 =~ NX1 + NX2 + NX3 + NX4 + NX5 + NX6 + NX7 + NX8 + NX9 + NX10 + NX11 + NX12
	   ETA1 =~ Y1 + Y2 + Y3

	   #REGRESSION
	   ETA1 ~ KSI1


## The R-syntax for the simulation of item-parcelling model (2I6P, 3I4P, 4I3P, and 6I2P)
### Generate the 2I6P dataset

    library("lavaan") #required for SEM analysis
    library("matrixcalc") #required for is.positive.definite function
    library("psych")  #required for reliability function

    #Use the generated item-level data 
    
    n=100   #sample size
    nrep=5000  #number of replications
    beta.parcel <- matrix(NA, nrow = nrep, ncol = 1)	#create space to store parameter estimates values
    FIT.parcel <- matrix(NA, nrow = nrep, ncol = 7) #create space to store model fit values
    alphaP <- matrix(NA, nrow = nrep, ncol = 1) #create space to store reliability values
    alphaY <-  matrix(NA, nrow = nrep, ncol = 1) #create space to store reliability values
    
    set.seed(12345)   
    
    for (i in 1:nrep){                    
    data.norm = data.f(n)	
    
    #BODY######################
		parcel1<-data.norm[,c(1:2)]
			Meanparcel1<-rowSums(parcel1/2)
		parcel2<-data.norm[,c(3:4)]
			Meanparcel2<-rowSums(parcel2/2)
		parcel3<-data.norm[,c(5:6)]
			Meanparcel3<-rowSums(parcel3/2)
		parcel4<-data.norm[,c(7:8)]
			Meanparcel4<-rowSums(parcel4/2)
		parcel5<-data.norm[,c(9:10)]
			Meanparcel5<-rowSums(parcel5/2)
		parcel6<-data.norm[,c(11:12)]
			Meanparcel6<-rowSums(parcel6/2)
		Y1 <- data.norm[,13]
		Y2 <- data.norm[,14]
		Y3 <- data.norm[,15]

    data.parcel=data.frame(Meanparcel1,Meanparcel2,Meanparcel3,Meanparcel4,Meanparcel5,Meanparcel6,Y1,Y2,Y3)

    #Fit the 2I6P model
	  	model.norm <- '
			KSI1 =~ Meanparcel1 + Meanparcel2 + Meanparcel3 + Meanparcel4 + Meanparcel5 + Meanparcel6
			ETA1 =~ Y1 + Y2 + Y3
		#REGRESSION
		ETA1 ~ KSI1
			'	
   		fit.parcel <- cfa(model=model.norm, data = data.parcel,estimator="ML")
    ###################################
    
    #store the results of parameter estimates
    beta.parcel[i, 1:1] = parameterEstimates (fit.parcel, ci = FALSE, standardized = TRUE)[10:10, 4] 
    #store the results of model fit
    FIT.parcel[i,1:7] = fitMeasures(fit.parcel, c("chisq","GFI", "AGFI", "RMSEA", "NFI", "NNFI","CFI")) 
    #store the results of reliability 
    alphaP[i] = alpha(data.parcel[1:6])$total[1]
    alphaY[i] = alpha(data.parcel[7:9])$total[1]
    } 
    
### Generate the 3I4P dataset (replace in the BODY)

    parcel1<-data.norm[,c(1:3)]
			Meanparcel1<-rowSums(parcel1/3)
    parcel2<-data.norm[,c(4:6)]
			Meanparcel2<-rowSums(parcel2/3)
	parcel3<-data.norm[,c(7:9)]
			Meanparcel3<-rowSums(parcel3/3)
	parcel4<-data.norm[,c(10:12)]
			Meanparcel4<-rowSums(parcel4/3)
		Y1 <- data.norm[,13]
		Y2 <- data.norm[,14]
		Y3 <- data.norm[,15]

		data.parcel=data.frame(Meanparcel1,Meanparcel2,Meanparcel3,Meanparcel4,Y1,Y2,Y3)

     #Fit the 3I4P model
		  model.norm <- '
			KSI1 =~ Meanparcel1 + Meanparcel2 + Meanparcel3 + Meanparcel4
			ETA1 =~ Y1 + Y2 + Y3
		
		#REGRESSION
		ETA1 ~ KSI1
			'	
 ### Generate the 4I3P dataset (replace in the BODY)
      
    parcel1<-data.norm[,c(1:4)]
			Meanparcel1<-rowSums(parcel1/4)
	parcel2<-data.norm[,c(5:8)]
			Meanparcel2<-rowSums(parcel2/4)
    parcel3<-data.norm[,c(9:12)]
			Meanparcel3<-rowSums(parcel3/4)
		Y1 <- data.norm[,13]
		Y2 <- data.norm[,14]
		Y3 <- data.norm[,15]

		data.parcel=data.frame(Meanparcel1,Meanparcel2,Meanparcel3,Y1,Y2,Y3)

     #Fit the 4I3P model
		model.norm <- '
			KSI1 =~ Meanparcel1 + Meanparcel2 + Meanparcel3 
			ETA1 =~ Y1 + Y2 + Y3
		
		#REGRESSION
		ETA1 ~ KSI1
			'	
### Generate the 6I2P dataset (replace in the BODY)
    
    parcel1<-data.norm[,c(1:6)]
			Meanparcel1<-rowSums(parcel1/6)
	parcel2<-data.norm[,c(7:12)]
			Meanparcel2<-rowSums(parcel2/6)
		Y1 <- data.norm[,13]
		Y2 <- data.norm[,14]
		Y3 <- data.norm[,15]

		data.parcel=data.frame(Meanparcel1,Meanparcel2,Y1,Y2,Y3)
  
     #Fit the 6I2P model
		  model.norm <- '
			KSI1 =~ Meanparcel1 + Meanparcel2
			ETA1 =~ Y1 + Y2 + Y3
		
		#REGRESSION
		ETA1 ~ KSI1
			'	

      
