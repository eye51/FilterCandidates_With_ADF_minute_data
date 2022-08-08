#
#   Read Index and scan for stationarity.

##############################################################################################
# First define the ADF function. 
##############################################################################################

BastiaanADF <-function (x, type = c("nc", "c", "ct"), lags = 1)

# Based on the ADF function of the Fseries package by Diethelm Wuertz

{
    if (ncol(as.matrix(x)) > 1) 
        stop("x is not a vector or univariate time series")
    if (any(is.na(x))) 
        stop("NAs in x")
    if (lags < 0) 
        stop("lags negative")
    doprint = FALSE
    CALL = match.call()
    DNAME = deparse(substitute(x))
    type = type[1]
    x.name = deparse(substitute(x))
    lags = lags + 1
    y = diff(x)
    n = length(y)
    z = embed(y, lags)
    y.diff = z[, 1]
    y.lag.1 = x[lags:n]
    tt = lags:n
    if (lags > 1) {
        y.diff.lag = z[, 2:lags]
        if (type == "nc") {
            res = lm(y.diff ~ y.lag.1 - 1 + y.diff.lag)
        }
        if (type == "c") {
            res = lm(y.diff ~ y.lag.1 + 1 + y.diff.lag)
        }
        if (type == "ct") {
            res = lm(y.diff ~ y.lag.1 + 1 + tt + y.diff.lag)
        }
    }
    else {
        if (type == "nc") {
            res = lm(y.diff ~ y.lag.1 - 1)
        }
        if (type == "c") {
            res = lm(y.diff ~ y.lag.1 + 1)
        }
        if (type == "ct") {
            res = lm(y.diff ~ y.lag.1 + 1 + tt)
        }
    }
    res.sum = summary(res)
    if (doprint) 
        print(res.sum)
    if (type == "nc") 
        coefNum = 1
    else coefNum = 2
    STAT = res.sum$coefficients[coefNum, 1]/res.sum$coefficients[coefNum, 
        2]
    if (type == "nc") 
        table = cbind(c(-2.66, -2.26, -1.95, -1.6, +0.92, +1.33, 
            +1.7, +2.16), c(-2.62, -2.25, -1.95, -1.61, +0.91, 
            +1.31, +1.66, +2.08), c(-2.6, -2.24, -1.95, -1.61, 
            +0.9, +1.29, +1.64, +2.03), c(-2.58, -2.23, -1.95, 
            -1.62, +0.89, +1.29, +1.63, +2.01), c(-2.58, -2.23, 
            -1.95, -1.62, +0.89, +1.28, +1.62, +2), c(-2.58, 
            -2.23, -1.95, -1.62, +0.89, +1.28, +1.62, +2))
    if (type == "c") 
        table = cbind(c(-3.75, -3.33, -3, -2.63, -0.37, +0, +0.34, 
            +0.72), c(-3.58, -3.22, -2.93, -2.6, -0.4, -0.03, 
            +0.29, +0.66), c(-3.51, -3.17, -2.89, -2.58, -0.42, 
            -0.05, +0.26, +0.63), c(-3.46, -3.14, -2.88, -2.57, 
            -0.42, -0.06, +0.24, +0.62), c(-3.44, -3.13, -2.87, 
            -2.57, -0.43, -0.07, +0.24, +0.61), c(-3.43, -3.12, 
            -2.86, -2.57, -0.44, -0.07, +0.23, +0.6))
    if (type == "ct") 
        table = cbind(c(-4.38, -3.95, -3.6, -3.24, -1.14, -0.8, 
            -0.5, -0.15), c(-4.15, -3.8, -3.5, -3.18, -1.19, 
            -0.87, -0.58, -0.24), c(-4.04, -3.73, -3.45, -3.15, 
            -1.22, -0.9, -0.62, -0.28), c(-3.99, -3.69, -3.43, 
            -3.13, -1.23, -0.92, -0.64, -0.31), c(-3.98, -3.68, 
            -3.42, -3.13, -1.24, -0.93, -0.65, -0.32), c(-3.96, 
            -3.66, -3.41, -3.12, -1.25, -0.94, -0.66, -0.33))
    table = t(table)
    tablen = dim(table)[2]
    tableT = c(25, 50, 100, 250, 500, 1e+05)
    tablep = c(0.01, 0.025, 0.05, 0.1, 0.9, 0.95, 0.975, 0.99)
    tableipl = numeric(tablen)
    for (i in (1:tablen)) tableipl[i] = approx(tableT, table[, 
        i], n, rule = 2)$y
    PVAL = approx(tableipl, tablep, STAT, rule = 2)$y
    if (is.na(approx(tableipl, tablep, STAT, rule = 1)$y)) {
        if (PVAL == min(tablep)) {
            warning("p-value smaller than printed p-value")
        }
        else {
            warning("p-value greater than printed p-value")
        }
    }
    PARAMETER = lags - 1
    names(PARAMETER) = "Lag order"
    METHOD = "Augmented Dickey-Fuller Test"
    names(STAT) = "Dickey-Fuller"
    test = list(statistic = STAT, parameter = PARAMETER, p.value = PVAL, 
        method = METHOD, data.name = DNAME)
    class(test) = c("list", "htest")
    title = test$method
    description = date()
    new("fURTEST", call = CALL, data = as.data.frame(x), 
        test = test, title = as.character(title), description = as.character(description))
}

##############################################################################################
# Combine 3 time series and calculate ADF. 
##############################################################################################


BastiaanMin3Series <-function (start,timeseries)
{


	y<-as.vector(timeseries)
	
	NoElements<-length(y)/3


	x1<-y[1:NoElements]
	x2<-y[(NoElements+1):(2*NoElements)]
	x3<-y[(2*NoElements+1):(3*NoElements)]

	a<-start[1]
	b<-start[2]
	

	c<- -a-b
	
	noemer<-numeric(NoElements)
	teller<-numeric(NoElements)

	noemer<-0
	teller<-0

	if (a < 0) noemer=noemer-a*x1 else teller=teller+a*x1
	if (b < 0) noemer=noemer-b*x2 else teller=teller+b*x2
	if (c < 0) noemer=noemer-c*x3 else teller=teller+c*x3


	e<-noemer/teller

	hh<-BastiaanADF(e,"c",1)

	stat<-hh@test$statistic

	return (stat)
}

##############################################################################################
# Combine 3 time series and calculate ADF. 
##############################################################################################


BastiaanMin3Series2 <-function (start,timeseries)
{


	y<-as.vector(timeseries)
	
	NoElements<-length(y)/3


	x1<-y[1:NoElements]
	x2<-y[(NoElements+1):(2*NoElements)]
	x3<-y[(2*NoElements+1):(3*NoElements)]

	a<-start[1]
	b<-start[2]
	

	c<- -a-b
	
	noemer<-numeric(NoElements)
	teller<-numeric(NoElements)

	noemer<-0
	teller<-0

	if (a < 0) noemer=noemer-a*x1 else teller=teller+a*x1
	if (b < 0) noemer=noemer-b*x2 else teller=teller+b*x2
	if (c < 0) noemer=noemer-c*x3 else teller=teller+c*x3


	e<-teller/noemer

	hh<-BastiaanADF(e,"c",1)

	stat<-hh@test$statistic

	return (stat)
}

##############################################################################################
#				End of own defined functions
##############################################################################################

#
# Read file containing all the index members:




filenameAllCandidates <- paste("C:/R-programms/Bastiaan_projects/FilterCandidates_with_ADF/candidates/all_candidates.txt",sep="")
filenameAllCandidates2 <- paste("C:/R-programms/Bastiaan_projects/FilterCandidates_with_ADF/candidates/all_candidates2.txt",sep="")
filenameOverzichCandidates <- paste("C:/R-programms/Bastiaan_projects/FilterCandidates_with_ADF/candidates/overzicht_candidates.txt",sep="")
filenameOverzichCandidates2 <- paste("C:/R-programms/Bastiaan_projects/FilterCandidates_with_ADF/candidates/overzicht_candidates2.txt",sep="")


cat ("",file=filenameAllCandidates,sep="",append=FALSE)
cat ("",file=filenameOverzichCandidates,sep="",append=FALSE)

cIndexNames <- scan (file="C:/R-programms/Bastiaan_projects/FilterCandidates_with_ADF/Russell1000_22-12-2004/Index.txt",what=character(0),sep="\n",quiet=TRUE)

counter<-0
iNumberOfStocks=(length(cIndexNames)-1)

print (iNumberOfStocks)
print (iNumberOfStocks)

#iNumberOfStock=10

for (i in 1:(iNumberOfStocks-2))


{

    print (i)
    print (iNumberOfStocks)

       stock1 <- paste(cIndexNames[i])
       filename <- paste("C:/R-programms/Bastiaan_projects/FilterCandidates_with_ADF/Russell1000_22-12-2004/",cIndexNames[i],".txt",sep="")

                                                                                                                                                                        
        x1 <- scan(file=filename,quiet=TRUE)   
        NumberObs1=length(x1)-4
        IC1=x1[3]                                            # Bloomberg Indusry code
        Capitalization1=x1[5]                                # cap. of company, use only > 500
        AvgStocks1=x1[6]       

        print (stock1)
        StockData1<- x1[7:length(x1)]                              # time-series data starts at 5th position
        LogStocData1<-log(StockData1)                              # use rendement for determining stationarity
        
        verhouding1=abs(max(StockData1)/min(StockData1))

        if ((Capitalization1 >= 2500) && (AvgStocks1>200000))
        {
        
            for (j in (i+1):(iNumberOfStocks-1))
             
            
            {
		
	       stock2 <- paste(cIndexNames[j])
             filename2 <- paste("C:/R-programms/Bastiaan_projects/FilterCandidates_with_ADF/Russell1000_22-12-2004/",cIndexNames[j],".txt",sep="")

                x2 <- scan(file=filename2,quiet=TRUE)                                                                 
                
                NumberObs2=length(x2)-4
                IC2=x2[3]                                         # Bloomberg Indusry code
                Capitalization2=x2[5]                                # cap. of company, use only > 2500
                AvgStocks2=x2[6]       
               
                StockData2<- x2[7:length(x2)]                              # time-series data starts at 5th position
                LogStocData2<-log(StockData2)                               # use rendement for determining stationarity
                verhouding2=abs(max(StockData2)/min(StockData2))
                
                RatioVerhouding=verhouding1/verhouding2;
                
                if (RatioVerhouding < 1)
                
                {
                
                    RatioVerhouding=verhouding2/verhouding1;
                
                }
                
                
                if ((Capitalization2 >=2500) && (NumberObs1==NumberObs2) && (AvgStocks2>200000)
                    && (IC1==IC2) && (RatioVerhouding < 3) && (IC1!=52593))
                {
                
      		     for (k in (j+1):iNumberOfStocks)
             
            
           			 {
			       stock3 <- paste(cIndexNames[k])
           			filename3 <- paste("C:/R-programms/Bastiaan_projects/FilterCandidates_with_ADF/Russell1000_22-12-2004/",cIndexNames[k],".txt",sep="")

                		x3 <- scan(file=filename3,quiet=TRUE)                                                                 
                
                		NumberObs3=length(x3)-4
                		IC3=x3[3]                                         # Bloomberg Indusry code
                		Capitalization3=x3[5]                                # cap. of company, use only > 2500
                		AvgStocks3=x3[6]       
               
                		StockData3<- x3[7:length(x3)]                              # time-series data starts at 5th position
                		LogStocData3<-log(StockData3)                               # use rendement for determining stationarity
                		verhouding3=abs(max(StockData3)/min(StockData3))
                
                		RatioVerhouding1=verhouding1/verhouding3;
                		RatioVerhouding2=verhouding2/verhouding3;
                
                		if (RatioVerhouding1 < 1)
                
                		{
                
              		      RatioVerhouding1=verhouding3/verhouding1;
                
                		}
                		if (RatioVerhouding2 < 1)
                
                		{
                
              		      RatioVerhouding2=verhouding3/verhouding2;
                
                		}
                
                
	                if ((Capitalization3 >=2500) && (NumberObs1==NumberObs3) && (AvgStocks3>200000)
      		              && (IC1==IC3) && (RatioVerhouding < 3) && (IC1!=52593))


				{


					combi=cbind(StockData1,StockData2,StockData3)
					start<-numeric(2)
					start[1]=0.5
					start[2]=0.5


					suppressWarnings(min1<-nlm(BastiaanMin3Series,start,timeseries=combi))
					suppressWarnings(min2<-nlm(BastiaanMin3Series2,start,timeseries=combi))
               
					crit<- -4.3
					stat1<-min1$minimum
					stat2<-min2$minimum


					if ((stat1 <crit) && (stat2 <crit))
					{

					print (stock1)
					print (stock2)
					print (stock3)



					noemer1<-numeric(NumberObs1)
					teller1<-numeric(NumberObs1)

					noemer1<-0
					teller1<-0

					a<-min1$estimate[1]
					b<-min1$estimate[2]
					c<- -a-b

					if (a < 0) noemer1=noemer1-a*StockData1 else teller1=teller1+a*StockData1
					if (b < 0) noemer1=noemer1-b*StockData2 else teller1=teller1+b*StockData2
					if (c < 0) noemer1=noemer1-c*StockData3 else teller1=teller1+c*StockData3


					testData1<-teller1/noemer1

					noemer2<-numeric(NumberObs1)
					teller2<-numeric(NumberObs1)

					noemer2<-0
					teller2<-0
					a<-min2$estimate[1]
					b<-min2$estimate[2]
					c<- -a-b


					if (a < 0) noemer2=noemer2-a*StockData1 else teller2=teller2+a*StockData1
					if (b < 0) noemer2=noemer2-b*StockData2 else teller2=teller2+b*StockData2
					if (c < 0) noemer2=noemer2-c*StockData3 else teller2=teller2+c*StockData3


					testData2<-teller2/noemer2

                            	print ("all 3 stable")
                            	print 

                           	candidate1=strsplit (cIndexNames[i]," US EQUITY")
	                        candidate2=strsplit (cIndexNames[j]," US EQUITY")
	                        candidate3=strsplit (cIndexNames[k]," US EQUITY")
                            
					print (candidate1)
					print (candidate2)
					print (candidate3)

####						Write time series in single file			####

                           
                        	    filenameCandidate=paste("C:/R-programms/Bastiaan_projects/FilterCandidates_with_ADF/candidates/",candidate1[[1]],"_",candidate2[[1]],"_",candidate3[[1]],"1.txt",sep="")
                    	          cat (testData1,file=filenameCandidate,sep=",")

                        	    filenameCandidate=paste("C:/R-programms/Bastiaan_projects/FilterCandidates_with_ADF/candidates/",candidate1[[1]],"_",candidate2[[1]],"_",candidate3[[1]],"2.txt",sep="")
                    	          cat (testData2,file=filenameCandidate,sep=",")
                            

####						Write time series in overview files			####


            	                cat (candidate1[[1]],"_",candidate2[[1]],"_",candidate3[[1]],",",file=filenameAllCandidates,sep="",append=TRUE)
      	                      cat (testData1,file=filenameAllCandidates,sep=",",append=TRUE)
	                            cat ("\n",file=filenameAllCandidates,sep="",append=TRUE)

            	                cat (candidate1[[1]],"_",candidate2[[1]],"_",candidate3[[1]],",",file=filenameAllCandidates2,sep="",append=TRUE)
      	                      cat (testData2,file=filenameAllCandidates2,sep=",",append=TRUE)
	                            cat ("\n",file=filenameAllCandidates2,sep="",append=TRUE)



####						Write combination in overview files			####

                            
                  	          cat ("1st=,",file=filenameOverzichCandidates,sep="",append=TRUE)
            	                cat (min1$estimate[1],file=filenameOverzichCandidates,sep=",",append=TRUE)
      	                      cat (",2nd=,",file=filenameOverzichCandidates,sep="",append=TRUE)
  	                            cat (min1$estimate[2],file=filenameOverzichCandidates,sep=",",append=TRUE)
					    cat (",3rd=,",file=filenameOverzichCandidates,sep="",append=TRUE)
	            	          cat (-min1$estimate[1]-min1$estimate[2],file=filenameOverzichCandidates,sep=",",append=TRUE)
	                            cat (",",candidate1[[1]],"_",candidate2[[1]],"_",candidate3[[1]],",",file=filenameOverzichCandidates,sep="",append=TRUE)
	                            cat (IC1,"\n",file=filenameOverzichCandidates,sep="",append=TRUE)


                            
                  	          cat ("1st=,",file=filenameOverzichCandidates2,sep="",append=TRUE)
            	                cat (min2$estimate[1],file=filenameOverzichCandidates2,sep=",",append=TRUE)
      	                      cat (",2nd=,",file=filenameOverzichCandidates2,sep="",append=TRUE)
  	                            cat (min2$estimate[2],file=filenameOverzichCandidates2,sep=",",append=TRUE)
					    cat (",3rd=,",file=filenameOverzichCandidates2,sep="",append=TRUE)
	            	          cat (-min2$estimate[1]-min2$estimate[2],file=filenameOverzichCandidates2,sep=",",append=TRUE)
	                            cat (",",candidate1[[1]],"_",candidate2[[1]],"_",candidate3[[1]],",",file=filenameOverzichCandidates2,sep="",append=TRUE)
	                            cat (IC1,"\n",file=filenameOverzichCandidates2,sep="",append=TRUE)


                            }
                    
			}

	
			}

    
                }                                              
                          
                
                }                                              # end of IF : both stocks satisfy intial selection criteria
                
            }                                                  # loop over secundary stock

        }                                                      # stock 1 satisfies capitalization criteria
}
