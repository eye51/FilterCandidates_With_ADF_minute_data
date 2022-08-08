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
#				End of own defined functions
##############################################################################################

#
# Read file containing all the index members:


HoursTraded=7		# data downloaded from 9:30 till 17:00
MinutesTraded=30
TickSize=1			# 1 minute data

NumberOfTicks=(HoursTraded*60/TickSize)+(MinutesTraded/TickSize)+1

print (NumberOfTicks)

filenameAllCandidates <- paste("C:/R-programms/Bastiaan_projects/FilterCandidates_With_ADF_minute_data/candidates/all_candidates.txt",sep="")
filenameOverzichCandidates <- paste("C:/R-programms/Bastiaan_projects/FilterCandidates_With_ADF_minute_data/candidates/overzicht_candidates.txt",sep="")

cat ("",file=filenameAllCandidates,sep="",append=FALSE)
cat ("",file=filenameOverzichCandidates,sep="",append=FALSE)

cIndexNames <- scan (file="C:/R-programms/Bastiaan_projects/FilterCandidates_With_ADF_minute_data/dax/Index.txt",what=character(0),sep="\n",quiet=TRUE)

counter<-0
iNumberOfStocks=(length(cIndexNames)-1)

print (iNumberOfStocks)

for (i in 1:(iNumberOfStocks-1))


{

    print (i)
       stock1 <- paste(cIndexNames[i])
       filename <- paste("C:/R-programms/Bastiaan_projects/FilterCandidates_With_ADF_minute_data/dax/",cIndexNames[i],".txt",sep="")
                          
	input1<-read.table (filename ,sep=",",header=FALSE) # format date,best bid, best ask, traded

#       x1 <- scan(file=filename,quiet=TRUE)   

	x1<-input1[[4]]	# traded price

        NumberObs1=length(x1)
	  NumberOfDays1=NumberObs1/NumberOfTicks

        print (stock1)
        StockData1<- x1[1:length(x1)]                              # time-series data starts at 5th position
        LogStocData1<-log(StockData1)                              # use rendement for determining stationarity
        
        verhouding1=abs(max(StockData1)/min(StockData1))

                
        for (j in (i+1):iNumberOfStocks)
             
            
            {

             filename2 <- paste("C:/R-programms/Bastiaan_projects/FilterCandidates_With_ADF_minute_data/dax/",cIndexNames[j],".txt",sep="")

	 	input2<-read.table (filename2 ,sep=",",header=FALSE) # format date,best bid, best ask, traded
		x2<-input2[[4]]
#               x2 <- scan(file=filename2,quiet=TRUE)                                                                 

                
                NumberObs2=length(x2)
              
                StockData2<- x2[1:length(x2)]                              # time-series data starts at 5th position
                LogStocData2<-log(StockData2)                               # use rendement for determining stationarity
                
			 if ((NumberObs1==NumberObs2))	# time series of the same size -> Start testing


			{
                
                
			TestStatistic<-numeric(NumberOfDays1) # define array for ADF results
			NoDaysStable<-0

			for (k in 1:NumberOfDays1)
			{
			
	                testData=LogStocData1[((k-1)*NumberOfDays1+1):(k*NumberOfDays1+1)]-LogStocData2[((k-1)*NumberOfDays1+1):(k*NumberOfDays1+1)]
                
      	          suppressWarnings(hh<-BastiaanADF (testData,"c",1))
                
          	      crit<- -3.3
          	      TestStatistic[k]<-hh@test$statistic

     		     	      if (TestStatistic[k]< crit)
          		      {
					NoDaysStable<-NoDaysStable+1

				}


			}	# loop over all days

			print (paste(cIndexNames[i],cIndexNames[j],NoDaysStable))


			if (NoDaysStable>(0.25*NumberOfDays1))

                      {                                   # satify all stability criteria..
                            print ("all 3 stable")
                             candidate1=strsplit (cIndexNames[i]," NA EQUITY")
                            candidate2=strsplit (cIndexNames[j]," NA EQUITY")
                            
                            print (candidate1)
                            print (candidate2)
					print (NoDaysStable)
					print ((0.75*NumberOfDays1))
	 	                testData=LogStocData1[1:NumberObs1]-LogStocData2[1:NumberObs1]
                          
                            filenameCandidate=paste("C:/R-programms/Bastiaan_projects/FilterCandidates_With_ADF_minute_data/candidates/",candidate1[[1]],"_",candidate2[[1]],".txt",sep="")
                            cat (testData,file=filenameCandidate,sep=",")
                            

                            cat (candidate1[[1]],"_",candidate2[[1]],",",file=filenameAllCandidates,sep="",append=TRUE)
                            cat (testData,file=filenameAllCandidates,sep=",",append=TRUE)
                            cat ("\n",file=filenameAllCandidates,sep="",append=TRUE)
                            
                            
                            
                       }	# satify all stability criteria..

                
            }                                                  # loop over secundary stock

        }                                                      # stock 1 satisfies capitalization criteria
}
