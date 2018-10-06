#' A function to find two optimal cutpoints
#'
#' @title  Find two optimal cutpoints using optimal equal-HR method
#' @description Use optimal equal-HR method to determine two optimal cutpoints of a continuous predictor that has a U-shape relationship with survival outcomes based on Cox regression model.
#' @param  cox_pspline_fit  Cox model with psplined x, e.g. coxph(Surv(t,d)~pspline(x,df=0,caic=T),data=test).
#' @param  data a dataframe contain survival outcome and a continuous variable which needs to find two optimal cutpoints.
#' @param  nquantile an integer; the default value is 100, which means using the 100-quantiles of log relative hazard to find cutpoints.
#' @param  exclude a decimals; it is used for excluding extreme values of log relative hazardthe. The default value is 0.05, which log relative hazard values smaller than 5th percentile or larger than 95th percentile are excluded.
#' @param  eps a decimals; the default value is 0,01. It restrict the difference between the log relative hazard values of two cadidate cutpoints to be less than 0.01.
#' @param  shape a string; equals "U" or "inverseU"
#' @export
#' @examples
#'
#' \donttest{
#' #### Example 1. Find two optimal cutpoints in an univariate Cox model
#' # Fit an univariate Cox model with pspline
#' require(survival)
#' result <- coxph(Surv(t,d)~pspline(x,df=0,caic=TRUE),data=test)
#' # Visualize the relationship
#' # Explore whether there is a U-shaped relationship between x and log relative hazard
#' termplot(result,se=TRUE,col.term=1,ylab='log relative hazard')
#' # Find two opitmal cutpoints using optimal equal-HR method.
#' cuts <- findcutpoints(cox_pspline_fit = result, data = test, shape='U')
#' cuts$optimal # output two optimal cutpoints
#'
#' #### Example 2. Find two optimal cutpoints in a multivariate Cox model
#' # Fit a multivariate Cox model with pspline
#' # The independent variable which is need to find cutpoints should be placed before other covariates.
#' # To find cutpoints of x, Surv(t,d)~pspline(x)+x1 should be used instead of Surv(t,d)~x1+pspline(x)
#' require(survival)
#' result <- coxph(Surv(t,d)~pspline(x,df=0,caic=TRUE)+x1,data=test)
#' # The rest procedure is the same as example 1
#' # Visualize the relationship
#' # Explore whether there is a U-shaped relationship between x and log relative hazard
#' termplot(result,se=TRUE,col.term=1,ylab='log relative hazard')
#' # Find two opitmal cutpoints of the first independent variable.
#' cuts <- findcutpoints(cox_pspline_fit = result, data = test, shape='U')
#' cuts$optimal # output two optimal cutpoints
#' }

findcutpoints <- function(cox_pspline_fit,data,nquantile=100,exclude=0.05,eps=0.01, shape="U"){
  #For R CMD check
  x <- NULL

  #Check the parameters
  if(missing(cox_pspline_fit)|!('coxph'%in%class(cox_pspline_fit))){
    stop('"cox_pspline_fit" is missing or incorrect')
  }
  if(missing(data)|class(data)!="data.frame"){
    stop('"data" is missing or incorrect')
  }
  if(shape!="U" & shape !="inverseU" ){stop('Invalid value for "shape" parameter')}

  #rename the variable names in dataset, (t,d,x,x1,x2,...)
  f <- cox_pspline_fit$formula
  variablenames <- all.vars(f)
  data <- data[,variablenames]
  if(length(variablenames)<3){
    stop("There are not enough variables in dataset. At least (times, cersor, x) are required.")
  } else{
    colnames(data) <- c('t','d','x')
    confounders <- NA
    if (length(variablenames)>3){
      for (i in 4:length(variablenames)){
        colnames(data)[i] <- paste("x", i-3, sep="")
      }
      confounders <- colnames(data)[4:length(variablenames)]
    }
  }


  #adjpiont:
  #find one point equal to ycut[i]
  #find two points(smaller & larger than ycut[i]) to determine the pseudo cutpoint
  #None of both: no cutpoints
  adjpoint <- function(data,y){
    lesslarger <- NA
    if (dim(data)[1]<=1){
      pseudocut  <- NA
    }else{
      for ( j in c(1:(dim(data)[1]-1))){
        if((data[j,'y']<y & data[j+1,'y']> y) |
           (data[j,'y']>y & data[j+1,'y']< y)){
          lesslarger <- j
        }
      }
      equal <- NA
      for ( j in c(1:dim(data)[1])){
        if(data[j,'y']== y){
          equal <- j
        }
      }
      pseudocut <- NA
      if (!is.na(equal)){
        pseudocut  <- data[equal,'x']
      }else{
        if(!is.na(lesslarger)){
          pseudocut  <- (data[lesslarger+1,'x']-data[lesslarger,'x'])*(y-data[lesslarger,'y'])/(data[lesslarger+1,'y']-data[lesslarger,'y'])+data[lesslarger,'x']
        }else{
          pseudocut  <- NA
        }
      }
    }
    return(pseudocut)
  }

  #get x and estimated y
  ptemp <- termplot(cox_pspline_fit, se=TRUE, plot=FALSE)
  xterm  <- ptemp[[1]]
  PI_fit <- xterm$y

  #Define the minimum y value as turning point of U shape , maximum y value as turning point of inverse U shape
  turningpoint_index <- NA
  if (shape=="U"){
    turningpoint_index <- which(PI_fit == min(PI_fit))
  }else{
    if(shape=="inverseU"){
      turningpoint_index <- which(PI_fit == max(PI_fit))
    }
  }


  #quantiles of y after excluding some extreme values
  ycut <- quantile(PI_fit, probs = seq((0+exclude), (1-exclude),1/nquantile))

  #cut_spline : a matrix contains cutpoints results
  cut_spline <-  matrix(NA,length(ycut),7)
  colnames(cut_spline)=c('Quantile','Lindex','Rindex','Lcutpoint','Rcutpoint','AIC','#medianrange')

  # find all candidate pairs of cutpoints
  for (i in 1:length(ycut)){
    y <- as.numeric(ycut[i])
    cut_spline[i,"Quantile"] <- names(ycut)[i]
    L <- which.min((PI_fit[0:turningpoint_index]-y)^2)
    R <- which.min((PI_fit[turningpoint_index:length(PI_fit)]-y)^2) + turningpoint_index - 1
    #control that two cutpoints have approximate estimated y value
    if(abs(PI_fit[L]-PI_fit[R])>eps){
      #sort by x, split by turningpoint_index
      adjust   <- xterm[order(xterm$x),c('x','y')]
      split    <- which(adjust$y == min(PI_fit))
      adjust_L <- adjust[1:split,]
      adjust_R <- adjust[split:dim(adjust)[1],]

      cut_spline[i,"Lcutpoint"] <- round(adjpoint(adjust_L,y),4)
      cut_spline[i,"Rcutpoint"] <- round(adjpoint(adjust_R,y),4)

    } else {
      cut_spline[i,"Lindex"] <- L
      cut_spline[i,"Rindex"] <- R

      cut_spline[i,"Lcutpoint"] <- xterm$x[L]
      cut_spline[i,"Rcutpoint"] <- xterm$x[R]
    }

    # fit Cox model with discrete x_c
    xl <- as.numeric(cut_spline[i,"Lcutpoint"])
    xr <- as.numeric(cut_spline[i,"Rcutpoint"])

    #
    if(!is.na(cut_spline[i,"Lcutpoint"]) & !is.na(cut_spline[i,"Rcutpoint"])){
      datatemp <- data
      datatemp <- within(datatemp,{
        x_c <- NA
        x_c[x < xl] <- 'C1'
        x_c[x >= xl & x <= xr] <- 'C2'
        x_c[x > xr] <- 'C3'
      }
      )
      datatemp$x_c <- relevel(factor(datatemp$x_c),ref='C2')

      #Cox model formula
      if(sum(is.na(confounders))){
        f_cox <- Surv(t,d)~factor(x_c)
      }else{
        f_cox <- as.formula(paste('Surv(t,d)~factor(x_c)','+', paste(confounders, collapse=" + ")))
      }

      coxHR <- survival::coxph(f_cox ,method='breslow',data=datatemp)
      cut_spline[i,"AIC"] <- as.numeric(extractAIC(coxHR)[2])
      cut_spline[i,'#medianrange'] <- sum(datatemp$x_c=='C2',na.rm = T)
    }
  }

  #multiple pairs of optimal cutpoints: use the median value
  number <-  which(cut_spline[,"AIC"]==min(cut_spline[,"AIC"], na.rm = T))
  multi  <-  length(number)
  if (multi>1){
    if(multi%%2 ==0){
      optimal <- c('Cutpoint_L'= (as.numeric(cut_spline[number[ceiling(median(multi))],"Lcutpoint"])+ as.numeric(cut_spline[number[floor(median(multi))],"Lcutpoint"]))/2,
                   'Cutpoint_R'= (as.numeric(cut_spline[number[ceiling(median(multi))],"Rcutpoint"])+ as.numeric(cut_spline[number[floor(median(multi))],"Rcutpoint"]))/2)
    } else {
      optimal <- c('Cutpoint_L'= as.numeric(cut_spline[number[median(multi)],"Lcutpoint"]),
                   'Cutpoint_R'= as.numeric(cut_spline[number[median(multi)],"Rcutpoint"]))
    }
  }else{
    optimal <- c('Cutpoint_L'= as.numeric(cut_spline[number,"Lcutpoint"]),'Cutpoint_R'=as.numeric(cut_spline[number,"Rcutpoint"]))
  }

  #return results
  cut_points<- list('allcuts'=cut_spline,'optimal'=optimal)
  return(cut_points)
}
