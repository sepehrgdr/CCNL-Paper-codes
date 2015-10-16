# CCNL-Paper-codes
# Continuous Cross Nested Logit

#Define Dataframe
nobs<-1000
df <- data.frame(ID=1:nobs, stringsAsFactors=FALSE)

#OD distance is Gamma distributed
dist <- rgamma(nobs,2,0.25)

df$distance<-dist

#OD Freeflow Speed is Discrete Uniform distributed
fspeed <- sample(c(40,50,60,70,80,90),nobs,replace=T)

df$ODfspeed<-fspeed

#OD Peak Speed is Discrete Uniform distributed
pspeed <- sample(c(15,20,25,30,35,40),nobs,replace=T)

df$ODpspeed<-pspeed

#Delay
delay <- 1 - (pspeed/fspeed)

df$delay<-delay

#Continuous Time and Cost Function
#Function defining travel speed function, use similar formula from paper: Time of Day Modeling in a Tour-Based Context: The Tel-Aviv Experience 
defineTravelTimeFunction<- function(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16,a17,a18,a19){
    fTime <- function(x, id){
        #        a*x^4+b*x^3+c*x^2+d*x+e
        re <- x
        delay<-df$delay[id]
        distance<-df$distance[id]
        fspeed<-df$ODfspeed[id]
        for (i in 1:length(x)){
            re[i]=distance/(fspeed*exp(a1+a2*log(distance)+a3*delay+a4*delay*exp(sin(x[i]*2*pi/24))+a5*delay*exp((sin(x[i]*2*pi/24)))^2+a6*delay*exp((sin(x[i]*2*pi/24)))^3+a7*delay*exp((sin(x[i]*2*pi/24)))^4
                                       + a8*delay*exp(cos(x[i]*2*pi/24))+a9*delay*exp((cos(x[i]*2*pi/24)))^2+a10*delay*exp((cos(x[i]*2*pi/24)))^3+a11*delay*exp((cos(x[i]*2*pi/24)))^4
                                       +a12*delay*exp(sin(x[i]*4*pi/24))+a13*delay*exp((sin(x[i]*4*pi/24)))^2+a14*delay*exp((sin(x[i]*4*pi/24)))^3+a15*delay*exp((sin(x[i]*4*pi/24)))^4
                                       + a16*delay*exp(cos(x[i]*4*pi/24))+a17*delay*exp((cos(x[i]*4*pi/24)))^2+a18*delay*exp((cos(x[i]*4*pi/24)))^3+a19*delay*exp((cos(x[i]*4*pi/24)))^4))
        }
        re
    }
    fTime
}


#Step function
defineTravelCostFunction<- function(a1,a2){
    fCost <- function(x, id){
        #        a*x^4+b*x^3+c*x^2+d*x+e
        re <- x
        for (i in 1:length(x)){
            re[i] <- 0
            if (x[i] > 5 & x[i] < 10) {
                re[i] <- a1
            } else if ( x[i] > 15 & x[i] < 22) {
                re[i] <- a2
            }
            
        }
        re
    }
    fCost
}


#Assume function for travel time and TTR
# networktraveltime_mnl <- defineTravelTimeMatrix(0.0317,0.0072,-1.4865,-11.2967,12.5130,-5.6582,0.9182,7.0780,-5.6571,1.8769,-0.2359,-1.9606,1.7074,-0.5433,0.0630,2.3598,-0.7701,0.1102,0.0018)
networktraveltime_cl<-defineTravelTimeFunction(0.0317,0.0072,-1.4865,-11.2967,12.5130,-5.6582,0.9182,7.0780,-5.6571,1.8769,-0.2359,-1.9606,1.7074,-0.5433,0.0630,2.3598,-0.7701,0.1102,0.0018)
# networktravelcost_mnl<-defineTravelCostMatrix(1,2)
networktravelcost_cl<-defineTravelCostFunction(1,2)
#networktravelcost_cl(3,1)
# networktravelcost_mnl(base)

curve(networktraveltime_cl(x,100),from=0,to=24)

#CCNL - True Model
library(stats4)

#CCNL Density Function
CCNL.density_plot<-function(t,L,U,ID,beta_h,beta_pho,a,b,c,d,e,ff,g,h,i,j,k) {
    
  #function to calculate allocation parameter alpha
    alphaCal <- function(t1, t2, beta_h){
        res = t1
        index = (abs(t1-t2)<=beta_h)
        index1 = (abs(t1-t2)>beta_h)
        res[index] = (beta_h-abs(t1[index]-t2[index]))/(beta_h^2)
        res[index1] = 0
        res = Vectorize(res)
        res
    }
    
    #function to calculate utility function
    yCal <- function(t,ID){
        t[t<0] = 24+t[t<0]
        t[t>24] = t[t>0]-24
        V_t <- a + b*60*networktraveltime_cl(t,ID)+c*networktravelcost_cl(t,ID)+d*sin(t*2*pi/24)+e*sin(t*4*pi/24)+ff*sin(t*6*pi/24)+g*sin(t*8*pi/24)+h*cos(t*2*pi/24)+i*cos(t*4*pi/24)+j*cos(t*6*pi/24)+k*cos(t*8*pi/24)
        
    }
    #function to calculate time-nest utility
    tNetCal <- function(t,t_net,ID,beta_pho,beta_h){
        alpha_t = alphaCal(t,t_net,beta_h)
        y_t = yCal(t,ID)
        (alpha_t*exp(y_t))^beta_pho
    }
    
    tNetCal_vec = Vectorize(tNetCal)
    
    #function to calculate utility of tk in nest tnet
    utilityTCal <- function(t_net,t_k,ID,beta_pho,beta_h){
        nest_parameter = integrate(tNetCal_vec,lower = t_net-beta_h, upper = t_net+beta_h, t_net = t_net, ID=ID,beta_pho = beta_pho,beta_h=beta_h,subdivisions=10000)
        a = tNetCal_vec(t_k,t_net,ID, beta_pho,beta_h)
        a*(nest_parameter$value)^(1/beta_pho-1)	
    }
    
    utilityTCal_vec = Vectorize(utilityTCal)
    
   #function to calculate numerator of PDF 
    numerator <- function(t_k,ID,beta_pho,beta_h){
        res = integrate(utilityTCal_vec,lower = t_k-beta_h, upper = t_k+beta_h, t_k=t_k, ID=ID,beta_pho = beta_pho,beta_h=beta_h,subdivisions=10000)
        res$value
    }
    
    #two functions together calculating denominator of PDF
    denom_utility <- function(q,ID, beta_pho,beta_h){
        utility = integrate(tNetCal_vec,lower = q-beta_h, upper = q+beta_h, t_net = q,ID=ID,beta_pho = beta_pho,beta_h=beta_h,subdivisions=10000)
        a = utility$value
        a^(1/beta_pho)
    }
    
    denom_utility_vec = Vectorize(denom_utility)
    
    denom <- function(L,U,ID,beta_pho,beta_h ){
        res = integrate(denom_utility_vec,lower = L, upper = U,ID = ID, beta_pho = beta_pho,beta_h=beta_h,subdivisions=10000)
        res$value
    }
    
    prob <- numerator(t,ID,beta_pho,beta_h)/denom(L,U,ID,beta_pho,beta_h)
}

CCNL.density_plot_vec = Vectorize(CCNL.density_plot)

# CCNL CDF function
CCNL.cdf_plot<-function(t,L,U,beta_h,beta_pho, a,b,c,d,e,ff,g,h,i,j,k,ID) {
    prob<-NULL
    for (i in 1:length(t)) {
        prob<-c(prob,integrate(CCNL.density_plot_vec,lower=L,upper=t[i],L=L,U=U,ID=ID,beta_h=beta_h,beta_pho=beta_pho,a=a,b=b,c=c,d=d,e=e,ff=ff,g=g,h=h,i=i,j=j,k=k,subdivisions=10000)$value)
    }
    return(prob)
}

# Test pdf function
test_PDF = CCNL.density_plot(t=0,L=0,U=24,beta_h = 0.75,beta_pho = 2, a=0,b=-0.5,c=0,d=3.70,e=2.43,ff=2.46,g=0.55,h=-4.73,i=-0.51,j=2.40,k=1.34,ID=1)

# Plot pdf function
curve(CCNL.density_plot_vec(t=x,L=0,U=24,beta_h = 0.75,beta_pho = 2, a=0,b=-0.5,c=0,d=3.70,e=2.43,ff=2.46,g=0.55,h=-4.73,i=-0.51,j=2.40,k=1.34,ID=1),from=0.00,to=24.00,n=100,add=FALSE,xlab="t",ylab="Density")

# Test cdf function
start.time <- Sys.time()
test_CDF = CCNL.cdf_plot(t=24,L=0,U=24,beta_h = 0.75,beta_pho = 2, a=0,b=-0.5,c=0,d=3.70,e=2.43,ff=2.46,g=0.55,h=-4.73,i=-0.51,j=2.40,k=1.34,ID=2)
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken
