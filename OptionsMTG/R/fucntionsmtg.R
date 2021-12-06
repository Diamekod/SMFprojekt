#' Generating price trajectory
#'
#' Generates random price trajectory for a given time period
#'
#' @param par parameters in format as result of the estim function
#' @param t the time (in days) of generated trajectory
#' @param dt the time (in days) between checking prices
#' @param method 1 - computes price using martingale measure directly, 2 - computes prices from original measure and Radon-Nikodym derivative
#' @return For method 1 returns matrix of prices (columns~assets), for method 2 returns list of price matrix and R-N derivative
#' @examples
#' rn <- function(x) 1
#' p <- list(-0.000228,0.00956,matrix(1,1,1),2204,0,matrix(0.000108,1,1),rn,0.02)
#' tr <- gentraj(p,200,1)
#' @export
gentraj=function(par,t,dt,method=1){
  # Dla metody 1 zwraca macierz cen
  # Dla metody 2 zwraca listę - macierz cen wg miary inwestora
  # oraz pochodną RN wyliczoną na cały okres - iloczyn pochodnych wyliczonych dla zwrotów w poszczególnych chwilach
  # Metoda generowania wprost z miary mtg
  if(method==1){
    mu=par[[5]]
    si=par[[6]]
    pricestart=par[[4]]
    # s_{i+1}=s_i*exp(N)
    # generujemy N z odpowiedniego rozkładu, żeby było mtg
    amount=t/dt
    # Kolumna ma długość równą liczbie aktywów.
    res=matrix(0,ncol=length(mu),nrow=amount+1)
    res[1,]=pricestart
    rand=mvtnorm::rmvnorm(amount,mu,si)
    rand=apply(rand,2,cumsum)
    # print(colMeans(rand))
    rand=exp(rand)
    res[-1,]=t(t(rand)*pricestart)
    return(res)
  }
  if(method==2){
    mu=par[[1]]
    si=par[[2]]
    rho=par[[3]]
    pricestart=par[[4]]
    dRN=par[[7]]
    # s_{i+1}=s_i*exp(N)
    # generujemy N z odpowiedniego rozkładu, żeby było mtg
    amount=t/dt
    # Kolumna ma długość równą liczbie aktywów.
    res=matrix(0,ncol=length(mu),nrow=amount+1)
    res_dRN=1
    res[1,]=pricestart
    rand=mvtnorm::rmvnorm(amount,rep(0,length(mu)),rho)
    res_dRN=apply(rand,1,dRN)
    res_dRN=prod(res_dRN)
    rand=t(t(rand)*si+mu)
    rand=apply(rand,2,cumsum)
    # print(colMeans(rand))
    rand=exp(rand)
    res[-1,]=t(t(rand)*pricestart)
    return(list(res,res_dRN))
  }
  return(NULL)
}

#' Estimating option price
#'
#' Computes average option payoff
#'
#' @param par parameters in format as result of the estim function
#' @param t the time (in days) of generated trajectory
#' @param dt the time (in days) between checking prices
#' @param method 1 - computes price using martingale measure directly, 2 - computes prices from original measure and Radon-Nikodym derivative
#' @param payoff function: trajectory |-> option payoff
#' @param assets vector of assets to which trajectories should be applied payoff function
#' @param N number of repeats
#' @return Returns vector of option prices at t=0
#' @examples
#' rn <- function(x) 1
#' p <- list(-0.000228,0.00956,matrix(1,1,1),2204,0,matrix(0.000108,1,1),rn,0.02)
#' payoff2200 <- function(x) max(0,x[length(x)]-2200)
#' prices <- payapply(p,200,payoff=payoff2200,N=20000)
#' @export
payapply=function(par,t,dt=1,method=1,payoff,assets=1,N=10000){
  day=1/251
  res=c()
  r=par[[8]]
  if(method==1){
    for(i in 1:N){
      traj=gentraj(par,t,dt,method)
      pay=payoff(traj[,assets])
      res[i]=pay
    }
    res=res*exp(-r*t*day)
    return(res)
  }
  if(method==2){
    for(i in 1:N){
      traj=gentraj(par,t,dt,method)
      dRN=traj[[2]]
      tr=traj[[1]]
      pay=payoff(tr[,assets])
      res[i]=pay*dRN
      #print(i)
    }
    res=res*exp(-r*t*day)
    return(res)
  }
  return(NULL)
}

#' Generating parameters
#'
#' Generates parameters of martingale measure from historical parameters
#' Mostly used in estim function.
#'
#' @param mu average historical return
#' @param si sd of historical returns
#' @param rho covariance matrix (diagonal=1) of historical returns
#' @param spot vector of spot prices
#' @param r risk-free interest rate
#' @param dt the time (in days) between checking prices
#' @return
#' Returns list containing: mean of historical returns, sd of historical returns, covariance matrix (1 on diagonal),
#' spot prices, mean of returns for mtg measure, covariance matrix of returns in mtg measure, R-N derivative, risk-free interest rate
#' @export
calcmtg=function(mu,si,rho,spot,r,dt){
  day=1/250
  nassets=length(mu)
  # Parametry zwrotu w mierze mtg (po uproszczeniu):
  mu_mtg=rep(0,nassets)
  si_mtg=rho
  si_mtg=si_mtg*2*r*dt*day
  # Pochodna R-N wykorzystuje gęstość rozkładów normalnych:
  # standardowego i pochodzącego z miary mtg. Ale nie jest to rozkład z parametrami wyżej, tylko
  # Parametry miary Q do pochodnej RN:
  mu_qrn=-mu/si
  si_q=sqrt(2*r*dt*day)/si
  si_qrn=rho
  for(i in 1:nassets){
    si_qrn[i,]=si_qrn[i,]*si_q
    si_qrn[,i]=si_qrn[,i]*si_q
  }
  dqdp=function(x){
    mvtnorm::dmvnorm(x,mu_qrn,si_qrn)/mvtnorm::dmvnorm(x,sigma=rho)
  }
  res=list(mu,si,rho,spot,mu_mtg,si_mtg,dqdp,r)
  return(res)
}

#' Generating parameters
#'
#' Generates parameters of martingale measure from historical data
#'
#' @param pricedata matrix of historical prices (columns~assets)
#' @param r risk-free interest rate
#' @param t the time (in days) from which the parameters are computed
#' @param dt the time (in days) between checking prices
#' @return
#' Returns list containing: mean of historical returns, sd of historical returns, covariance matrix (1 on diagonal),
#' spot prices, mean of returns for mtg measure, covariance matrix of returns in mtg measure, R-N derivative, risk-free interest rate
#' @examples
#' data(wig)
#' par <- estim(wig_m,0.02)
#' par2<- estim(wigkgh_m,0.0135)
#' @export
estim=function(pricedata,r,t=0,dt=1){
  # data - każda kolumna zawiera wektor cen (zamknięcia) danego aktywa
  # t - w dniach okres z którego estymujemy
  # dt - interwał w dniach
  day=1/250
  N=ncol(pricedata)
  #print(N)
  end=nrow(pricedata)
  #print(end)
  amount=t*dt
  if(amount==0){
    amount = end-1
  }
  pricedata_tmp=as.matrix(pricedata[(end-amount):end,])
  end=nrow(pricedata_tmp)
  #print(pricedata_tmp)
  # mu_c - dryfy historyczne
  mu_c=c()
  si_c=c()
  # rho na razie zawiera tylko korelacje między aktywami
  # Przemnożenie przez zmienności, żeby mieć macierz do generowania nastąpi później
  rho=matrix(0,nrow=N,ncol=N)
  ret=matrix(0,ncol=N,nrow=nrow(pricedata)-1)
  for(i in 1:N){
    ret[,i]=log(pricedata[-1,i]/pricedata[-end,i])
    mu_c[i]=mean(ret[,i])
    si_c[i]=sd(ret[,i])
  }
  for(i in 1:N){
    for (j in i:N){
      rho[i,j]=rho[j,i]=cor(ret[,i],y=ret[,j])
    }
  }
  diag(rho)=1
  # Ceny spot:
  spotprice=pricedata[end,]
  res=calcmtg(mu_c,si_c,rho,spotprice,r,dt)
  return(res)
}
