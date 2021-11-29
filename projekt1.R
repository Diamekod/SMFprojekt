library(mvtnorm)
#=========Opis==================
# Wejście:
# - liczba używanych aktywów
# - csv z danymi (jeden na aktywo)
# - interwal czasu zbierania danych (w dniach gieldowych) (lub lista, jeśli różne)
# - czas z ktorego maja byc wyliczane parametry (od konca danych wstecz)
# - Funkcja payoffu od trajektorii :
#   - domyślnie od aktywa 1, lista aktywów będących argumentami payoffu jako parametr
# - Kilka domyślnych funkcji do generowania payoffu opcji call/put, z barierami

# Parametry opcji:
# - call/put
# - czas rozpoczcia opcji (w dniach od zakończenia danych)
# - czas zapadniecia opcji ---||---
# - strike
# - bariera down/up/in/out/-
# - stopa proc.

# Funkcje:
# - generowanie parametrów na podstawie danych: zwraca listę średnich, macierz cov, ceny spot oraz parametrów wyznaczających miarę martyngałową
#   przyjmując długość odcinka do estymacji (w dniach), interwał zbierania danych (domyślnie 1 dzień) i moment na który liczymy (domyślnie 0 - koniec danych), oraz r
# - generowanie jednej trajektorii dla zadanych parametrów - dwie metody.
# Uwaga: przy generowaniu trajektorii mamy losować z rozkładu normalnego (standardowego, potem innego), a następnie je modyfikować
# Matematycznie to jest to samo, co losowanie ze zmodyfikowanego rozkładu normalnego. Tak zrobimy.
# - Funkcja 'apply' biorąca payoff i estymatory, generuje n (10000) trajektorii, wylicza średni payoff zdyskontowany

# Parametry mają format:
# = dryf
# = historyczna zmienność
# = macierz kowariancji (miary inwestora)
# = ceny spot
# = średnie mtg (== 0)
# = macierz kowariancji mtg
# = pochodna R-N
# = r

r = 0.01349
#=========Funkcje==================
#=========gentraj==================
gentraj=function(par,t,dt,method){
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
    rand=rmvnorm(amount,mu,si)
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
    rand=rmvnorm(amount,rep(0,length(mu)),rho)
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
#=========payapply=================
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
#=========estim=======================
calcmtg=function(mu,si,rho,spot,r,dt,nassets){
  day=1/250
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
    dmvnorm(x,mu_qrn,si_qrn)/dmvnorm(x)
  }
  res=list(mu,si,rho,spot,mu_mtg,si_mtg,dqdp,r)
  return(res)
}
estim=function(pricedata,r,t=0,dt=1,t_spot=0){
  # data - każda kolumna zawiera wektor cen (zamknięcia) danego aktywa
  # t - w dniach okres z którego estymujemy
  # dt - interwał w dniach
  # t_spot - liczba (w dniach) od końca danych na kiedy estymujemy
  # wszystkie
  day=1/250
  N=ncol(pricedata)
  #print(N)
  end=nrow(pricedata)
  #print(end)
  amount=t*dt
  if(amount==0){
    amount = end-1
  }
  pricedata_tmp=as.matrix(pricedata[(end-t_spot*dt-amount):(end-t_spot*dt),])
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
    si_c=sd(ret[,i])
  }
  for(i in 1:N){
    for (j in i:N){
      rho[i,j]=rho[j,i]=cor(ret[,i],y=ret[,j])
    }
  }
  diag(rho)=1
  # Ceny spot:
  spotprice=pricedata[end,]
  res=calcmtg(mu_c,si_c,rho,spotprice,r,dt,N)
  return(res)
}
#=========Moje testy==================
payoffcall=function(strike){
  f=function(x){
    last=x[length(x)]
    return(max(c(last-strike,0)))
  }
}
calluao=function(strike,bar){
  f=function(x){
    last=x[length(x)]
    if(max(x)>=strike+bar){
      return(0)
    }
    return(max(c(last-strike,0)))
  }
}
calluai=function(strike,bar){
  f=function(x){
    last=x[length(x)]
    if(max(x)>=strike+bar){
      return(max(c(last-strike,0)))
    }
    return(0)
  }
}
callpar=function(strike,bar,t){
  f=function(x){
    last=x[length(x)]
    t_up=which(x>strike+bar)
    cnt=1
    n=length(t_up)-1
    if(n<t) return(0)
    #print(n)
    for(i in 1:n){
      if(t_up[i+1]==(t_up[i]+1)){
        cnt=cnt+1
      } else{
        cnt=1
      }
      if(cnt==t) return(max(c(last-strike,0)))
    }
    return(0)
  }
}
calllookback=function(strike){
  f=function(x){
    t=nrow(x)
    if(x[t,2]>x[1,2]){
      return(max(c(x[t,1]-min(x[,1]),0)))
    }
    return(0)
  }
}

payofftest=payoffcall(2200)
payoffpar=callpar(2200,200,10)
payoffuao=calluao(2200,200)
payoffuai=calluai(2200,200)
payofflb=calllookback(2200)

# Sposób wzięcia długości albo liczby wierszy:
# min(c(nrow(spot),length(spot)))
setwd("./")
wig20 <- read.csv(file = 'wig20_koniec_2803.csv', sep=',')
kghm <- read.csv(file = 'kgh_koniec_2803.csv', sep=',')
#125:311
r=0.0135
texp=round(3/4*250)
wig=wig20[124:311,5]
kgh=kghm[124:311,5]
wig_m=matrix(wig,ncol=1)
param=estim(pricedata=wig_m,r)

t1=gentraj(param,texp,1,1)
t2=gentraj(param,texp,1,2)

# EC up and out:
call2200=payapply(param,texp,payoff=payoffuao,method=1)
mean(call2200)

# Różnica w tempie zbieżności:
call2200_1=payapply(param,texp,payoff=payofftest,method=1,N=100)
call2200_2=payapply(param,texp,payoff=payofftest,method=2,N=100)

mean(call2200_1)
mean(call2200_2)

library(fOptions)
sbs=param[[2]]*sqrt(250)
GBSOption('c',wig[188],2200,texp/250,r,r,sbs)@price

wigkgh_m=matrix(c(wig,kgh),ncol=2)
param2=estim(pricedata=wigkgh_m,r)
t2=gentraj(param2,texp,1,1)
# EC w modelu dwuaktywowym:
mean(payapply(param2,texp,payoff=payofftest))
# Opcja lookback:
mean(payapply(param2,texp,payoff=payofflb,assets=c(1,2)))




#=========Stare testy==================
g=function(y){
  f=function(x){x+y}
  return(f)
}
b=g(3)
b(4)

si=matrix(c(1,0.5,0.5,1),ncol=2)
m=c(1,0.3)
r=rmvnorm(10000,m,si)
rmvnorm(3,0)
sum(abs(colMeans(r)-m))

a=matrix(1:6,ncol=2)
a
cor(a[,1],y=a[,2])
a*c(1,2)
# Przemnożenie przez stałe kolumnowo (tzn. kol. 1 przez stałą 1, itd.)
t(t(a)*c(1,2))
# Sumowanie kolumnowo:
apply(a,2,cumsum)

#przyrosty:
pr=wig[-1]/wig[-length(wig)]
pr=log(pr)
mu_c=mean(pr)
si_c=sd(pr)
mu_c*250
si_c*sqrt(250)
# parametry mtg:
r=0.0153
r_d=r/250
mu_mtg=0
si_mtg=sqrt(2*r_d)
si_mtg
# Inna wersja - zadajemy zmienność. Wtedy dryf robi się strasznie duży:
si_mtg2=0.01
mu_mtg2=(r_d-mu_c-(si_c^2*si_mtg2^2)/2)/si_c

p_c=list(0,0,0,wig20[311,5],mu_c,matrix(si_c^2,nrow=1,ncol=1))
p_mtg=list(0,0,0,wig20[311,5],mu_mtg,matrix(si_mtg^2,nrow=1,ncol=1))
p_mtg2=list(0,0,0,wig20[311,5],mu_mtg2,matrix(si_mtg2^2,nrow=1,ncol=1))

p=p_mtg
pay=c()
for(i in 1:5000){
  t=gentraj(p,180,1,1)[101,1]
  pay[i]=payofftest(t)
}
mean(pay*exp(-r*0.75))
te=gentraj(p,180,1,1)
gentraj(p,180,1,1)

  
cena <- wig20$Otwarcie
N = length(cena)
zwrot <- log(cena[-1]/cena[-N])

zwrot_rok <- zwrot[251:500]
dryf = mean(zwrot_rok)
zmiennosc = sd(zwrot_rok)

dt = 1/251
dryf_r = dryf/dt
zmiennosc_r = zmiennosc/sqrt(dt)



cenadzis <- 2238.85
t0 <- 64/251
T <- 1

# Funkcje do B-S:
d1=function(S,E,r,sig,Tt,t){
  (log(S/E)+(r+0.5*sig^2)*(Tt-t))/(sig*sqrt(Tt-t))
}
d2=function(S,E,r,sig,Tt,t){
  d1(S,E,r,sig,Tt,t)-sig*sqrt(Tt-t)
}
delta<-function(S,E,r,sig,Tt,t){
  pnorm(d1(S,E,r,sig,Tt,t))
}
cena_opcji_call<-function(S,E,r,sig,Tt,t=t0){
  S*pnorm(d1(S,E,r,sig,Tt,t))-E*exp(-r*(Tt-t))*pnorm(d2(S,E,r,sig,Tt,t))
}
delta_put<-function(S,E,r,sig,Tt,t){
  pnorm(d1(S,E,r,sig,Tt,t))-1
}
cena_opcji_put<-function(S,E,r,sig,Tt,t){
  S*(pnorm(d1(S,E,r,sig,Tt,t))-1)-E*exp(-r*(Tt-t))*(pnorm(d2(S,E,r,sig,Tt,t))-1)
}

cena <- wig20new$Otwarcie
N = length(cena)
zwrot <- log(cena[-1] / cena[-N])

zwrot_rok <- zwrot[1:187]
dryf = mean(zwrot_rok)
zmiennosc = sd(zwrot_rok)

dt = 1/251
dryf_r = dryf/dt
zmiennosc_r = zmiennosc/sqrt(dt)

r = 0.0153
cenadzis <- 2290.80
t0 <- 188/251
T <- 1
ro <- T-t0


