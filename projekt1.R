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
# - Funkcja 'apply' biorąca payoff i estymatory, generuje n (10000) trajektorii, wylicza średni payoff

# Parametry mają format:
# = dryf
# = macierz kowariancji
# = ceny spot
# = średnie mtg
# = macierz kowariancji mtg
# = pochodna R-N
# = r

r = 0.01349
#=========Funkcje==================
#=========gentraj==================
gentraj=function(par,t,dt,method){
  # Metoda generowania wprost z miary mtg
  if(method==1){
    mu=par[[4]]
    si=par[[5]]
    pricestart=par[[3]]
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
  return(0)
}
m=c(0.1/200,0.15/200)
m=c(0,0)
s1=0.16/200
s2=0.2/200
si=matrix(c(s1,0.5*sqrt(s1*s2),0.5*sqrt(s1*s2),s2),ncol=2)
spot=c(1000,1500)
t=matrix(0,nrow=1000,ncol=2)
p=list(0,0,spot,m,si)
for(i in 1:1000){
  t[i,]=gentraj(p,100,1,1)[101,]
}
colMeans(t)
exp(s1*100)

test=gentraj(p,100,1,1)[101,]
#=========payapply=================
payapply=function(par,t,dt=1,method=1,payoff,assets=1,r=0,N=10000){
  day=1/251
  res=c()
  for(i in 1:N){
    traj=gentraj(par,t,dt,method)
    pay=payoff(traj[,assets])
    res[i]=pay
  }
  res=res*exp(-r*t*day)
  return(res)
}
#=========estim=======================
estim=function(pricedata,t,dt,r,t_spot=0){
  # data - każda kolumna zawiera wektor cen (zamknięcia) danego aktywa
  # t - w dniach okres z którego estymujemy
  # dt - interwał w dniach
  # t_spot - liczba (w dniach) od końca danych na kiedy estymujemy
  # wszystkie 
  N=ncol(pricedata)
  end=nrow(pricedata)
  amount=t*dt
  pricedata=pricedata[(end-t_spot*dt-amount):(end-t_spot*dt),]
  end=nrow(pricedata)
  # mu_c - dryfy historyczne
  mu_c=c()
  si_c=c()
  # rho - macierz kowariancji do generowania miary inwestora
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
  # Właściwe ustalenie rho, żeby mozna było generować:
  for(i in 1:N){
    rho[i,]=rho[i,]*si_c[i]
    rho[,i]=rho[,i]*si_c[i]
  }
  # Ceny spot:
  spotprice=pricedata[end,]
  # Jak obliczamy martyngalowe średnie i zmienności:
  mu_mtg=0
  si_mtg=0
  
  si_mtg_matrix=matrix(0,ncol=N,nrow=N)
  
}
#=========Moje testy==================
payoffcall=function(strike){
  f=function(x){
    last=x[length(x)]
    return(max(c(last-strike,0)))
  }
}
payofftest=payoffcall(2200)
# Sposób wzięcia długości albo liczby wierszy:
# min(c(nrow(spot),length(spot)))

setwd("./")
wig20 <- read.csv(file = 'wig20_koniec_2803.csv', sep=',')
#125:311
wig=wig20[124:311,5]
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

p_c=list(0,0,wig20[311,5],mu_c,matrix(si_c^2,nrow=1,ncol=1))
p_mtg=list(0,0,wig20[311,5],mu_mtg,matrix(si_mtg^2,nrow=1,ncol=1))
p_mtg2=list(0,0,wig20[311,5],mu_mtg2,matrix(si_mtg2^2,nrow=1,ncol=1))

p=p_mtg
pay=c()
for(i in 1:5000){
  t=gentraj(p,180,1,1)[101,1]
  pay[i]=payofftest(t)
}
mean(pay*exp(-r*0.75))
te=gentraj(p,180,1,1)
gentraj(p,180,1,1)


wig20new <- read.csv(file = '1001.csv', sep=',')

day=1/251

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

#=========Stare testy==================
g=function(y){
  f=function(x){x+y}
  return(f)
}
b=g(3)
b(4)

optionprice=function(strike,traj)

  
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


