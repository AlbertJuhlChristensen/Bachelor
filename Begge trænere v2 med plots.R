library(tidyverse)
library(tictoc)
library(ggplot2)
# Oprettelse af den optimale strategifunktion, men hvor input også er en matrix, #
# for både hold A og hold B #
Strategi_funktion_matrix <- function(lamA,lamB, holdA,holdB,a_plus,a_minus,eps,dt) {
  matrix <- x_værdier <- mod_matrix <- matrix(0,nrow=21,ncol=1/dt + 1)
  matrix[11,1/dt +1] <- mod_matrix[11,1/dt +1] <- 1+eps
  matrix[12:21,1/dt+1] <-3
  ssh <- array(0, dim = c(21, 180, 3))
  matrix[21,] <- 3
  mod_matrix[1,] <- 3
  mod_matrix[1:10,1/dt + 1] <- 3
  x <- 0
  for (q in (1/dt):1) {
    for (p in 2:20) {
      lam1 <- holdA[p,q]
      lam2 <- holdB[p,q]
      if (2*matrix[p,q+1]-matrix[p-1,q+1]-matrix[p+1,q+1]<0) {
        x <- (2*matrix[p,q+1]-matrix[p-1,q+1]-matrix[p+1,q+1]) / (-abs(2*a_plus*(matrix[p-1,q+1]-matrix[p,q+1]-10^(-10))))
        if (a_plus*x^2+2*x+max(lam1+lam2,lamA+lamB)-1/dt>0) {
          d <- 2^2-4*a_plus*(max(lam1+lam2,lamA+lamB)-1/dt)
          højrerod <- (-2+sqrt(d))/(2*a_plus)
          x <- højrerod        }
      }
      else if (2*matrix[p,q+1]-matrix[p-1,q+1]-matrix[p+1,q+1]>0){
        x <- (2*matrix[p,q+1]-matrix[p-1,q+1]-matrix[p+1,q+1]) / (-abs(2*a_minus*(matrix[p-1,q+1]-matrix[p,q+1]-10^(-10))))
        x <- max(x,-lamA)
        if (lamB+x+a_minus*x^2 < 0) {
          d <- 1-4*a_minus*lamB
          højrerod <- (-1+sqrt(d))/(2*a_minus)
          x <- højrerod
        }
      }
      else {
        x <- 0
      }
      x_værdier[p,q] <- x
      lam1x <- lam1 + x
      lam2x <- lam2 + x
      if (x>0){
        lam2x <- lam2x + a_plus*x^2
      }
      if (x<0){
        lam2x <- lam2x + a_minus*x^2
      }
      p_plus <- (lam1x)*dt
      p_minus <- (lam2x)*dt
      p_0 <- max(0,1-p_plus-p_minus)
      ssh[p,q,1] <- p_plus
      ssh[p,q,2] <- p_0
      ssh[p,q,3] <- p_minus
      matrix[p,q] <- p_plus*matrix[p+1,q+1]+p_0*matrix[p,q+1]+p_minus*matrix[p-1,q+1]
      mod_matrix[p,q] <- p_plus*mod_matrix[p+1,q+1]+ p_0*mod_matrix[p,q+1]+p_minus*mod_matrix[p-1,q+1]
    }
  }
  return(list(XP=matrix,handling=x_værdier, mod_XP=mod_matrix, ssh=ssh))
}
# Funktion der kan spejle en matrix #
spejl <- function(x){
  spejl <- x
  spejl[1:nrow(x),] <- x[nrow(x):1,]
  return(spejl) }

#Oprettelse af XP for "klassisk" Poisson model #
ssh_pois <- function(lam1,lam2) {
  p_1 <- p_X <- p_2 <- 0
  for (i in 0:20){
    for (j in 0:20) {
      if (j < i)  {
        p_1 <-p_1+dpois(i,lam1)*dpois(j,lam2) } #ssh for hjemmemål>udemål
      if (j>i) {
        p_2  <- p_2 + dpois(i,lam1)*dpois(j,lam2)}
    }
  }
  p_X <- 1-p_1-p_2
  return(c(p_1,p_X,p_2))
}

XP_pois <- function(lam1,lam2){
  XP <- c(0,0)
  XP[1] <- 3*ssh_pois(lam1,lam2)[1]+ssh_pois(lam1,lam2)[2]
  XP[2] <- 3*ssh_pois(lam1,lam2)[3] + ssh_pois(lam1,lam2)[2]
  return(XP)
}
# Funktion der tager begge holds optimale handlinger og kan udregne forventet værdi # 
Strategi_funktion_x0 <- function(lam1,lam2,z,w,dt) {
  matrix <- x_værdier <- mod_matrix <- matrix(0, nrow=21,ncol =1/dt +1 )
  matrix[11,1/dt +1 ] <- mod_matrix[11,1/dt +1] <- 1
  matrix[12:21,1/dt +1] <- 3
  matrix[21,] <- 3
  mod_matrix[1,] <- 3
  mod_matrix[1:10,1/dt +1 ] <- 3
  for (q in (1/dt):1) {
    for (p in 2:20) {
      lam1_ny <- lam1 + z[p,q]
      lam2_ny <- lam2 + w[p,q]
      p_plus <- lam1_ny*dt
      p_minus <- lam2_ny*dt
      p_0 <- max(0,1-p_plus-p_minus)
      matrix[p,q] <- p_plus*matrix[p+1,q+1]+ p_0*matrix[p,q+1]+p_minus*matrix[p-1,q+1] 
      mod_matrix[p,q] <- p_plus*mod_matrix[p+1,q+1]+p_0*mod_matrix[p,q+1]+p_minus*mod_matrix[p-1,q+1]
    }
  }
  return(list(XPA=matrix, XPB=mod_matrix))
}
#Ligevægtsfunktion med to forskellige eps som input
XP_ligevægt <- function(lam1,lam2,a_plus, a_minus ,eps1,eps2, dt){
  x_ny <- y_ny <- xA <- xB <- yA <- yB <- matrix(0,nrow=21,ncol=1/dt +1 )
  holdA <- matrix(lam1, nrow=21, ncol=1/dt+1)
  holdB <- matrix(lam2, nrow=21, ncol=1/dt+1)
  for (i in 1:1000){
    A <- Strategi_funktion_matrix(lam1,lam2,holdA+spejl(xB),holdB + spejl(yB),a_plus,a_minus,eps1,dt)
    xA <- A$handling
    yA <- xA + ifelse(xA>0,a_plus*(xA)^2,a_minus*(xA)^2)
    XPA <- A$XP[11,1]
    XPA_mod <- A$mod_XP[11,1]
    B <- Strategi_funktion_matrix(lam2,lam1,holdB+spejl(yA),holdA+spejl(xA),a_plus,a_minus,eps2,dt)
    yB <- B$handling
    xB <- yB + ifelse(yB>0,a_plus*(yB)^2,a_minus*(yB)^2)
    XPB <- B$XP[11,1]
    XPB_mod <- B$mod_XP[11,1]
    #print(paste("XPA=",XPA,"XPB=",XPB))
    if (all.equal(x_ny,xA, tolerance=10^-6)==TRUE && all.equal(y_ny,yB, tolerance=10^-6)==TRUE) {
      break
    }
    x_ny <- xA
    y_ny <- yB
  }
  #print(i)
  ssh <- array(NA, dim = c(21, 1/dt+1, 3))  
  ssh[,,1] <- (holdA+xA+spejl(xB))*dt #p_plus
  ssh[,,3] <- (holdB+yA + spejl(yB))*dt # p_minus
  ssh[,,2] <- max(0, 1-ssh[,,1]-ssh[,,3])
  return(list(XPA=A$XP,XPA_mod=A$mod_XP, XPB=B$XP,XPB_mod=B$mod_XP, xA=xA, yA=yA, xB=xB, yB=yB, ssh=ssh, it=i))
}


# Oprettelse af en masse kombinationer for at teste om ligevægt opnås #
lam1_grid <- seq(0.5, 2.5, by = 1)
lam2_grid <- seq(0.5, 2.5, by = 1)
a_plus_grid <- seq(0.5, 4.5, by = 1)
a_minus_grid <- seq(0.5, 4.5, by = 1)
eps1_grid <- seq(-1,2, by=1)
eps2_grid <- seq(-1,2, by=1)

kombinationer <- expand.grid(lam1 = lam1_grid, lam2 = lam2_grid, 
                             a_plus = a_plus_grid, a_minus=a_minus_grid, eps1=eps1_grid, eps2=eps2_grid)

# Histogram over antallet af iterationer før ligevægt opnås #
kombinationer$it <- mapply(function(lam1, lam2, a_plus, a_minus, eps1,eps2) {
  værdi <- XP_ligevægt(lam1, lam2, a_plus, a_minus, eps1,eps2, 1/180)
  return(værdi$it)
            }, 
          lam1 = kombinationer$lam1, 
          lam2 = kombinationer$lam2, 
          a_plus = kombinationer$a_plus, 
          a_minus = kombinationer$a_minus,
          eps1 = kombinationer$eps1,
          eps2 = kombinationer$eps2 )

ggplot(data=count(kombinationer,it), aes(x=as.factor(it), y=n)) + geom_bar(stat="identity", fill="blue")+
  labs(x="Iterationer", y="Antal forekomster")+theme_bw()

# Plot for ligevægt sammenlignet med den simple model med eps=0 svarende til figur 10 #
data <- data.frame(XP=numeric(), lam_forskel=numeric(), metode=numeric(), a_forskel=numeric())
lam_sekvens <- seq(-2,2, by=1)
a_sekvens <- seq(0,4, by=0.5)
for (a in a_sekvens){
  a_plus <- 2.5 + a/2 
  a_minus <- 2.5 - a/2
  for (lam in lam_sekvens){
    lam1 <- 1.5 + lam/2
    lam2 <- 1.5 - lam/2
    ligevægt <- XP_ligevægt(lam1,lam2,a_plus,a_minus,0,0,1/180)$XPA[11,1]
    optimal <- Strategi_funktion_matrix(lam1,lam2,matrix(lam1,21,181),matrix(lam2,21,181),a_plus,a_minus,0,1/180)$XP[11,1]
    data <- rbind(data, 
                  data.frame(XP=ligevægt-XP_pois(lam1,lam2)[1],lam_forskel=lam1-lam2, metode="Begge trænere", a_forskel=a_plus-a_minus),
                  data.frame(XP=optimal-XP_pois(lam1,lam2)[1],lam_forskel=lam1-lam2, metode="Èn træner", a_forskel=a_plus-a_minus))
  }
}
ggplot(data=data, aes(x=a_forskel, y=XP, color=factor(lam_forskel), linetype=metode, group=interaction(metode,factor(lam_forskel)))) + 
  geom_line(size=0.7) + labs(x="Forskel i a", y="Gevinst i forventet point", color="Forskel i lambda", linetype="Optimering") +
  theme_bw()+theme(text=element_text(size=14))

# Plot for forskellige epsilonværdier - svarende til figur 11 #
data2 <- data.frame(XP=numeric(), lam_forskel=numeric(),eps1=numeric(),eps2=numeric(), lam_forskel=numeric())
eps1_grid <- seq(-1,1,by=0.1)
eps2_grid <- seq(-1,1, by=1)
lam_sekvens2 <- seq(-1,1, by=1)
for (eps1 in eps1_grid){
  for (eps2 in eps2_grid){
    for (lam in lam_sekvens2){
      lam1 <- 1.5 + lam/2
      lam2 <- 1.5 - lam/2
      kamp <- XP_ligevægt(lam1,lam2,2.5,2.5,eps1,eps2,1/180)
      z <- kamp$xA + spejl(kamp$xB)
      w <- kamp$yA + spejl(kamp$yB)
      XP <- Strategi_funktion_x0(lam1,lam2,z,w,1/180)
      data2 <- rbind(data2, 
                     data.frame(XP=XP$XPA[11,1]-XP_pois(lam1,lam2)[1], 
                                eps1=eps1,eps2=eps2, lam_forskel=lam1-lam2))
  }
  }
}

ggplot(data=data2, aes(x=eps1, y=XP, color=factor(lam_forskel), linetype=factor(eps2),group=interaction(factor(lam_forskel),factor(eps2))))+
  geom_line(size=0.7)+theme_bw()+scale_linetype_manual(values = c("solid", "twodash", "dotted"))+
  labs(x="Epsilon hold A", y="Gevinst i forventet point", color="Forskel i lambda", linetype="Epsilon hold B")+
  theme_bw()+theme(text=element_text(size=14))





