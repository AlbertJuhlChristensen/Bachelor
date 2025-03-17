library(tidyverse)
library(tictoc)
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
        if (a_plus*x^2+2*x+(max(lam1+lam2,lamA+lamB)-1/dt)>0) {
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

#Funktion der spejler#
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

Strategi_funktion_x0 <- function(lam1,lam2,eps,z,w) {
  matrix <- x_værdier <- mod_matrix <- matrix(0, nrow=21,ncol =1/dt +1 )
  matrix[11,1/dt +1 ] <- mod_matrix[11,1/dt +1] <- 1+eps
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
#Ligevægtsfunktion 3 forskellige metoder. Alle metoder stopper når begge hold foretager 
#det samme to på hinanden følgende iterationer
#Metode 1: Hold A handler optimalt med model. Hold B handler optimalt men får ikke gavn af den effekt
#der kommer ved at hold A satser. De starter fra scratch #
dt <- 1/180
XP_ligevægt1 <- function(lam1,lam2,a_plus,a_minus,eps,dt){
  x_ny <- y_ny <- x <- y <- matrix(0,nrow=21,ncol=1/dt +1  )
  hold1 <- matrix(lam1,nrow=21,ncol=1/dt +1  )
  hold2 <- matrix(lam2,nrow=21,ncol=1/dt +1  )
  for (i in 1:50){
    XPA <- Strategi_funktion_matrix(lam1, lam2, hold1,hold2+spejl(y),a_plus,a_minus,eps,dt)$XP
    XPA_mod <- Strategi_funktion_matrix(lam1, lam2, hold1,hold2+spejl(y),a_plus,a_minus,eps,dt)$mod_XP
    x <- Strategi_funktion_matrix(lam1, lam2, hold1,hold2+spejl(y),a_plus,a_minus,eps,dt)$handling
    XPB <- Strategi_funktion_matrix(lam2, lam1, hold2,hold1+spejl(x),a_plus,a_minus,eps,dt)$XP
    XPB_mod <- Strategi_funktion_matrix(lam2, lam1, hold2,hold1+spejl(x),a_plus,a_minus,eps,dt)$mod_XP
    y <- Strategi_funktion_matrix(lam2,lam1, hold2,hold1+spejl(x),a_plus,a_minus,eps,dt)$handling
    print(paste("Iteration", i, ": XPA=", XPA[11,1], "XPB=", XPB[11,1], "XPA_mod=",XPA_mod[11,1], "XPB_mod=", XPB_mod[11,1]))
    if (all.equal(x_ny,x)==TRUE && all.equal(y_ny,y)==TRUE) {
      break
    }
    x_ny <- x
    y_ny <- y
  }
  return(list(XPA=XPA,XPA_mod=XPA_mod, XPB=XPB,XPB_mod=XPB_mod, hanA=x, hanB=y))
}
# Illustration af scenarie, hvor det ses, at de klarer sig dårligere i deres egen optimering grundet, 
# de ikke får den konvekse effekt #
test1 <- XP_ligevægt1(1.5,1.5,2,2,0,1/180)
Strategi_funktion_x0(1.5,1.5,0,test1$hanA, spejl(test1$hanB))$XPA[11,1]


# Metode 2: Hold A handler optimalt, og hold B handler derefter optimalt hvor de starter med udgangspunkt i
# at deres egen scoringsintensitet er blevet ændret pga hold A's valg. Men når hold A handler igen, så starter 
# de selv forfra og handler udelukkende ud fra hvordan hold B har påvirket scoringsintensiteterne. 
XP_ligevægt2 <- function(lam1,lam2,a_plus, a_minus ,eps, dt){
  x_ny <- y_ny <- xA <- xB <- yA <- yB <- matrix(0,nrow=21,ncol=1/dt +1 )
  holdA <- matrix(lam1, nrow=21, ncol=1/dt+1)
  holdB <- matrix(lam2, nrow=21, ncol=1/dt+1)
  for (i in 1:1000){
    A <- Strategi_funktion_matrix(lam1,lam2,holdA+spejl(xB),holdB + spejl(yB),a_plus,a_minus,eps,dt)
    xA <- A$handling
    yA <- xA + ifelse(xA>0,a_plus*(xA)^2,a_minus*(xA)^2)
    XPA <- A$XP[11,1]
    XPA_mod <- A$mod_XP[11,1]
    B <- Strategi_funktion_matrix(lam2,lam1,holdB+spejl(yA),holdA+spejl(xA),a_plus,a_minus,eps,dt)
    yB <- B$handling
    xB <- yB + ifelse(yB>0,a_plus*(yB)^2,a_minus*(yB)^2)
    XPB <- B$XP[11,1]
    XPB_mod <- B$mod_XP[11,1]
    print(paste( "Iteration",i, ":XPA=", XPA, "XPB=", XPB, "XPA_mod=",XPA_mod, "XPB_mod=", XPB_mod, "ssh=",max(A$ssh[,,1]+A$ssh[,,2]+A$ssh[,,3]))) 
    if (all.equal(x_ny,xA)==TRUE && all.equal(y_ny,yB)==TRUE) {
      break
    }
    x_ny <- xA
    y_ny <- yB
  }
  ssh <- array(NA, dim = c(21, 1/dt+1, 3))  
  ssh[,,1] <- (holdA+xA+spejl(xB))*dt #p_plus
  ssh[,,3] <- (holdB+yA + spejl(yB))*dt # p_minus
  ssh[,,2] <- max(0, 1-ssh[,,1]-ssh[,,3])
    return(list(XPA=A$XP,XPA_mod=A$mod_XP, XPB=B$XP,XPB_mod=B$mod_XP, xA=xA, yA=yA, xB=xB, yB=yB, ssh=ssh, it=i))
}
#Resultat af denne metode for samme situation #
test2 <- XP_ligevægt2(1.5,1.5,2,2,0,1/180)
test2$XPA[11,1]
Strategi_funktion_x0(1.5,1.5,0,test2$xA+spejl(test2$xB),test2$yA+spejl(test2$yB))$XPA[11,1]

#Metode 3: Akkumulering #:  
XP_ligevægt3 <- function(lam1,lam2,a_plus, a_minus ,eps, dt){
  x_ny <- y_ny <- xA <- xB <- yA <- yB <- matrix(0,nrow=21,ncol=1/dt +1 )
  holdA <- matrix(lam1, nrow=21, ncol=1/dt+1)
  holdB <- matrix(lam2, nrow=21, ncol=1/dt+1)
  for (i in 1:100){
    A <- Strategi_funktion_matrix(lam1, lam2, holdA,spejl(holdB),a_plus,a_minus,eps,dt)
    xA <- A$handling
    yA <- A$handling + ifelse(A$handling>0,a_plus*(A$handling)^2,a_minus*(A$handling)^2)
    XPA <- A$XP[11,1]
    print(i)
    XPA_mod <- A$mod_XP[11,1]
    B <- Strategi_funktion_matrix(lam2, lam1, holdB+spejl(yA),spejl(holdA)+spejl(xA),a_plus,a_minus,eps,dt)
    yB <- B$handling
    xB <- B$handling + ifelse(B$handling>0,a_plus*(B$handling)^2,a_minus*(B$handling)^2)
    XPB <- B$XP[11,1]
    XPB_mod <- B$mod_XP[11,1]
    print(paste( "XPA=", XPA, "XPB=", XPB, "XPA_mod=",XPA_mod, "XPB_mod=", XPB_mod, "ssh=",max(A$ssh[,,1]+A$ssh[,,2]+A$ssh[,,3]))) 
    if (all.equal(xA,x_ny)==TRUE & all.equal(yA,y_ny)==TRUE) {
      break
    }
    holdA <- holdA + xA + spejl(xB)
    holdB <- holdB + spejl(yA) + yB
    x_ny <- xA
    y_ny <- yA
  }
  print(i)
  return(list(XPA=A$XP,XPA_mod=A$mod_XP, XPB=B$XP,XPB_mod=B$mod_XP, xA=xA, yA=yA, xB=xB, yB=yB, it=i))
}
# Som beskrevet i projektet kan det ses nedenfor at denne metode aldrig vil stoppe, da trænerne
# altid vil ændre taktik grundet den konvekse effekt. 
test3 <- XP_ligevægt3(1.5,1.5,2,2,0,1/180)



