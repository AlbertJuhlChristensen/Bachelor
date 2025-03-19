library(tictoc)
library(ggplot2)
#Funktion til at kunne tage rødder i andengradsligning#
rødderfkt <- function(a,b,c){
  d <- b^2-4*a*c
  if (d<0){
    return(c(NA,NA))
  }
  else {
    return((-b+c(-1,1)*sqrt(d))/(2*a))
  }
}
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
# Kode svarende til metode 2 - numerisk løsning #
Strategi_funktion1 <- function(lam1,lam2,a_plus,a_minus,eps,dt) {
  matrix <- x_værdier <- mod_matrix <- matrix(0,nrow=21,ncol =1/dt+1)
  matrix[11,1/dt+1] <- mod_matrix[11,1/dt+1] <- 1+eps
  matrix[12:21,1/dt+1] <-3
  matrix[21,] <- 3
  mod_matrix[1,] <- 3
  mod_matrix[1:10,1/dt+1] <- 3
  lambda_x <- c(0,0)
  nedre_grænse <- max(dt-lam1,rødderfkt(a_minus,1,lam2-dt)[2],na.rm = TRUE)
  øvre_grænse <- 1/dt - lam1
    for (q in (1/dt):1) {
      for (p in 2:20) {
        max_funktion <- function(x) {
          lambda_x[1] <- lam1 + x
          if (x>0) {
            lambda_x[2] <- lam2 + x + a_plus*x^2}
          else {
            lambda_x[2] <- lam2 + x + a_minus*x^2}
          p_0 <- exp(-(lambda_x[1]+lambda_x[2])*dt)
          p_plus <- lambda_x[1]/(lambda_x[1]+lambda_x[2]) * (1-p_0)
          p_minus <- 1-p_0-p_plus
          ssh <- c(p_plus,p_0,p_minus)
          max <- sum(ssh*c(matrix[p+1,q+1], matrix[p,q+1], matrix[p-1,q+1]))
          return(list(XP=max,ssh=ssh))
        }
        løsning <- optimise(function(x) max_funktion(x)$XP
                            ,c(nedre_grænse,øvre_grænse), 
                            maximum = TRUE)
        matrix[p,q] <- løsning$objective
        x <- x_værdier[p,q] <- løsning$maximum
        lambda_x[1] <- lam1 + x
        if (x>0) {
          lambda_x[2] <- lam2 + x + a_plus*x^2}
        else {
          lambda_x[2] <- lam2 + x + a_minus*x^2}
        p_0 <- exp(-(lambda_x[1]+lambda_x[2])*dt)
        p_plus <- lambda_x[1]/(lambda_x[1]+lambda_x[2]) * (1-p_0)
        p_minus <- max(0,1-p_0-p_plus)
        ssh <- c(p_plus,p_0,p_minus)
        mod_matrix[p,q] <- sum(ssh*c(mod_matrix[p+1,q+1], mod_matrix[p,q+1], mod_matrix[p-1,q+1]))
      }
    }
  return(list(XP=matrix,handling=x_værdier, mod_XP=mod_matrix, ssh=ssh))
}

# Kode uden optimal handling for at kunne tjekke efter #
Strategi_funktion_x0 <- function(lam1,lam2,eps) {
  matrix <- x_værdier <- mod_matrix <- matrix(0,nrow=21,ncol =91)
  matrix[11,91] <- mod_matrix[11,91] <- 1+eps
  matrix[12:21,91] <-3
  matrix[21,] <- 3
  mod_matrix[1,] <- 3
  mod_matrix[1:10,91] <- 3
  lambda_x <- c(lam1,lam2)
  for (q in 90:1) {
    for (p in 2:20) {
      p_0 <- exp(-(lam1+lam2)*dt)
      p_plus <- lam1/(lam1+lam2) * (1-p_0)
      p_minus <- max(0,1-p_0-p_plus)
      ssh <- c(p_plus,p_0,p_minus)
      matrix[p,q] <- sum(ssh*c(matrix[p+1,q+1], matrix[p,q+1], matrix[p-1,q+1])) 
      mod_matrix[p,q] <- sum(ssh*c(mod_matrix[p+1,q+1], mod_matrix[p,q+1], mod_matrix[p-1,q+1]))
    }
  }
    return(list(XP=matrix, mod_XP=mod_matrix))
  }

# Kode svarende til metode 1 - manuelt forsøg på optimal x-værdi #
Strategi_funktion2 <- function(lam1,lam2,a_plus,a_minus,eps,dt) {
  matrix <- x_værdier <- mod_matrix <- matrix(0,nrow=21,ncol =1/dt+1)
  matrix[11,1/dt+1] <- mod_matrix[11,1/dt+1] <- 1+eps
  matrix[12:21,1/dt+1] <-3
  matrix[21,] <- 3
  mod_matrix[1,] <- 3
  mod_matrix[1:10,1/dt+1] <- 3
  x_vektor <- test <- (-50:100)/10
  for (q in (1/dt):1) {
    for (p in 2:20) {
      for (n in 1:length(x_vektor)) {
      x <- x_vektor[n]
      lam1x <- lam1 + x
      lam2x <- lam2 + x
      if (x>0) {
        lam2x <- lam2x+a_plus*x^2
      }
      if (x<0){
        lam2x <- lam2x + a_minus*x^2
      }
      p_plus <- (lam1x)*dt; straf <- (p_plus>1)+(p_plus<0)
      p_minus <- (lam2x)*dt ; straf <- straf +(p_minus>1)+(p_minus<0)
      p_0 <- 1-p_plus-p_minus ; straf <- straf + (p_0>1) + (p_0<0)
      test[n] <- p_plus*matrix[p+1,q+1]+p_0*matrix[p,q+1]+p_minus*matrix[p-1,q+1] - 4*straf
      }
      x <- x_værdier[p,q] <- x_vektor[which.max(test)]
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
      matrix[p,q] <- p_plus*matrix[p+1,q+1]+p_0*matrix[p,q+1]+p_minus*matrix[p-1,q+1]
      mod_matrix[p,q] <- p_plus*mod_matrix[p+1,q+1]+ p_0*mod_matrix[p,q+1]+p_minus*mod_matrix[p-1,q+1]
    }
  }
  return(list(XP=matrix,handling=x_værdier, mod_XP=mod_matrix))
}
# Kode svarende til metode 3 - analytisk løsning #
Strategi_funktion3 <- function(lam1,lam2,a_plus,a_minus,eps,dt) {
  matrix <- x_værdier <- mod_matrix <- matrix(0,nrow=21,ncol=1/dt+1)
  matrix[11,1/dt+1] <- mod_matrix[11,1/dt+1] <- 1+eps
  matrix[12:21,1/dt+1] <-3
  matrix[21,] <- 3
  mod_matrix[1,] <- 3
  mod_matrix[1:10,1/dt+1] <- 3
  x <- 0
  for (q in (1/dt):1) {
    for (p in 2:20) {
      if (2*matrix[p,q+1]-matrix[p-1,q+1]-matrix[p+1,q+1]<0) {
        x <- (2*matrix[p,q+1]-matrix[p-1,q+1]-matrix[p+1,q+1]) / (-abs(2*a_plus*(matrix[p-1,q+1]-matrix[p,q+1]-10^(-10))))
        if (a_plus*x^2+2*x+(lam1+lam2-1/dt)>0) {
          d <- 2^2-4*a_plus*(lam1+lam2-1/dt)
          højrerod <- (-2+sqrt(d))/(2*a_plus)
          x <- højrerod
        }
      }
      else if (2*matrix[p,q+1]-matrix[p-1,q+1]-matrix[p+1,q+1]>0){
        x <- (2*matrix[p,q+1]-matrix[p-1,q+1]-matrix[p+1,q+1]) / (-abs(2*a_minus*(matrix[p-1,q+1]-matrix[p,q+1]-10^(-10))))
        x <- max(x,-lam1)
        if (lam2+x+a_minus*x^2 < 0) {
          d <- 1-4*a_minus*lam2
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
      matrix[p,q] <- p_plus*matrix[p+1,q+1]+p_0*matrix[p,q+1]+p_minus*matrix[p-1,q+1]
      mod_matrix[p,q] <- p_plus*mod_matrix[p+1,q+1]+ p_0*mod_matrix[p,q+1]+p_minus*mod_matrix[p-1,q+1]
    }
  }
  return(list(XP=matrix,handling=x_værdier, mod_XP=mod_matrix))
}
# Tre funktioner til at undersøge, hvor tidskrævende koderne er i ms. Kørt 100 gange hver og tages gennemsnit #
tid_funktion1 <- function(lam1,lam2,a_plus,a_minus,eps,dt){
  tic()
for (i in 1:100){
  Strategi_funktion1(lam1,lam2,a_plus,a_minus,eps,dt)$XP[11,1]
}
  tid <- toc(quiet = TRUE)
  return(1000*(tid$toc-tid$tic)/100)
}

tid_funktion2 <- function(lam1,lam2,a_plus,a_minus,eps,dt){
  tic()
  for (i in 1:100){
    Strategi_funktion2(lam1,lam2,a_plus,a_minus,eps,dt)$XP[11,1]
  }
  tid <- toc(quiet = TRUE)
  return(1000*(tid$toc-tid$tic)/100)
}

tid_funktion3 <- function(lam1,lam2,a_plus,a_minus,eps,dt){
  tic()
  for (i in 1:100){
    Strategi_funktion3(lam1,lam2,a_plus,a_minus,eps,dt)$XP[11,1]
  }
  tid <- toc(quiet = TRUE)
  return(1000*(tid$toc-tid$tic)/100)
}

#Plot forventet point for både hold A og hold B givet dt hen af x-aksen - svarende til figur 3 #

tid_skridt <- seq(45,360,by= 45)
Data <- rbind(
  data.frame(skridt = tid_skridt, XP = sapply(tid_skridt, function(t) Strategi_funktion1(1.5,1.5,2,2,0,1/t)$XP[11,1]), funktion = "Optimise", gruppe="A"),
  data.frame(skridt = tid_skridt, XP = sapply(tid_skridt, function(t) Strategi_funktion2(1.5,1.5,2,2,0,1/t)$XP[11,1]), funktion = "Manuelt", gruppe="A"),
  data.frame(skridt = tid_skridt, XP = sapply(tid_skridt, function(t) Strategi_funktion3(1.5,1.5,2,2,0,1/t)$XP[11,1]), funktion = "Analytisk", gruppe="A"),
  data.frame(skridt = tid_skridt, XP = sapply(tid_skridt, function(t) Strategi_funktion1(1.5,1.5,2,2,0,1/t)$mod_XP[11,1]), funktion = "Optimise", gruppe="B"),
  data.frame(skridt = tid_skridt, XP = sapply(tid_skridt, function(t) Strategi_funktion2(1.5,1.5,2,2,0,1/t)$mod_XP[11,1]), funktion = "Manuelt", gruppe="B"),
  data.frame(skridt = tid_skridt, XP = sapply(tid_skridt, function(t) Strategi_funktion3(1.5,1.5,2,2,0,1/t)$mod_XP[11,1]), funktion = "Analytisk", gruppe="B")
)
ggplot(data=Data, aes(x=skridt,y=XP, color=funktion, group=interaction(funktion,gruppe), linetype=gruppe))+
  geom_line()+geom_point() + geom_hline(data=data.frame(XP=XP_pois(1.5,1.5)),aes(yintercept=XP), linetype="dashed")+
  labs(x="Antal skridt", y="Forventet point", color="Metode", linetype="Hold")+theme_bw()+theme(text=element_text(size=14))

#Plot over tid ud af x-aksen og tidsskridt ud af y-aksen - svarende til figur 4 #
tid_data <- rbind(
  data.frame(skridt = tid_skridt, tid = sapply(tid_skridt, function(t) tid_funktion1(1.6,1.4,2,2,0,1/t)), funktion = "Optimise"),
  data.frame(skridt = tid_skridt, tid = sapply(tid_skridt, function(t) tid_funktion2(1.6,1.4,2,2,0,1/t)), funktion = "Manuelt"),
  data.frame(skridt = tid_skridt, tid = sapply(tid_skridt, function(t) tid_funktion3(1.6,1.4,2,2,0,1/t)), funktion = "Analytisk")
)
ggplot(data=tid_data, aes(x=skridt,y=tid, color=funktion, group=funktion))+
  geom_line()+geom_point() + scale_y_log10()+
  labs(x="Antal skridt", y="Tid i ms", color="Metode")+theme_bw() + theme(text=element_text(size=14))

