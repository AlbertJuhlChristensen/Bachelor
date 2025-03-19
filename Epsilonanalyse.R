#Den optimale funktion til optimal handling #
Strategi_funktion <- function(lam1,lam2,a_plus,a_minus,eps,dt) {
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
        if (a_plus*x^2+2*x+(lam1+lam2-1/dt)>0 ) {
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
# Funktion der kan regne forventet point ud fra givne x-værdier i hver knude#
Forventet_point <- function(lam1,lam2,a_plus,a_minus,dt,handling){
  matrix <- matrix(0,nrow=21,ncol=1/dt +1)
  matrix[11,1/dt +1] <- 1
  matrix[12:21,1/dt +1] <- 3
  matrix[21,] <- 3
  for (q in (1/dt):1) {
    for (p in 2:20) {
      x <- handling[p,q]
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
      p_0 <- 1-p_plus-p_minus
      matrix[p,q] <- p_plus*matrix[p+1,q+1]+p_0*matrix[p,q+1]+p_minus*matrix[p-1,q+1]
    }}
  return(matrix)
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
# Oprettelse af forskellige sekvenser og plot svarende til figur 8 #
eps_værdi <- seq(-1,2, by=0.2)
lam1_værdier <- seq(1,3, by=0.2)
lam2_værdier <- seq(3,1, by=-0.2)
lam_værdier <- seq(-2,2, by=0.5)
data <- data.frame(eps=numeric(),XP=numeric(),lam_forskel=numeric())
for (eps in eps_værdi){
  for (lam_diff in lam_værdier) {
    lam1 <- 2 + lam_diff/2
    lam2 <- 2 - lam_diff/2
    max <- Strategi_funktion(lam1,lam2,2,2,eps,1/180)$handling
    XP_værdi <- Forventet_point(lam1,lam2,2,2,1/180,max)
    data <- rbind(data, data.frame(eps=eps, XP=XP_værdi[11,1]-XP_pois(lam1,lam2)[1], lam_forskel=lam1-lam2))
  }
}

ggplot(data=data, aes(x=eps,y=XP, color=factor(lam_forskel), group=factor(lam_forskel)))+geom_line()+
  labs(x="Epsilon", color="Forskel i lambda", y="Forskel på forventet point ift. uden model")+theme_bw() +
  theme(text=element_text(size=14))

# Analyse af den konvekse faktor oveni og plot svarende til figur 9 #
a_værdier <- seq(0,4,by=2)
lam_værdier2 <- seq(-1,1, by=1)
eps_værdi2 <- seq(-1,2,by=0.5)
data2 <- data.frame(eps=numeric(),XP=numeric(),lam_forskel=numeric(), a_forskel=numeric())
for (eps in eps_værdi2){
  for (lam_diff in lam_værdier2) {
    for (a in a_værdier) {
      a_plus <- 2.5 + a/2
      a_minus <- 2.5 - a/2
      lam1 <- 1.5 + lam_diff/2
      lam2 <- 1.5 - lam_diff/2
      max <- Strategi_funktion(lam1,lam2,a_plus,a_minus,eps,1/180)$handling
      XP_værdi <- Forventet_point(lam1,lam2,a_plus,a_minus,1/180,max)
      data2 <- rbind(data2, data.frame(eps=eps, XP=XP_værdi[11,1]-XP_pois(lam1,lam2)[1], 
                                     lam_forskel=lam1-lam2, a_forskel=a_plus-a_minus))
    }
  }
}

ggplot(data=data2, aes(x=eps,y=XP, color=factor(a_forskel), linetype=factor(lam_forskel),group=interaction(lam_forskel,a_forskel)))+
   geom_line(size=0.8)+theme_bw() + scale_linetype_manual(values = c("dotted", "solid", "twodash")) + 
  labs(x="Epsilon", y="Gevinst i forventet point", color="Forskel i a", linetype="Forskel i lambda") +
  theme(text=element_text(size=14))



