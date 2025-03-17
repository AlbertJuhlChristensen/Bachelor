library(ggplot2)
#Den optimalestrategifunktion #
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
# Funktion der kan plotte 3 grafer i Facet-wrap for forskellige tidspunker#
graf_funktion_X <- function(minutter,målforskel,lam1,lam2,a_plus,a_minus,eps1,eps2,eps3,dt) {
  dataliste <- lapply(minutter, function(minut) {
    data.frame(
      forskel = rep((-målforskel):målforskel, 3),
      XP = c(
        Strategi_funktion(lam1, lam2, a_plus, a_minus, eps1, dt)$handling[(11 - målforskel):(11 + målforskel), minut],
        Strategi_funktion(lam1, lam2, a_plus, a_minus, eps2, dt)$handling[(11 - målforskel):(11 + målforskel), minut],
        Strategi_funktion(lam1, lam2, a_plus, a_minus, eps3, dt)$handling[(11 - målforskel):(11 + målforskel), minut]
      ),
      eps = factor(rep(c(eps1, eps2, eps3), each = 2 * målforskel + 1)),
      minut = as.factor(minut)
    )
  })
  data <- do.call(rbind,dataliste)
  minut_titel <- function(x) {
    paste("Minut",2*x-2)
  }
  graf <- ggplot(data, aes(y=forskel,x=XP, color=eps)) + geom_line() + geom_point() + 
    facet_wrap(~minut,ncol=3,labeller=labeller(minut= function(x) paste("Minut",as.numeric(x)/2)),  scales="free_x") + 
    labs(y = "Målforskel", x = "Handling", color = "Epsilon")+
    theme_bw() + theme(legend.position = "right") + theme(text=element_text(size=14))
  return(graf)
}
# Kørsel af denne vil give plot svarende til figur 5 #
graf_funktion_X(c(2*1,2*45,2*88),2,1.6,1.4,2,2,-1,0,1, 1/180)

#Oprettelse af sekvens for lambda og plot svarende til figur 6 #
dt <- 1/180
lam_sekvens <- seq(-2,2, by=0.01)

data <- data.frame(lam_forskel = as.numeric(), XP_diff = as.numeric(), a_værdier= as.numeric())
for (a in c(0.5,1,2,4,8)){
  for (lam_diff in lam_sekvens){
    lam1 <- 1.5 + lam_diff/2
    lam2 <- 1.5 - lam_diff/2
    XP_model <- Strategi_funktion(lam1,lam2,a,a,0,dt)$XP[11,1]
    XP_0 <- XP_pois(lam1,lam2)[1]
    data <- rbind(data, data.frame(lam_forskel = lam1-lam2, XP_diff=XP_model-XP_0, a_værdier= a))
    }
  }
ggplot(data=data, aes(x=lam_forskel, y=XP_diff, color=as.factor(a_værdier))) + 
  geom_line(size=0.8) + labs(x="Forskel i lambda", y= "Gevinst i forventet point", color="Værdi af a")+
  theme_bw() + theme(text=element_text(size=14))

# Oprettelse af sekvens af a, og plot svarende til figur 7 #
data2 <- data.frame(lam_forskel = as.numeric(), XP_diff = as.numeric(), a_forskel= as.numeric())
lam_sekvens2 <- seq(-2,2, by=0.5)
a_sekvens <- seq(-4,4, by=0.01)
for (a_diff in a_sekvens){
  a_plus <- 2.5 + a_diff/2
  a_minus <- 2.5 - a_diff/2
  for (lam_diff in lam_sekvens2){
    lam1 <- 1.5 + lam_diff/2
    lam2 <- 1.5 - lam_diff/2
    XP_model <- Strategi_funktion(lam1,lam2,a_plus,a_minus,0,dt)$XP[11,1]
    XP_0 <- XP_pois(lam1,lam2)[1]
    data2 <- rbind(data2, data.frame(lam_forskel = lam1-lam2, XP_diff=XP_model-XP_0, a_forskel=a_plus-a_minus))
  }
}

ggplot(data=data2, aes(x=a_forskel, y=XP_diff, color=as.factor(lam_forskel)))+
  geom_line(size=0.8)+ labs(x="Forskel i a", y="Gevinst i forventet point", color="Forskel i lambda")+
  theme_bw()+ theme(text=element_text(size=14))

