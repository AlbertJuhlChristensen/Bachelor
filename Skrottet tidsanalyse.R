#Skrottet tidsanalyse
lambda_værdier <- list(
  c(1.1,1.9),
  c(1.3,1.7),
  c(1.5,1.5),
  c(1.7,1.3),
  c(1.9,1.1)
)
data1<- data2<- data3 <- data.frame(tid=numeric(), Forventet_point=numeric(),lam_forskel=numeric())
for (lam in lambda_værdier) {
  tid1<- tid_funktion1(lam[1],lam[2],2,2,0)
  tid2 <- tid_funktion2(lam[1],lam[2],2,2,0)
  tid3 <- tid_funktion3(lam[1],lam[2],2,2,0)
  XP1 <- Strategi_funktion1(lam[1],lam[2],2,2,0)$XP[11,1]
  XP2 <- Strategi_funktion2(lam[1],lam[2],2,2,0)$XP[11,1]
  XP3 <- Strategi_funktion3(lam[1],lam[2],2,2,0)$XP[11,1]
  data1 <- rbind(data1,data.frame(tid = tid1, Forventet_point = XP1, lam_forskel = lam[1]-lam[2]))
  data2 <- rbind(data2,data.frame(tid = tid2, Forventet_point = XP2, lam_forskel = lam[1]-lam[2]))
  data3 <- rbind(data3,data.frame(tid = tid3, Forventet_point = XP3, lam_forskel = lam[1]-lam[2]))
}
Uden_model <- data.frame(XP=numeric(),lam_forskel=numeric())
for (lam in lambda_værdier){
  linje <- XP_pois(lam[1],lam[2])[1]
  Uden_model <- rbind(Uden_model,data.frame(lam_forskel=lam[1]-lam[2],XP=linje))
}
data1$ID <- "Optimise"
data2$ID <- "Manuelt"
data3$ID <- "Analytisk"
samlet_data <- rbind(data1,data2,data3)
ggplot(samlet_data, aes(x=tid, y=Forventet_point, color=factor(lam_forskel), shape=ID)) + 
  geom_point()+
  geom_hline(data=Uden_model, aes(yintercept = XP, color=factor(lam_forskel)),linetype="dashed")+
  labs(x="Tid i ms", y="Forventet point",color="Forskel i lam", shape="Metode")+ggtitle("Tidsanalyse") + theme_bw()
XP_pois(1.5,1.5)
