#Plot_loss

df <- data.frame(Case = c("A","B","C"),
                 PL = c(0.551,0.545,0.544),
                 OL = c(1.125,1.237,1.232))




ggplot(df,aes(x=Case,y=PL)) +
  geom_segment( aes(x=Case, xend=Case, y=0, yend=PL),size=1.5,color="#1E90FF") +
  #  geom_point( aes(x=Color, y=lw), color="#FF4500", size=6 ) +
  geom_point( aes(x=Case, y=PL), color="#FF4500", size=4 ) +
  #  geom_errorbar(aes(ymin=lw, ymax=up), size=1, width=.25,color="#965C85") +
  coord_flip() +
  ylab("variable importance")  +
  xlab("") +
  theme_bw() +
  theme( axis.text.y   = element_text(size=12))
