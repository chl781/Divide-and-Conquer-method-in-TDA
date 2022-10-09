# Plot for all data in the PD for 2 data sets

PlotRepeat1 <- function(S1,S2){
  S=matrix(0,nrow=nrow(S1)+nrow(S2),
           ncol=3)
  S[,1:2]=rbind(S1,S2)
  S<-data.frame(S)
  colnames(S) = c('X','Y','index')
  S[,3]=c(rep('true',nrow(S1)),rep('Combined',nrow(S2)))
  linecolors <- c("#714C02", "#01587A")
  fillcolors <- c("#9D6C06", "#077DAA")
  
  # partially transparent points by setting `alpha = 0.5`
  ggplot(S, aes(X, Y, colour=index, fill = index)) +
    geom_point(position=position_jitter(h=0.001, w=0.01),
               shape = 21, alpha = 0.5, size = 3) +
    scale_color_manual(values=linecolors) +
    scale_fill_manual(values=fillcolors) +
    theme_bw()+
    geom_abline(slope = 1,intercept = 0)+
    xlim(0,2*max(S[,1]))+
    ylim(0,2*max(S[,2]))
}