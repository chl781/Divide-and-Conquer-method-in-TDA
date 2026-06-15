GetDiagram <- function(xdiag0,xlab0="Birth",ylab0="Death",labels=c(expression(H[0]), expression(H[1]))){
  xdiag <- tibble(Homology=as.character(xdiag0[,1]),
                  birth=xdiag0[,2],
                  death=xdiag0[,3])
  ggplot(xdiag, aes(x=birth,y=death,color=Homology,shape=Homology)) +
    geom_point(size=4) +
    geom_abline(slope=1,intercept=0) +
    xlab(xlab0) + ylab(ylab0) +
    scale_color_discrete(name = NULL, labels = labels) +
    scale_shape_discrete(name = NULL, labels = labels) +
    theme_bw() +
    theme(text = element_text(size=20),
          legend.title = element_blank(),
          legend.position=c(.75,.4),
          legend.background = element_blank(),
          legend.box.background = element_rect(color = "black")) +
    xlim(range(xdiag0[,2:3])) +
    ylim(range(xdiag0[,2:3])) +
    coord_fixed() #fixes the aspect ratio
}