#Load packages
library(ggplot2)
library(grid)

args = commandArgs(TRUE)
collate.qc = args[1] # 'qc/collate.qc.tsv'

#Read in data
Quals <- read.delim(collate.qc, header=TRUE, stringsAsFactors=TRUE, row.names="Name")

previous_theme <- theme_set(theme_bw())

#make plots
a <- ggplot(data=Quals, aes(x=PCT_MRNA_BASES))
a <- a + geom_histogram() + xlab(expression("mRNA Base %")) +  ylab(NULL)
a <- a + theme(
  axis.line=element_line(colour="black", size=0.25),
  axis.title.x = element_text(size=8),
  axis.text.y = element_text(size=8),
  axis.text.x = element_text(size=8, angle=45,vjust=.5),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.background = element_rect(colour = NA),
  panel.border=element_blank())
#a



b <- ggplot(data=Quals, aes(x=PCT_MAPPED))
b <- b + geom_histogram() + xlab(expression("% Mapped")) +  ylab(NULL)
b <- b + theme(
  axis.line=element_line(colour="black", size=0.25),
  axis.title.x = element_text(size=8),
  axis.text.y = element_text(size=8),
  axis.text.x = element_text(size=8, angle=45,vjust=.5),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.background = element_rect(colour = NA),
  panel.border=element_blank())
#b

c <- ggplot(data=Quals, aes(x=PCT_INTRONIC_BASES))
c <- c + geom_histogram() + xlab(expression("Intronic Base %")) +  ylab(NULL)
c <- c + theme(
  axis.line=element_line(colour="black", size=0.25),
  axis.title.x = element_text(size=8),
  axis.text.y = element_text(size=8),
  axis.text.x = element_text(size=8, angle=45,vjust=.5),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.background = element_rect(colour = NA),
  panel.border=element_blank())
#c


d <- ggplot(data=Quals, aes(x=MEDIAN_3PRIME_BIAS))
d <- d + geom_histogram() + xlab(expression("Median 3’UTR Bias")) +  ylab(NULL)
d <- d + theme(
  axis.line=element_line(colour="black", size=0.25),
  axis.title.x = element_text(size=8),
  axis.text.y = element_text(size=8),
  axis.text.x = element_text(size=8, angle=45,vjust=.5),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.background = element_rect(colour = NA),
  panel.border=element_blank())
#d

e <- ggplot(data=Quals, aes(x=MEDIAN_5PRIME_BIAS))
e <- e + geom_histogram() + xlab(expression("Median 5’UTR Bias")) +  ylab(NULL)
e <- e + theme(
  axis.line=element_line(colour="black", size=0.25),
  axis.title.x = element_text(size=8),
  axis.text.y = element_text(size=8),
  axis.text.x = element_text(size=8, angle=45,vjust=.5),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.background = element_rect(colour = NA),
  panel.border=element_blank())
#e

f <- ggplot(data=Quals, aes(x=PCT_USABLE_BASES))
f <- f + geom_histogram() + xlab(expression("Usable Base %")) +  ylab(NULL)
f <- f + theme(
  axis.line=element_line(colour="black", size=0.25),
  axis.title.x = element_text(size=8),
  axis.text.y = element_text(size=8),
  axis.text.x = element_text(size=8, angle=45,vjust=.5),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.background = element_rect(colour = NA),
  panel.border=element_blank())
#f

#output graphs
jpeg("qc/qc_histogram.jpg", quality=100, height= 2.8 , width=8, units= "in", res=300, pointsize=8)
grid.newpage()
pushViewport(viewport(layout=grid.layout(1,6)))
vplayout <- function(x,y)
  viewport(layout.pos.row=x, layout.pos.col=y)
print(b, vp = vplayout(1,1))
print(a, vp = vplayout(1,2))
print(c, vp = vplayout(1,3))
print(f, vp = vplayout(1,4))
print(e, vp = vplayout(1,5))
print(d, vp = vplayout(1,6))
dev.off()
