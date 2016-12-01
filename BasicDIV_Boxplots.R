
all_data[all_data[,"dose"]==0,"trt"] <- "Control"
cntrl_data <- subset(all_data, trt=="Control")

box_param <- function(param, name) {
  pdf(paste(name, "cntrl_boxplots.pdf", sep="_"))
  boxplot(get(param) ~ DIV, data=cntrl_data, boxwex=0.5, outline=FALSE, ylab=name, xlab="DIV", col=c("lightcyan4","darkslategrey","darkslategray4","darkslategray3","darkslategray2"), notch=TRUE, cex.lab=1.6, cex.axis=1.4, las=1)
  dev.off()
}

box_param("burst.per.min", "Burst Rate")
box_param("meanfiringrate", "Mean firing rate")
box_param("r", "Mean correlation")
box_param("ns.percent.of.spikes.in.ns", "% in network spikes")
box_param("mi", "Mutual information")
box_param("per.spikes.in.burst", "% in bursts")