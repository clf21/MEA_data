
# plotting individual treatment
gp <- ggplot(data=subset(all_data, trt=="Acetaminophen"), aes(x=DIV, y=mi, group=dose, color=dose)) + geom_point() + stat_summary(size=1.2, geom="line") + scale_colour_gradient(trans="log", low="blue", high="red")
gp


gp <- ggplot(data=subset(all_data, trt=="Haloperidol"), aes(x=DIV, y=r, group=dose, color=dose)) + geom_point() + stat_summary(size=1.2, geom="line") + scale_colour_gradient(trans="log", low="blue", high="red")
gp

gp <- ggplot(data=subset(all_data, trt=="Haloperidol"), aes(x=DIV, y=mi, group=dose, color=dose)) + geom_point() + stat_summary(size=1.2, geom="line") + scale_colour_gradient(trans="log", low="blue", high="red")
gp


gp <- ggplot(data=subset(all_data, trt=="Sodium Orthovanadate"), aes(x=DIV, y=r, group=dose, color=dose)) + geom_point() + stat_summary(size=1.2, geom="line") + scale_colour_gradient(trans="log", low="blue", high="red")
gp

gp <- ggplot(data=subset(all_data, trt=="Sodium Orthovanadate"), aes(x=DIV, y=mi, group=dose, color=dose)) + geom_point() + stat_summary(size=1.2, geom="line") + scale_colour_gradient(trans="log", low="blue", high="red")
gp

