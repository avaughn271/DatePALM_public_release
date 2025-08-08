df = read.csv("RBC.csv")

df = df[which(df$ancestry !=  "EastAsian"), ]
df = df[which(df$ancestry !=  "African"), ]
#df = df[which(df$ancestry !=  "Yamnaya"), ]

for (i in 1:nrow(df)) {
  temp = df[i,1]
  temp = gsub("_indbootstrap", "", temp)
  df[i,1] =  ""
}
for (i in 1:nrow(df)) {
  if (df$ancestry[i] == "Yamnaya") { df$ancestry[i] = "Steppe" }
  temp = df[i,1]
  temp = gsub("_indbootstrap", "", temp)
  df[i,1] = gsub("_", " ", temp)
}
#df = df[df$ancestry != "Steppe", ]

#df = df[which(df$ancestry !=  "EastAsian"), ]
names(df)[2] = "Ancestry"
library(ggplot2)
ntraits = length(unique(df$phenotype))

shading <- data.frame(xstart = seq(0.5,ntraits + 0.5 - 1,1),  
                      xend = seq(1.5,ntraits + 1.5 - 1,1), 
                      col =  c(  2   )    )   
# Plot
ymin = min( df$mean - df$yerr1) 
ymax =   max( df$mean + df$yerr2)

names(df)[1] = "Phenotype"

 
dodge <- position_dodge(width = 0.87)

 p = ggplot(df, aes(x = Phenotype, y = mean, color = Ancestry)) +
  geom_rect(data=shading, aes(xmin=xstart,xmax=xend, 
                              ymin = ymin -10,   ymax  = ymax  +!0), 
            fill = c(  "white"  )  ,    alpha = 0.2,    inherit.aes = FALSE)  + geom_hline(yintercept = 0, linetype="dashed", alpha=0.5) + 
  geom_point(size=2, position=position_dodge(width=0.87)) +  geom_errorbar(aes(ymin = mean - yerr1, ymax = mean + yerr2),  position = position_dodge(width = 0.87), width = 0.5) +
  theme_classic() +  theme(axis.text = element_text(size=10),  
                           axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1) , plot.title = element_text( hjust = 0.5)) +
  coord_cartesian(ylim = c(ymin, ymax)) +  labs(y = "Ancestral Risk Score (z-score)") +
  scale_color_manual(values=c( "#56B4E9", "#009E73", "#E69F00",  
                                 "royalblue4", "#CC79A7")) + labs(x = "Red Blood Cell Count")+ theme(legend.position="none")

 
pdf(("RiskScore_RBC.pdf"), width = 3, height = 8)
print(p)
dev.off()

