df = read.csv("NonRBC.csv")

df = df[which(df$ancestry !=  "EastAsian"), ]
df = df[which(df$ancestry !=  "African"), ]
#df = df[which(df$ancestry !=  "Yamnaya"), ]

df = df[which(df[[1]] !=  "Elevated_Blood_Glucose"), ]

for (i in 1:nrow(df)) {
  temp = df[i,1]
  temp = gsub("_indbootstrap", "", temp)
  df[i,1] = gsub("_", " ", temp)
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
                      col =  c(  rep(c(2,1), floor(ntraits/2)    )    )   )
# Plot
ymin = min( df$mean - df$yerr1) 
ymax =   max( df$mean + df$yerr2)

names(df)[1] = "Phenotype"


vv = c()
for (i in 1:nrow(df)) {
  l = df$mean - df$yerr1
  h = df$mean + df$yerr2
  if (h[i] < 0 |  l[i] > 0) {
    print(df[i,])
    vv = c(vv, i)
  }
}

stars_df <- df[vv, c("Phenotype", "Ancestry", "mean", "yerr2")]
stars_df$label <- "*"
stars_df$y_star <- stars_df$mean    # Add vertical offset for star

dodge <- position_dodge(width = 0.87)

 p = ggplot(df, aes(x = Phenotype, y = mean, color = Ancestry)) +
  geom_rect(data=shading, aes(xmin=xstart,xmax=xend, 
                              ymin = ymin -10,   ymax  = ymax  +!0), 
            fill = c(rep(c("grey","white"), floor(ntraits / 2)  )  )  ,    alpha = 0.2,    inherit.aes = FALSE)  + geom_hline(yintercept = 0, linetype="dashed", alpha=0.5) + 
  geom_point(size=2, position=position_dodge(width=0.87)) +  geom_errorbar(aes(ymin = mean - yerr1, ymax = mean + yerr2),  position = position_dodge(width = 0.87), width = 0.5) +
  theme_classic() +  theme(axis.text = element_text(size=10),  
                           axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1) , plot.title = element_text( hjust = 0.5)) +
  coord_cartesian(ylim = c(ymin, ymax)) +  labs(y = "Ancestral Risk Score (z-score)") +
  scale_color_manual(values=c( "#56B4E9", "#009E73", "#E69F00",  
                                 "royalblue4", "#CC79A7"))

 # 2. Create a new data frame for your specific points
 special_points <- data.frame(
   x = c(5-0.345, 9, 11-0.345, 15-0.345, 17+0.345),
  y= c( -1.923418-0.1 ,-2.236729 -0.1 ,-2.067708 -0.1 ,-1.597757-0.1  ,2.276339+0.1)
 )
 
 # 3. Add the special points on top
 p = p + geom_point(data = special_points, aes(x = x, y = y), 
                color = "black", size = c(1.2,1.2,1.2,1.2,1.4),fill="black", shape = c(25,25,25,25,17))  # triangle shape

pdf(("RiskScore_bases.pdf"), width = 11.5, height = 5)
print(p)
dev.off()

