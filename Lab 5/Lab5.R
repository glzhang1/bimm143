#' ---
#' title: "BIMM 143 Class 05 F graphics intro"
#' author: "Gloria Zhang"
#' date: "Jan 24th, 2019"
#' ---



#boxplot
x <- rnorm(1000,0)
boxplot(x)

boxplot (x,horizontal = TRUE)

#hands on session 2
#line plot
read.table("bimm143_05_rstats/weight_chart.txt",header= TRUE)
weight <- read.table("bimm143_05_rstats/weight_chart.txt",header= TRUE)
plot(weight ,pch=15, cex=1.5, lwd=2, ylim=c(2,10), xlab="Age(months)", ylab="Weight (kg)", main="Baby Weight with Age")
#View(weight)
plot(weight[,1], weight[,2],typ="o")
plot(weight$Age, weight$Weight, typ= "o",pch=15, cex=1.5, lwd=2, ylim=c(2,10), xlab="Age(months)", ylab="Weight (kg)", main="Baby Weight with Age",col=c("blue","red"))
plot(weight$Age, weight$Weight, typ= "o",pch=15, cex=1:3, lwd=2, ylim=c(2,10), xlab="Age(months)", ylab="Weight (kg)", main="Baby Weight with Age",col=c("blue","red"))

# barplot
barplot(VADeaths, beside = TRUE)
 # Input our feature count data
mouse <- read.table("bimm143_05_rstats/feature_counts.txt", header= TRUE, sep="\t")
#View(feat)
#?barplot
barplot(mouse$Count, horiz= TRUE, names.arg= mouse$Feature, las=2)

#change margin so we can see labels
# default parameters
par()$mar
par(mar= c(5.1, 11.1, 4.1, 2.1))
barplot(mouse$Count, horiz= TRUE, names.arg= mouse$Feature, las=2)

#colors
barplot(mouse$Count, horiz= TRUE, names.arg= mouse$Feature, las=2, col=cm.colors(5))
mf <- read.table("bimm143_05_rstats/male_female_counts.txt", header= TRUE, sep="\t")
#View(mf)
barplot(mf$Count, names.arg= mf$Sample, las=2, col= c("blue1","red2"))
par(mar= c(3.1, 6.1, 4.1, 2.1))

#expression data
e <- read.table("bimm143_05_rstats/up_down_expression.txt", header= TRUE, sep="\t")
#View(e)

 #how many genes
nrow(e)

#how many up, down and all around?
table(e$State)
plot(e$Condition1,e$Condition2, col=e$State)

#palette
palette(c("red","light gray","blue"))
plot(e$Condition1,e$Condition2, col=e$State)