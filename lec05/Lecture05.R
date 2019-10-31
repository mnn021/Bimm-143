# Lecture 5: Data visualization 
x <- rnorm(1000)

# Some summary stats
mean(x)
sd(x)

summary(x)
boxplot(x)

hist(x)
rug(x)

# The most important measurement should exploit the highest ranked encoding possible

# Lab: Section 2 Scatterplots
# let's read out input files first 

read.table("bimm143_05_rstats/weight_chart.txt", header = TRUE)
baby <- read.table("bimm143_05_rstats/weight_chart.txt", header = TRUE)

# Plotting weight_chart.txt , Section 2A

plot(baby$Age, baby$Weight, type = "o", pch = 15, cex = 1.5, lwd = 2, ylim = c(2,10), 
     xlab = "Months", ylab = "kg", main = "Baby weight with Age")

plot(1:5, pch = 1:5, cex = 1:5)

# Barry's plot: 
plot(baby$Age, baby$Weight, typ="o", 
     pch=15, cex=1.5, lwd=2, ylim=c(2,10), 
     xlab="Age (months)", ylab="Weight (kg)", 
     main="Baby weight with age", col="blue")


# Section 2B
read.delim("bimm143_05_rstats/feature_counts.txt")

# can also do...
mouselec5ex <- read.table("bimm143_05_rstats/feature_counts.txt", sep="\t", header=TRUE)
mouselec5 <- read.delim("bimm143_05_rstats/feature_counts.txt")
barplot(mouselec5$Count)


barplot(mouselec5$Count, names.arg = mouselec5$Feature, horiz = TRUE, ylab = "", 
        main = "Number of features in the mouse GRCm38 genome", las = 1, xlim = c(0,80000))
# par(mar) controls the numerical vector of the form c(bottom, left, top, right)
# do par(mar) before putting the box plot
par()$mar
par(mar = c(3.1 , 11.1, 4.1, 2))
barplot(mouselec5$Count, names.arg = mouselec5$Feature, horiz = TRUE, ylab = "", 
        main = "Number of features in the mouse GRCm38 genome", las = 1, xlim = c(0,80000))


#Barry's 
par(mar=c(5, 11, 2, 2))
barplot(mouselec5$Count, names.arg = mouselec5$Feature, horiz = TRUE, las = 1, 
        main = "Number of features in the mouse GRCm38 genome", xlim = c(0,80000))


# Section 2C optional --> already did a histogram
x <- c(rnorm(10000),rnorm(10000)+4)
hist(x, breaks=80)


# Section 3A
read.csv2("bimm143_05_rstats/male_female_counts.txt", header = TRUE)
read.delim("bimm143_05_rstats/male_female_counts.txt", header = TRUE)

mflec5 <- read.delim("bimm143_05_rstats/male_female_counts.txt", header = TRUE)
barplot(mflec5$Count, names.arg = mflec5$Sample, col = rainbow(nrow(mflec5)), ylab = "Count", las = 2)

# or...
barplot(mflec5$Count, names.arg = mflec5$Sample, col = rainbow(10), ylab = "Count", las = 2)

#col is for columns and las changes label orientation

par(mar = c(3, 5, 4, 7))
barplot(mflec5$Count, names.arg = mflec5$Sample, col = c("blue2", "red2"), 
        ylab = "Count", las = 2)


# Section 3B
geneslec5 <- read.delim("bimm143_05_rstats/up_down_expression.txt")
nrow(geneslec5)

table(geneslec5$State)
plot(geneslec5$Condition1, geneslec5$Condition2, col=geneslec5$State, xlab = "Expression Condition 1", 
     ylab = "Expression Condition 2")
palette()

# Run levels() on the State column. 
levels(geneslec5$State)

palette(c("blue", "grey", "red"))
plot(geneslec5$Condition1, geneslec5$Condition2, col=geneslec5$State, xlab = "Expression Condition 1", 
     ylab = "Expression Condition 2")


# Section 3C plotting expression vs gene regulation
read.delim("bimm143_05_rstats/expression_methylation.txt")
methlec5 <- read.delim("bimm143_05_rstats/expression_methylation.txt")
nrow(methlec5)
plot(methlec5$gene.meth, methlec5$expression)

densCols(methlec5$gene.meth, methlec5$expression)
dcols <- densCols(methlec5$gene.meth, methlec5$expression)

# Plot changing the plot character ('pch') to a solid character
plot(methlec5$gene.meth, methlec5$expression, col = dcols, pch = 20)

# Find the indices of genese with exp. levels above 0
inds <- methlec5$expression > 0

# Plot just these points
plot(methlec5$gene.meth[inds], methlec5$expression[inds])

## Make a color vector for these genes and plot
dcols <- densCols(methlec5$gene.meth[inds], methlec5$expression[inds])
plot(methlec5$gene.meth[inds], methlec5$expression[inds], col = dcols, pch = 20)

dcols.custom <- densCols(methlec5$gene.meth[inds], methlec5$expression[inds],
                         colramp = colorRampPalette(c("blue2",
                                                      "green2",
                                                      "red2",
                                                      "yellow")) )

plot(methlec5$gene.meth[inds], methlec5$expression[inds], 
     col = dcols.custom, pch = 20)
