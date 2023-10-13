my_packages <- c('ggplot2', 'popbio')
new_packages <- my_packages[!(my_packages %in% installed.packages()[,'Package'])]
if(length(new_packages)) install.packages(new_packages, dependencies = T)

library(ggplot2)
library(popbio)

#Juvenile survival: 0.463 (95% CI 0.404–0.524)
#Yearling survival: 0.510 (95% CI 0.445–0.574)
#Adult 2+ survival: 0.559 (95% CI 0.499–0.618)


# first enter the values into a dataframe 
survival <- data.frame(stage=factor(c('Juvenile','Yearling','Adult'), levels=c('Juvenile','Yearling','Adult')), estimate=c(0.463, 0.510, 0.559), lcl=c(0.404, 0.445, 0.499), ucl=c(0.524, 0.574, 0.618))

# then plot by stage
ggplot(survival, aes(stage, estimate, ymin=lcl, ymax=ucl)) + 
  geom_errorbar(width=0.2) + geom_point() + ylim(0,1)

nestdata <- read.table("~/Downloads/MPM practical/gjeroynest.txt", header = TRUE, sep = '\t')
head(nestdata)

ClutchNo <- mean(nestdata$clutchno)
HatchingSuc <- mean(nestdata$hatchingsuc)
FledglingNo <- mean(nestdata$chickno)

(ClutchNo * HatchingSuc * FledglingNo) / 2

# save our estimates of the vital rates
R <- (ClutchNo * HatchingSuc * FledglingNo) / 2
Phi.juv <- survival$estimate[survival$stage=='Juvenile'] 
Phi.yr <- survival$estimate[survival$stage=='Yearling'] 
Phi.ad <- survival$estimate[survival$stage=='Adult'] 

# Juvenile to Juvenile: Phi.juv * R
# Yearling to Juvenile: Phi.yr * R
# Adult to Juvenile: Phi.ad * R
# Juvenile to Yearling: Phi.juv
# Yearling to Yearling: 0 
# Adult to Yearling: 0 
# Juvenile to Adult: 0
# Yearling to Adult: Phi.yr
# Adult to Adult: Phi.ad

# put the transition probabilities into a vector 
sparrowMPM <- c(Phi.juv * R, Phi.yr * R, Phi.ad * R, Phi.juv, 0, 0, 0, Phi.yr, Phi.ad)

# save that vector as a matrix, specifying the number of rows and columns
# use the byrow=TRUE argument to tell R that the first the elements of the vector correspond to the first row of the matrix 
sparrowMPM <- matrix(sparrowMPM, nrow=3, ncol=3, byrow=T)
sparrowMPM

lambda(sparrowMPM)

# project over 15 years
t <- 15
# start with 50 juveniles, 20 yearlings and 30 adults
n0 <- c(50,20,30)

# project dynamics 
projection <- pop.projection(sparrowMPM, n0, iterations = t)
projected <- data.frame(time=1:15, N=projection$pop.sizes)

# plot projected pop size over time
ggplot(projected, aes(time, N)) + 
  geom_line() + ylim(0,150) + ylab('Projected N')

popest <- read.table("~/Downloads/MPM practical/popest.txt", header = TRUE, sep = '\t')
head(popest)

# plot N over time
ggplot(popest, aes(year, N)) + 
  geom_line() + ylim(0,200) + ylab('Observed N')

stages <- c('Juv','Yr','Ad')
colnames(sparrowMPM) <- stages
rownames(sparrowMPM) <- stages

stable.stage(sparrowMPM)
reproductive.value(sparrowMPM)

# list the vital rates
sparrow.param <- list(Phi.juv = Phi.juv, Phi.yr = Phi.yr, Phi.ad = Phi.ad, R = R)

# give the matrix equation 
sparrow.equation <- expression(Phi.juv * R, Phi.yr * R, Phi.ad * R, Phi.juv, 0, 0, 0, Phi.yr, Phi.ad)

# run the sensitivity analysis
sens <- vitalsens(sparrow.equation, sparrow.param)
sens

# plot elasticity of the vital rates 
sens$vitalrate <- factor(c('Phi.juv', 'Phi.yr', 'Phi.ad', 'R'), levels = c('Phi.juv', 'Phi.yr', 'Phi.ad', 'R'))
ggplot(sens, aes(vitalrate, elasticity)) + 
  geom_bar(stat = 'identity') 
