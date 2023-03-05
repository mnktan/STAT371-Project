### STAT371 Project ###

library(car)
library(sandwich)
library(lmtest)
library(MASS)
install.packages("SignifReg")
library(SignifReg)
options(scipen=999999)
options(scipen=1)
concrete <- read.csv("Concrete-Data.csv")

### part 1 ###
# Choosing variables to be used for the model: do independence test and
# exclude dependent variables
concrete_mod <- lm(ConcreteCompressiveStrength ~ Cement + BlastFurnaceSlag + FlyAsh +
                                                 Water + Superplasticizer + CoarseAggregate +
                                                 FineAggregate + Age, data = concrete)
summary(concrete_mod)
linearHypothesis(concrete_mod, c("CoarseAggregate = 0", "FineAggregate = 0"))
coeftest(concrete_mod, vcov = vcovHC(concrete_mod, "HC1"))

# Diagnostics (4-plot), leverage, using robust regression

# 4-plot function from slides 15 case study 
diagnostics_plots <- function(model, data, index) { 
  n <- nrow(data)
  par(mfrow = c(2,2))
  # quantile-quantile plot
  plot(model, which = 2, pch = 1, cex = 0.55, col = "dark green", cex.main = 0.8)
  
  # studentized residual vs index 
  plot(x = index, y = rstudent(model),
       main = "Studentized residuals vs index", xlab ="Index",
       ylab = "Studentized residuals", pch = 1, cex = 0.55, col = "dark green", cex.main = 0.8)
  
  # Residuals vs predicted  
  plot(x = fitted.values(model), y = resid(model),
       main = "Residuals vs predicted values", 
       xlab = "Fitted values", ylab = "Residuals", pch = 1, cex = 0.55, col = "dark green", cex.main = 0.8)
  
  # Index vs leverage
  plot(x = index, y = hatvalues(model), xlab = "Index", ylab = "Leverage",
       main = "Leverage vs index plot", pch = 1, cex = 0.55, col = "dark green", cex.main = 0.8)
}

# index from 1-1030
cindex <- seq(1, 1030)

diagnostics_plots(concrete_mod, data = concrete, index = cindex) 

# quantile-quantile plot ONLY
plot(concrete_mod, which = 2, pch = 19, cex = 0.8)

# find large leverage points and residuals
# largest leverage is on index 67
max_lev <- hatvalues(concrete_mod)[hatvalues(concrete_mod) == max(hatvalues(concrete_mod))]

# get two largest residuals
# indexes respectively (382, 384)
con_res <- sort(resid(concrete_mod))
con_res_lar1 <- con_res[1030]
con_res_lar2 <- con_res[1029]


### part 3 ###
# interaction
# check for interaction between coarse, age, etc. and every other variables
# the interaction var will be product of two cont vars
# so ex: coarsecement = coarse * cement

### function to test ###
inter_test <- function(int1, int2) {
  inter_var <- concrete[int1][,] * concrete[int2][,]
  inter_dat <- data.frame(concrete, inter_var)
  inter_mod <- lm(ConcreteCompressiveStrength ~ Cement + BlastFurnaceSlag + FlyAsh +
                            Water + Superplasticizer + CoarseAggregate +
                            FineAggregate + Age + inter_var, 
                          data = inter_dat)
  #summary(inter_mod)
  coeftest(inter_mod, vcov = vcovHC(inter_mod, "HC1"))
}

# top 3 inters
inter_test("Superplasticizer", "Age")
inter_test("Water", "Age")
inter_test("FlyAsh", "Age")


# beat 0.6142
# high sig: cement, blast, fly, water, super, age
# highest: flyash and age, water and age, super and age (HIGHEST)
# cement, blastfurnaceslag, flyash, water, superplasticizer, coarseaggregate
# fineaggregate, age

concrete <- read.csv("Concrete-Data.csv")

# rho value for contour plot
rho <- function(age, super) {
  2.2658697 + (0.0876651*age) + (-0.3834723*super) + (0.0194690*age*super)
}

xage <- seq(min(concrete$Age), max(concrete$Age), length = 1000)
xsup <- seq(0, max(concrete$Superplasticizer), length = 1000)
Rho <- outer(xage, xsup, "rho")

image(xage, xsup, Rho, col = heat.colors(1000), main = "Contour Plot of function", 
      xlab = "Age", ylab = "SuperPlasticizer")
contour(xage, xsup, Rho, add = T, levels = c(0,10,20,30,40,50,60,70,80,90,100))


# other contour plot
install.packages("plotly")
library(plotly)
plot.new()
plot_ly(x = xage, y = xsup, Rho, type = "contour")



### IF INTERACTION EXISTS DO THE FOLLOWING TO INTERPRET###
# 1) do contour plot where all variables not involved are 0
#    so only the two variables and interaction variable should be
#    regressed. See how much y changes.
# 2) (OPTIONAL) See how much R2 increases from adding the interaction variable
# 3) add interaction to model

### part 5 ###
# robust regression
# model selection 
# we will present the results of forward and backward selection
independent_vars <- c("Cement", "BlastFurnaceSlag", "FlyAsh",
                      "Water", "Superplasticizer", "CoarseAggregate",
                      "FineAggregate", "Age") 

independent_vars_wint <- c("Cement", "BlastFurnaceSlag", "FlyAsh",
                      "Water", "Superplasticizer", "CoarseAggregate",
                      "FineAggregate", "Age","asuper") 


# contains indp vars, edit them if change is needed
asuper <- concrete$Age * concrete$Superplasticizer
concrete_int_mod <- data.frame(concrete, asuper)
concrete_ind <- concrete_int_mod[1:1030,c("ConcreteCompressiveStrength",independent_vars_wint)]


# non-interaction model
model.concrete.basic<-lm(ConcreteCompressiveStrength~1,data=concrete_ind) # Null model 
model.concrete.full<-lm(ConcreteCompressiveStrength~.,data=concrete_ind) # Full model 

# forward selection using AIC
select.forward.AIC <- stepAIC(model.concrete.basic, 
                              scope=formula(model.concrete.full), 
                              direction="forward") 
select.forward.AIC$anova

# backward selection using AIC
select.backward.AIC <- stepAIC(model.concrete.full, 
                               scope=formula(model.concrete.full), 
                               direction="backward") 
select.backward.AIC$anova


#### !!!!ROBUST REGRESSION!!!! ####
concrete_mod_wint <- lm(ConcreteCompressiveStrength ~ asuper + Cement + BlastFurnaceSlag + 
                          FlyAsh + Water + Superplasticizer + CoarseAggregate +
                          FineAggregate + Age, 
                        data = concrete_ind)
coeftest(concrete_mod_wint, vcov = vcovHC(concrete_mod_wint, "HC1"))
#### !!!!ROBUST REGRESSION!!!! ####



concrete_new_mod <- lm(ConcreteCompressiveStrength ~ asuper + Cement + BlastFurnaceSlag + 
                       FlyAsh + Water + Superplasticizer + Age, 
                       data = concrete_ind)
summary(concrete_new_mod)
coeftest(concrete_new_mod)

diagnostics_plots(concrete_new_mod, data = concrete, index = cindex) 

# 99% confidence interval for estimate of Superplasticizer on 1022 DF
cval <- qt(0.995, 1022)
cval95 <- qt(0.975, 1022)

super_99ci <- c(-0.4068613 - cval*0.0876906, -0.4068613 + cval*0.0876906)
super_99ci

super_95ci <- c(-0.4068613 - cval95*0.0876906, -0.4068613 + cval95*0.0876906)
super_95ci
