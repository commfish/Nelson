# notes ----
# Nelson Alaska sockeye forecast
# author's contact: sarah.power@alaska.gov
# date November 2019

# load ----
library(tidyverse)  # for data manipulation
library(fpp2) #for time series
library(RDS) #for time series prediction
library(ggrepel) #for graphing
library(ggpubr)#for graphing
library(grid)#for graphing
library(broom)#for cleaning up data, used in predction
library(caret)#used for cross validation 
library(Metrics)#used for cross validation 
#source('code/functions.R')
options(scipen = 999)
set.seed(100) # for reproducible results

# data ----
# oage = ocean age = years in the ocean
#spawnerstm4 spawners that produced that year class.
#run = the fish that returned  in either the harvest or the escapement for that river system
data <- read_csv('data/nelson.csv') %>%
  filter(year > 2002) #%>%

# Only the most recent years should have NAs
data[!complete.cases(data),]
tail(data,13)

d_spawn <- data %>% 
  filter(year > 2002) %>%
  select(year, spawnerstm4, oage_2)

# Only the most recent years should have NAs
d_spawn[!complete.cases(d_spawn),] 

(last_yr <- max(data$year, na.rm =TRUE))

tail(d_spawn)
# check this particularly that this_yr-2 is capturing the right data point once all the data is entered
(new_data <- d_spawn %>%
  filter(year == last_yr-1) %>%
  select(spawnerstm4))
d_spawn <- d_spawn[complete.cases(d_spawn),] 

#more processing data long
data_l <- data %>%
  select(starts_with("oage")) %>%
  gather(key = oage, value = fish, c(oage_1:oage_4)) 

#check missing values There should be some for the older ages since those fish haven't returned yet
data_l[!complete.cases(data_l),]

#IF only the most recent years for age classes are missing remove them.   ... other wise figure out why they are missing!
data_l <- na.omit(data_l)  

data_l <- data_l %>%
  separate(oage, c("do_delete", "oage")) %>%
  select(-do_delete) 

data_l$oage <- as_factor(data_l$oage)


# median return by size ----

# For ages 1, 3 and 4, use the median from the last 10 years
(quant <- data_l %>%
   group_by(oage) %>%
   summarize(lwr90 = quantile(tail(fish, 10), c(0.10)), 
             est = median(tail(fish, 10)),
             upr90 = quantile(tail(fish, 10), c(.90))))

# time series analysis ----

data_ts <- data %>%
  select(run)

#check missing values
data_ts[!complete.cases(data_ts),]

#remove missing values for time sereies analysis
data_ts <- data_ts[complete.cases(data_ts),] %>%
  as.ts()

# into a time series data set
# analysis 
#autoplot(data_ts)
#ggAcf(data_ts)
#ggAcf(diff(data_ts))

# Ran a variety of time series models, checking diagnostics along the way.

#naive (last years value)
fc <- fcn <- naive(data_ts, h = 1)
#simple exponential smoothing
fc <- fcses <- ses(data_ts, h = 1)
#holt
fc <- fch <- holt(y = data_ts, h = 1,  exponential = FALSE) 

#holt for the 2020 forecast don't use this method since residuals don't pass the Ljung-Box test.
fc <- fch <- holt(y = data_ts, h = 1,  exponential = TRUE) # Getting an error
#damped holt
#fc <- holt(y = data_ts, h = 1, level = c(80, 80),  damped = TRUE, lambda = "auto") #can specify different Prediction intervals if needed. 
fc <- holt(y = data_ts, h = 1, damped = TRUE, lambda = "auto") 


summary(fc)
dev.off()
checkresiduals(fc)
#accuracy(fc)
autoplot(fc) + autolayer(fitted(fc))

#Checking we have the right number of things also used inspection to check that
#the most recent values are accurate.
length(fitted(fc))

#For the 2020 forecast the model with the smallest RMSE and MAPE is either the
#SES with alpha = .9999 which is virtually the Naive estimate: Just use last
#years value. The following has an even smaller RMSE, especially cross
#validated.

# analysis ----

my_exp_summary <- function (data, lev = NULL, model = NULL) {
  c(RMSE = sqrt(mean((expm1(data$obs) - expm1(data$pred)^2))),
    Rsquared = summary(lm(pred ~ obs, data))$r.squared,
    MAE = mae(expm1(data$obs), expm1(data$pred)),
    MAPEEXP = mape(expm1(data$obs), expm1(data$pred)))
}
my_summary <- function (data, lev = NULL, model = NULL) {
  c(RMSE = sqrt(mean((data$obs -data$pred)^2)),
    Rsquared = summary(lm(pred ~ obs, data))$r.squared,
    MAE = mae((data$obs), (data$pred)),
    MAPE = mape((data$obs), (data$pred)))
}

#model 1 ----
# ocean age 2 fish are regressed on spawners for the ocean age 2 prodgeny.
lm2s <- lm(oage_2 ~ spawnerstm4 , data = d_spawn)
summary(lm2s)
#to annotate the graph need library(grid)
rsq <- summary(lm2s)$adj.r.squared
pvalue <- summary(lm2s)$coefficients[2,4]
rp <- paste0("adj.r^2 = ", round(rsq,2), "  pvalue = ", round(pvalue, 3))

dev.off()
layout(matrix(c(1,2,3,4),2,2)) # optional 4 graphs/page
plot(lm2s) # check for normality

newpoint <- broom::augment(lm2s, newdata = new_data)
(pred <- predict(lm2s, newdata = new_data, interval = "prediction", level = 0.80))
lwr <- pred[2] # needed for ggplot
upr <- pred[3]
predict(lm2s, newdata = new_data, interval = "prediction", level = 0.80)

#Use to make 95% CI and PI 
minspawnerstm4 <- min(data$spawnerstm4, na.rm = TRUE)
maxspawnerstm4 <- max(data$spawnerstm4, na.rm = TRUE)
predx <- data.frame(spawnerstm4 = seq(from = minspawnerstm4, to = maxspawnerstm4, by = (maxspawnerstm4-minspawnerstm4)/19))

# ... confidence interval
conf.int <- cbind(predx, predict(lm2s, newdata = predx, interval = "confidence", level = 0.80))

# ... prediction interval
pred.int <- cbind(predx, predict(lm2s, newdata = predx, interval = "prediction", level = 0.80))

g.pred <- ggplot(pred.int, aes(x = spawnerstm4, y = fit)) +
  geom_smooth(data = pred.int, aes(ymin = lwr, ymax = upr), stat = "identity") + # prediction interval
  geom_smooth(data = conf.int, aes(ymin = lwr, ymax = upr), stat = "identity") + #confidence interval
  geom_point(data = data, aes(x = spawnerstm4, y = oage_2)) + #plots all the points
  geom_text_repel(data = data, aes(x = spawnerstm4, y = oage_2, label = year)) +
  geom_point(data = newpoint, aes(y = .fitted), size = 3, color = "red") + # adds this years new point
  geom_text_repel(data = newpoint, aes(x = spawnerstm4, y = .fitted, label = round(.fitted, 0 )), adj = 1) +  
  #annotate("text", label = rp, x = 205000, y = 550000) + 
  stat_regline_equation(label.x = 150000, label.y = 600000) +
  theme_bw() +
  #coord_cartesian(ylim = c(0, 500000), xlim = c(0, 850000)) +
  xlab("Spawners time minus 4 years") +
  ylab("ocean age 2") #+ #ggtitle("oage_2 vs spawnerstm4")
g.pred  
dev.off()
ggsave(filename = paste0("figures/oage_2_spawnerstm4_2003up", ".png", sep = ""), device = png(), width = 7, height = 9, units = "in", dpi = 300)

#Repeated K- fold Cross validation ----

# there should be no NAs
d_spawn[!complete.cases(d_spawn),]

#IF only the most recent years for age classes are missing remove them.   ...
#other wise figure out why they are missing!
data_cv <- data %>%
  select(spawnerstm4, oage_2) %>%
  na.omit()  #can't have NA's for cross validation.

data_cv[!complete.cases(data_cv),]
#data <- data_cv
# define training control use one fo the following lines
train_control <- trainControl(method = "repeatedcv", number = 2, repeats = 8, summaryFunction = my_summary)

# train the model Warning messages are okay.
model <- train(oage_2 ~ spawnerstm4, data = data_cv, trControl=train_control, method="lm")
# summarize result
print(model)

#forecast summary ----
pred
quant
quant[2,3] <- pred[1]
quant[2,2] <- pred[2]
quant[2,4] <- pred[3]

#check it matches worksheet.
(lwr <- sum(quant$lwr90[1], quant$lwr90[3:4]))
(est <- sum(quant$est[1], quant$est[3:4]))
(upr <- sum(quant$upr90[1],  quant$upr90[3:4]))

lwr <- sum(quant$lwr90[1], pred[2], quant$lwr90[3:4])
est <- sum(quant$est[1], pred[1], quant$est[3:4])
upr <- sum(quant$upr90[1], pred[3], quant$upr90[3:4])

quant %>%
  summarize(lwr = sum(lwr),
            est = sum(est),
            upr = sum(upr)) 

(nelson_f <- data.frame(est, lwr, upr))
nelson_f$est

#additional for report ----
escapement_goal <- 158000
(harvest_est <- nelson_f$est - escapement_goal )

