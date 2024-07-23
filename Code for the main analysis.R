rm(list=ls())
library(tidyverse)
library(foreign)
library(haven)
library(pROC)
library(splines)
library(growthrates)
library(parallel)
library(doParallel)
library(SuperLearner)
library(gbm)
library(rms)
library(tableone)
library(missRanger)
library(shapviz)
library(data.table)
library(cutpointr)
library(broom)

NTproBNP = read.xport("E:/NHANES/SSBNP_A.XPT") %>% filter(!is.na('SSBNP'))
cbc1999 = read.xport("E:/NHANES/LAB25.XPT")
cbc2001 = read.xport("E:/NHANES/L25_B.XPT")
cbc2003 = read.xport("E:/NHANES/L25_C.XPT")
SoB1999 = read.xport("E:/NHANES/CDQ.XPT") %>% dplyr::select(SEQN, CDQ010)
SoB2001 = read.xport("E:/NHANES/CDQ_B.XPT") %>% dplyr::select(SEQN, CDQ010)
SoB2003 = read.xport("E:/NHANES/CDQ_C.XPT") %>% dplyr::select(SEQN, CDQ010)
hb = rbind(cbc1999, cbc2001, cbc2003) %>% dplyr::select(SEQN, LBXHGB)
SoB.on.stair = rbind(SoB1999, SoB2001, SoB2003) %>% 
  mutate(CDQ010 = ifelse(is.na(CDQ010)|CDQ010 >= 7, NA, CDQ010))

NHANES = read_dta("E:/NHANES/NHANES 1999-2018.dta")

D4= NHANES %>% 
  inner_join(NTproBNP %>% dplyr::select(SEQN, SSBNP), by = c("seqn" = "SEQN")) %>% 
  left_join(hb, by = c("seqn" = "SEQN")) %>% 
  left_join(SoB.on.stair, by = c("seqn" = "SEQN")) %>% 
  filter(ridageyr >= 20 & ridageyr <= 79)

ggplot(D4, aes(x = ridageyr, y = SSBNP)) + geom_smooth() + geom_jitter() + ylim(0,500)

attributes_nh = D4 %>% 
  map_df(~tibble(var_name = names(.),
                 label = attr(., "label")),
         .id = "var_name") %>%
  select(var_name, label)

#exclude those who have prior_cvd
D4$prior_cvd <- with(D4, ifelse(mcq160b==1 | mcq160c==1 | mcq160d==1 | 
                                  mcq160e==1 | mcq160f==1, 1, 0))
D4 <- D4 %>% filter(prior_cvd == 0)
 
#smoking_status 0:Never Smoked at least 100 cigarettes in life,  1:Former, 2: Current,
D4$Smoking <- with(D4, ifelse(smq020==2,0,(ifelse(smq040==1|smq040==2,2,
                                                         1))))

#Marital_status
D4$Marital_status <- with(D4, ifelse(dmdmartl==1, 1, 0))

#education_status
D4$Education_status <- with(D4, case_when(dmdeduc2 == 1 ~ 1,
                                          dmdeduc2 == 2 ~ 2,
                                          dmdeduc2 == 3 ~ 3,
                                          dmdeduc2 %in% c(4,5)  ~ 4,
                                          TRUE ~ NA))

#Hypertension
D4$HTN <- with(D4, case_when(bpq020 == 1 ~ 1,
                             bpq020 == 2 ~ 0,
                             TRUE ~ NA))
#DM Drugs
D4 <- D4 %>% mutate(DM = case_when(diq010 == 1 ~ 1,
                                         diq010 %in% c(2,3) ~ 0,
                                         TRUE ~ NA))

#ethnicity_group
D4$Ethnicity <- with(D4, ifelse(ridreth1==3, 1, ifelse(ridreth1==4, 2, ifelse(ridreth1==1, 3, 0))))

#Family income
D4 <- D4 %>% mutate(Family_income= case_when(indfminc %in% c(1,2,3,4,13) ~ 1, #$0-19,999
                                                       indfminc %in% c(5,6,7,8) ~ 2, #$20,000-54,999
                                                       indfminc %in% c(9,10,11) ~ 3,  #$55,000-
                                                       TRUE ~ NA))
#SBP
D4 <- D4 %>% mutate(SBP = ifelse(bpxsar != 0, bpxsar, NA))

#HDL
D4 <- D4 %>% mutate(HDL = ifelse(sddsrvyr == 3, lbxhdd, lbdhdl))

#eGFR
#gender: 1 male, 2, female
D4 <- D4 %>% 
  mutate(cre = ifelse(sddsrvyr == 2, lbdscr, lbxscr), 
         age = ridageyr, gender = riagendr,
         CKDeq_A = ifelse(gender == 1, 0.9, 0.7),
         CKDeq_B = ifelse(cre <= CKDeq_A,ifelse(gender == 1, -0.302, -0.241) ,-1.2),
         eGFR = 142*((cre/CKDeq_A)^CKDeq_B)*(0.9938^age)*(1+0.012*(gender-1)))


#Rename with keeping original variable names
D4 <- D4 %>% mutate(Female= riagendr - 1,
             Age = ridageyr,
             BMI = bmxbmi,
             Alb = lbxsal,
             HbA1c = lbxgh,
             AST = lbxsassi,
             ALT = lbxsatsi,
             BUN = lbxsbu,
             Hb = LBXHGB,
             TC = lbxtc,
             TG = lbxstr,
             Na = lbxsnasi,
             K = lbxsksi,
             Statin = ifelse(is.na(statin), 0,1)) 

#Convert to categorical variables
D4 <- D4 %>% mutate_at(.vars = c("sddsrvyr","Female","Smoking", "HTN", "DM", "Marital_status",
                                 "Ethnicity", "Statin", "Family_income", "Education_status"), .funs = as.factor)

####Prediction model
names<-c("sddsrvyr", "Female","Ethnicity", "Marital_status", "Education_status",
          "Family_income","Smoking", 
          "HTN","DM", "Statin") 
ses <- c("Marital_status", "Education_status", "Family_income")
 
cov<-c("Age", "BMI","SBP","eGFR", "HbA1c", "AST", "ALT", "BUN", "Hb",
       "Alb", "TC", "HDL", "Na", "K") 

D4$outcome <- with(D4, ifelse(SSBNP>125.00,1,0))
select<-c(names, cov, "outcome", "CDQ010")
dselect <- D4[select]

cores <- detectCores(logical = FALSE)
registerDoParallel(cores = cores)
set.seed(1234)
impDF <- missRanger(dselect, . -outcome - CDQ010
                       ~ . -outcome - CDQ010,
                       maxiter = 10,
                       num.trees = 100, verbose = 0)
newd <- impDF
summary(newd)

# Analysis function -------------------------------------------------------
drop.SES = FALSE
use.sob = FALSE
thsld.obs.plv = FALSE
minimum = FALSE

analysis <- function(newd, drop.SES = FALSE, use.sob = FALSE, thsld.obs.plv = FALSE, minimum = FALSE) {
  if (use.sob == TRUE){
    newd <- newd %>% filter(CDQ010 == 2)
  }
  newd <- newd %>% select (-CDQ010)

  if (drop.SES == TRUE){
    newd <- newd %>% dplyr::select(-all_of(ses))
  }
  if(minimum == TRUE){
    newd <- newd %>% dplyr::select(sddsrvyr, Age, Female, BMI, SBP, HTN, DM, Statin, outcome)
  }
  
filename.suffix <- c("")
filename.suffix <- ifelse(drop.SES == TRUE, paste(filename.suffix, ".withoutSES", sep=""), filename.suffix)
filename.suffix <- ifelse(use.sob == TRUE, paste(filename.suffix, ".noSOB", sep=""), filename.suffix)
filename.suffix <- ifelse(thsld.obs.plv == TRUE, paste(filename.suffix, ".obs.plv", sep=""), filename.suffix)
filename.suffix <- ifelse(minimum == TRUE, paste(filename.suffix, ".minimum", sep=""), filename.suffix)
  
d_train <- newd[newd$sddsrvyr %in% c(1,2), ]  
d_test  <- newd[newd$sddsrvyr==3, ] 
d_train<-select(d_train, subset=-sddsrvyr)
d_test<-select(d_test, subset=-sddsrvyr)

#################################################################
## Reference model -- logistic regression with triage variable ##
#################################################################

# Fit logistic regression
fit_logistic<- glm(outcome ~ .-outcome,
                   family=binomial(link="logit"), 
                   data=d_train) 

# Prediction in test set
pred_logistic_prob <- predict(fit_logistic, 
                              newdata = d_test, type = "response")

# ROC 
(roc_logistic <- pROC::roc(response = d_test$outcome,
                     predictor= pred_logistic_prob) )
ci.auc(d_test$outcome, pred_logistic_prob) 
plot(roc_logistic, legacy.axes=TRUE)


# Fit logistic regression WITH SPLINE
fit_spline<- glm(outcome ~ . - outcome + rcs(SBP, knots=3) + rcs(BMI, knots = 3),
                 family=binomial(link="logit"), 
                 data=d_train) 

# Prediction in test set
pred_spline_prob <- predict(fit_spline, 
                            newdata = d_test, type = "response")

# ROC 
(roc_spline <- pROC::roc(response = d_test$outcome,
                   predictor= pred_spline_prob) )
ci.auc(d_test$outcome, pred_spline_prob) 

###################################################
## Logistic regression with Lasso regularization ##
###################################################
# Split into training/tnamesest sets

library("glmnet")
d_train2<-d_train
d_test2<-d_test
# Create matrices of training set
y <- d_train2 %>% select(outcome) %>% 
  data.matrix() 
x <- d_train2 %>% select(-outcome) %>%
  data.matrix()

# Prediction in test set
# Create matrices of test set
y_test <- d_test2 %>% select(outcome)  %>% data.matrix() 
x_test <- d_test2 %>% select(-outcome) %>% data.matrix()

# Sanity check -- no NA alllowed
sum(is.na(y)); sum(is.na(x))

# Fit Lasso in training set with cross validation to identify the best lambda 
set.seed(1)
fit_lasso_cv <- cv.glmnet(x, y, 
                          family = "binomial",
                          type.measure = "mse", nfolds =10,
                          standardize = TRUE) 

plot(fit_lasso_cv)
fit_lasso_cv$lambda.min # minimal lambda
coef(fit_lasso_cv, s = "lambda.min") # coefficients for min lambda

#Prediction
pred_lasso_prob <- predict(fit_lasso_cv, newx = x_test, 
                           s = "lambda.min", type="response") %>% as.vector()
hist(pred_lasso_prob)
# ROC
(roc_lasso <- pROC::roc(response  = y_test, 
                  predictor = pred_lasso_prob) )
ci.auc(y_test, pred_lasso_prob)

###################
## Random forest ##
###################
# Split into training/test sets
library("caret")
set.seed(1, "L'Ecuyer-CMRG")
cl<-makeCluster(16)
registerDoParallel(cl)
d_train$outcome<-factor(d_train$outcome) %>% factor(labels = c("No", "Yes"))
# Set Training Control
myTrainingControl <- trainControl(method = "cv", 
                                  number = 10, 
                                  savePredictions = TRUE, 
                                  classProbs = TRUE, 
                                  verboseIter = FALSE)
# Train RF
fit_RF <- caret::train(outcome ~ .,   
                       data = d_train, 
                       method = "ranger", 
                       #tuneLength = 3,     
                       importance = "permutation",
                       trControl = myTrainingControl)
stopCluster(cl)
print(fit_RF)


# Prediction in test set
pred_RF_prob <- predict(fit_RF, d_test, type="prob")
hist(pred_RF_prob$Yes)

# ROC
(roc_RF <- pROC::roc(response = d_test$outcome,
               predictor= pred_RF_prob$Yes) )
ci.auc(d_test$outcome, pred_RF_prob$Yes)

###########################################################
## Gradient boosted decision tree (GBDT) -- xgboost:tree ##
###########################################################
# Split into training/test sets
library("xgboost")
set.seed(1, "L'Ecuyer-CMRG")
cl<-makeCluster(16)
registerDoParallel(cl)
# Set Training Control
d_train$outcome<-factor(d_train$outcome) %>% factor(labels = c("No", "Yes"))
myTrainingControl <- trainControl(method = "cv", 
                                  number = 10, 
                                  savePredictions = TRUE, 
                                  classProbs = TRUE, 
                                  verboseIter = FALSE) 

fit_xgbTree <- caret::train(outcome ~ ., 
                            data = d_train, 
                            method = "xgbTree", 
                            #tuneLength = 3,
                            trControl = myTrainingControl)
stopCluster(cl)
print(fit_xgbTree)

# Prediction in test set
pred_xgb_prob <- predict(fit_xgbTree, d_test, type="prob")

# ROC
(roc_xgb <- pROC::roc(response = d_test$outcome,
                predictor= pred_xgb_prob$Yes) )
ci.auc(d_test$outcome, pred_xgb_prob$Yes)


####################
## Super learner  ##
####################
SL.library<- c("SL.glm", "SL.glmnet",  "SL.gbm", "SL.randomForest")
folds=10 

library("glmnet")
outcome<-as.numeric(d_train$outcome)-1

set.seed(1, "L'Ecuyer-CMRG")
cl<-makeCluster(16)
registerDoParallel(cl)

fitY<-SuperLearner(Y=outcome, X=subset(d_train, select=-outcome), family="binomial",
                   method="method.AUC",
                   SL.library=SL.library,
                   cvControl=list(V=folds))
stopCluster(cl)
fitY


# Obtain the predicted probability of the outcome from SL
pred_sl_prob<-predict(fitY, newdata=d_test %>% dplyr::select(-outcome))$pred


# ROC
(roc_sl <- pROC::roc(response = d_test$outcome,
               predictor= pred_sl_prob) )
ci.auc(d_test$outcome, pred_sl_prob)
plot(roc_sl, legacy.axes=TRUE)

##################################################################################
svg(filename = paste("E:/NHANES/ROC_", filename.suffix, ".svg", sep = ""))
plot(roc_logistic, legacy.axes=TRUE, col="black", lty="longdash", xlim = c(1, 0), ylim= c(0, 1))
plot(roc_spline, legacy.axes=TRUE, col="gray", lty="longdash", add=TRUE)
plot(roc_lasso, legacy.axes=TRUE, col="blue", lty="longdash", add=TRUE) 
plot(roc_RF, legacy.axes=TRUE, col="purple", lty="longdash", add=TRUE)
plot(roc_xgb, legacy.axes=TRUE, col='green', lty="longdash", add=TRUE)
plot(roc_sl, legacy.axes=TRUE, col="red", lty="longdash", add=TRUE)
legend(0.4, 0.5, legend=c("Logistic", "Logistic+Spline", "Lasso", "RF", "GBM", "SuperLearner"),
       col=c("black", "gray", "blue", "purple", "green", "red"), lty=2:2, cex=0.7)

dev.off()

ggsave(filename = paste("E:/NHANES/VarImp_", filename.suffix, ".svg", sep = ""), plot = VarImp, dpi = 320, device = "svg")

if (thsld.obs.plv == FALSE){
  
# Confusion Matrix 
prob_log<- predict(fit_logistic, newdata = d_test, type = "response")
out<-(as.numeric(d_test$outcome)-1)

dtemp<-as.data.frame(cbind(prob_log, out))
cp <- cutpointr(dtemp, prob_log,  out, method = maximize_metric, metric = youden)
prevalence<-cp$optimal_cutpoint

d_test$outcome<-factor(d_test$outcome) %>% factor(labels = c("No", "Yes"))
pred_logistic_class <- ifelse(pred_logistic_prob >=  prevalence, 1, 0) %>% factor(labels = c("No", "Yes"))
CFM_logistic <- caret::confusionMatrix(data      = pred_logistic_class,   
                       reference = d_test$outcome,
                       mode      = "sens_spec",
                       positive  = "Yes")

# Confusion Matrix 

prob_spline<- predict(fit_spline, newdata = d_test, type = "response")
out<-(as.numeric(d_test$outcome)-1)

dtemp<-as.data.frame(cbind(prob_spline, out))
cp <- cutpointr(dtemp, prob_spline,  out, method = maximize_metric, metric = youden)
prevalence<-cp$optimal_cutpoint

d_test$outcome<-factor(d_test$outcome) %>% factor(labels = c("No", "Yes"))
pred_spline_class <- ifelse(pred_spline_prob >=  prevalence, 1, 0) %>% factor(labels = c("No", "Yes"))
CFM_spline <- caret::confusionMatrix(data      = pred_spline_class,   
                       reference = d_test$outcome,
                       mode      = "sens_spec",
                       positive  = "Yes")

# Confusion matrix
prob_lasso<- predict(fit_lasso_cv, newx = x_test, 
                     s = "lambda.min", type="response") %>% as.vector()
out<-(as.numeric(d_test$outcome)-1)

dtemp<-as.data.frame(cbind(prob_lasso, out))
cp <- cutpointr(dtemp, prob_lasso, out, method = maximize_metric, metric = youden)
prevalence<-cp$optimal_cutpoint

d_test$outcome<-factor(d_test$outcome) %>% factor(labels = c("No", "Yes"))
pred_lasso_class <- ifelse(pred_lasso_prob >=  prevalence, 1, 0) %>% factor(labels = c("No", "Yes"))
CFM_lasso <- caret::confusionMatrix(data      = pred_lasso_class, 
                       reference = d_test$outcome,
                       mode      = "sens_spec",
                       positive  = "Yes")


# Confusion Matrix
prob_rf<- predict(fit_RF, d_test, type="prob") %>% 
  dplyr::select(Yes) %>% 
  as_vector()
out<-(as.numeric(d_test$outcome)-1)

dtemp<-as.data.frame(cbind(prob_rf, out))
cp <- cutpointr(dtemp, prob_rf, out, method = maximize_metric, metric = youden)
prevalence<-cp$optimal_cutpoint
d_test$outcome<-factor(d_test$outcome) %>% factor(labels = c("No", "Yes"))
pred_RF_class <- ifelse(pred_RF_prob$Yes >= prevalence, 1, 0) %>% factor(labels = c("No", "Yes"))
CFM_RF <- confusionMatrix(data      = pred_RF_class,    
                reference = d_test$outcome,
                mode      = "sens_spec",
                positive  = "Yes")


# Confusion Matrix
prob_gbm<- predict(fit_xgbTree, d_test, type="prob")%>% 
  dplyr::select(Yes) %>% 
  as_vector()
out<-(as.numeric(d_test$outcome)-1)

dtemp<-as.data.frame(cbind(prob_gbm, out))
cp <- cutpointr(dtemp, prob_gbm,  out, method = maximize_metric, metric = youden)
prevalence<-cp$optimal_cutpoint

d_test$outcome<-factor(d_test$outcome) %>% factor(labels = c("No", "Yes"))
pred_xgb_class <- ifelse(pred_xgb_prob$Yes >= prevalence, 1, 0) %>% factor(labels = c("No", "Yes"))
CFM_xgb <- confusionMatrix(data      = pred_xgb_class,    
                reference = d_test$outcome,
                mode      = "sens_spec",
                positive  = "Yes")

# Confusion Matrix

prob_sl<- predict(fitY, newdata=d_test %>% dplyr::select(-outcome))$pred
out<-(as.numeric(d_test$outcome)-1)

dtemp<-as.data.frame(cbind(prob_sl, out))
cp <- cutpointr(dtemp, prob_sl,  out, method = maximize_metric, metric = youden)
#cp <- cutpointr(dtemp, prob_sl,  out, method = oc_youden_kernel)
prevalence<-cp$optimal_cutpoint
d_test$outcome<-factor(d_test$outcome) %>% factor(labels = c("No", "Yes"))
pred_sl_class <- ifelse(pred_sl_prob >= prevalence, 1, 0) %>% factor(labels = c("No", "Yes"))
CFM_sl <- confusionMatrix(data      = pred_sl_class,    
                reference = d_test$outcome,
                mode      = "sens_spec",
                positive  = "Yes")

CFM = rbind(CFM_logistic$byClass,
      CFM_spline$byClass,
      CFM_lasso$byClass,
      CFM_RF$byClass,
      CFM_xgb$byClass,
      CFM_sl$byClass)[,1:4] 
rownames(CFM) <- c("Logistic", "Spline", "Lasso", "RF", "Gradient boosting", "SuperLearner")
ROCs <- rbind(ci.auc(d_test$outcome, pred_logistic_prob)[c(2,1,3)],
              ci.auc(d_test$outcome, pred_spline_prob)[c(2,1,3)],
              ci.auc(d_test$outcome, pred_lasso_prob)[c(2,1,3)],
              ci.auc(d_test$outcome, pred_RF_prob$Yes)[c(2,1,3)],
              ci.auc(d_test$outcome, pred_xgb_prob$Yes)[c(2,1,3)],
              ci.auc(d_test$outcome, pred_sl_prob)[c(2,1,3)])
rownames(ROCs) <- c("Logistic", "Spline", "Lasso", "RF", "Gradient boosting", "SuperLearner")
colnames(ROCs) <- c("AUC", "95%CI_L", "95%CI_H")

result.table = cbind(ROCs, CFM)
result.table

#DeLong's test
roc.test(roc_logistic, roc_spline, method = "delong")
roc.test(roc_logistic, roc_lasso, method = "delong")
roc.test(roc_logistic, roc_RF, method = "delong")
roc.test(roc_logistic, roc_xgb, method = "delong")
roc.test(roc_logistic, roc_sl, method = "delong")


write.csv(result.table, paste("E:/NHANES/results", filename.suffix, ".csv", sep = ""))
} else if (thsld.obs.plv == TRUE) {
  #Confusion matrix with observed-prevalence cutoff
  d_test$outcome<-factor(d_test$outcome) %>% factor(labels = c("No", "Yes"))
  prevalence = mean(newd$outcome)
  prevalence
  prob_log<- predict(fit_logistic, newdata = d_test, type = "response")
  pred_logistic_class <- ifelse(pred_logistic_prob >=  prevalence, 1, 0) %>% factor(labels = c("No", "Yes"))
  CFM_logistic_obs.prevalence <- caret::confusionMatrix(data      = pred_logistic_class,   
                                                        reference = d_test$outcome,
                                                        mode      = "sens_spec",
                                                        positive  = "Yes")
  
  prob_spline<- predict(fit_spline, newdata = d_test, type = "response")
  pred_spline_class <- ifelse(pred_spline_prob >=  prevalence, 1, 0) %>% factor(labels = c("No", "Yes"))
  CFM_spline_obs.prevalence <- caret::confusionMatrix(data      = pred_spline_class,   
                                                      reference = d_test$outcome,
                                                      mode      = "sens_spec",
                                                      positive  = "Yes")
  
  prob_lasso<- predict(fit_lasso_cv, newx = x_test, 
                       s = "lambda.min", type="response") %>% as.vector()
  pred_lasso_class <- ifelse(pred_lasso_prob >=  prevalence, 1, 0) %>% factor(labels = c("No", "Yes"))
  CFM_lasso_obs.prevalence <- caret::confusionMatrix(data      = pred_lasso_class, 
                                                     reference = d_test$outcome,
                                                     mode      = "sens_spec",
                                                     positive  = "Yes")
  
  prob_rf<- predict(fit_RF, d_test, type="prob") %>% 
    dplyr::select(Yes) %>% 
    as_vector()
  pred_RF_class <- ifelse(pred_RF_prob$Yes >= prevalence, 1, 0) %>% factor(labels = c("No", "Yes"))
  CFM_RF_obs.prevalence <- confusionMatrix(data      = pred_RF_class,    
                                           reference = d_test$outcome,
                                           mode      = "sens_spec",
                                           positive  = "Yes")
  
  prob_gbm<- predict(fit_xgbTree, d_test, type="prob")%>% 
    dplyr::select(Yes) %>% 
    as_vector()
  pred_xgb_class <- ifelse(pred_xgb_prob$Yes >= prevalence, 1, 0) %>% factor(labels = c("No", "Yes"))
  CFM_xgb_obs.prevalence <- confusionMatrix(data      = pred_xgb_class,    
                                            reference = d_test$outcome,
                                            mode      = "sens_spec",
                                            positive  = "Yes")
  
  prob_sl<- predict(fitY, newdata=d_test %>% dplyr::select(-outcome))$pred
  pred_sl_class <- ifelse(pred_sl_prob >= prevalence, 1, 0) %>% factor(labels = c("No", "Yes"))
  CFM_sl_obs.prevalence <- confusionMatrix(data      = pred_sl_class,    
                                           reference = d_test$outcome,
                                           mode      = "sens_spec",
                                           positive  = "Yes")
  
  CFM = rbind(CFM_logistic_obs.prevalence$byClass,
              CFM_spline_obs.prevalence$byClass,
              CFM_lasso_obs.prevalence$byClass,
              CFM_RF_obs.prevalence$byClass,
              CFM_xgb_obs.prevalence$byClass,
              CFM_sl_obs.prevalence$byClass)[,1:4] 
  rownames(CFM) <- c("Logistic", "Spline", "Lasso", "RF", "Gradient boosting", "SuperLearner")
  write.csv(CFM, paste("E:/NHANES/results", filename.suffix, ".csv", sep = ""))
  
}

}

#main analysis
analysis(impDF, drop.SES = FALSE, use.sob = FALSE, thsld.obs.plv = FALSE)

#use observed prevalence
analysis(impDF, drop.SES = FALSE, use.sob = FALSE, thsld.obs.plv = TRUE)

#drop SES
analysis(impDF, drop.SES = TRUE, use.sob = FALSE, thsld.obs.plv = FALSE)

#free from SoB
analysis(impDF, drop.SES = FALSE, use.sob = TRUE, thsld.obs.plv = FALSE)

#minimum
analysis(impDF, minimum = TRUE)