library(readxl)
Mayo_Clinic_PBC_trial_data_restricted <- read_excel("C:/Users/Acer/Desktop/Mayo Clinic PBC trial data restricted.xlsx")
delta = as.numeric(Mayo_Clinic_PBC_trial_data_restricted$status==2)
Mayo_Clinic_PBC_trial_data_restricted$delta = delta

findpct = function(a){
  overall = mean(a) #overall
  drug = mean(a[which(treatment==1)]) #drug
  control = mean(a[which(treatment==2)]) #control
  val = c(overall,drug,control)
  names(val) = c("overall","drug","control")
  return(val)
}

findpct(ascites)
findpct(hepatomegaly)
findpct(spiders)
findpct(edema)
hist(ascites)
hist(hepatomegaly)
hist(spiders)
hist(edema)

library(km.ci)
library(KMsurv)
library(survival)
library(My.stepwise)
attach(Mayo_Clinic_PBC_trial_data_restricted)
#The following gives us which within the data are not censored. The rest are considered indpendently censored at the time given.
summary(Mayo_Clinic_PBC_trial_data_restricted)
surv.trt1 = Surv(time[treatment==1],delta[treatment==1])
surv.trt2 = Surv(time[treatment==2],delta[treatment==2])

fit1 = survfit(surv.trt1~1)
fit2 = survfit(surv.trt2~1)
band.1 = km.ci(fit1,method="logep")
plot(band.1,main="Survival Function: Treatment")
band.2 = km.ci(fit2,method="logep")
plot(band.2,main="Survival Function: Placebo")
sv = Surv(time,delta)

#Before fitting the Cox model we have to deal with some missing data. We might try fitting cholesterol and platelets for the missing data
#based on some of the other values such as age, sex, and bilirubin. For this missing data we will simply 
my.variable.list = c("sex","age","ascites","hepatomegaly","spiders","edema","bilirubin","albumin","alkaline","prothrombin")

trunc.cholest = Mayo_Clinic_PBC_trial_data_restricted[which(Mayo_Clinic_PBC_trial_data_restricted$cholesterol!="."),]
trunc.cholest$cholesterol = as.numeric(trunc.cholest$cholesterol)
trunc.platelets = Mayo_Clinic_PBC_trial_data_restricted[which(Mayo_Clinic_PBC_trial_data_restricted$platelets!="."),]
trunc.platelets$platelets = as.numeric(trunc.platelets$platelets)

#My.stepwise.lm(Y="cholesterol",variable.list = my.variable.list,data=trunc.cholest)
#Gives us a model that includes bilirubin, intercept, edema, age, alkaline, and prothrombin
model.cholest = lm(cholesterol~1+bilirubin+edema+age+alkaline+prothrombin,data=trunc.cholest)

Mayo_Clinic_PBC_trial_data_restricted$cholesterol[which(Mayo_Clinic_PBC_trial_data_restricted$cholesterol==".")] = predict(model.cholest,Mayo_Clinic_PBC_trial_data_restricted[which(Mayo_Clinic_PBC_trial_data_restricted$cholesterol=="."),])

Mayo_Clinic_PBC_trial_data_restricted$cholesterol = as.numeric(Mayo_Clinic_PBC_trial_data_restricted$cholesterol)


My.stepwise.lm(Y="platelets",my.variable.list,data=trunc.platelets)

model.platelets = lm(platelets~1+edema+hepatomegaly+alkaline+prothrombin+albumin+sex,trunc.platelets)

Mayo_Clinic_PBC_trial_data_restricted$platelets[which(Mayo_Clinic_PBC_trial_data_restricted$platelets==".")] = predict(model.platelets,Mayo_Clinic_PBC_trial_data_restricted[which(Mayo_Clinic_PBC_trial_data_restricted$platelets=="."),])

Mayo_Clinic_PBC_trial_data_restricted$platelets = as.numeric(Mayo_Clinic_PBC_trial_data_restricted$platelets)

Mayo_Clinic_PBC_trial_data_restricted$older= as.numeric(age > median(age))
Mayo_Clinic_PBC_trial_data_restricted$high.bili = as.numeric(bilirubin > 1.9)
Mayo_Clinic_PBC_trial_data_restricted$low.alb = as.numeric(albumin < 3.4)
Mayo_Clinic_PBC_trial_data_restricted$plate.abnormal = as.numeric(platelets < 150000 | platelets > 450000)
Mayo_Clinic_PBC_trial_data_restricted$prothrombin.over = as.numeric(prothrombin>14)
Mayo_Clinic_PBC_trial_data_restricted$prothrombin.under = as.numeric(prothrombin < 10)

test1 = survdiff(sv~treatment, rho = 0)

#Tests for difference weighted to earlier and later
test2 = survdiff(sv~treatment, rho = 1)
test3 = survdiff(sv~treatment, rho = -1)
test4 = survdiff(sv~treatment+strata(older))
testsex = survdiff(sv~treatment+strata(sex))
test5 = survdiff(sv~treatment+strata(high.bili))
test6 = survdiff(sv~treatment+strata(low.alb))
test7 = survdiff(sv~treatment+strata(plate.abnormal))
test8 = survdiff(sv~treatment+strata(prothrombin.over))
test9 = survdiff(sv~treatment+strata(ascites))
test10 = survdiff(sv~treatment+strata(hepatomegaly))
test11 = survdiff(sv~treatment+strata(spiders))
test12  = survdiff(sv~treatment+strata(edema))


test10 

my.variable.list = c("low.alb","plate.abnormal","prothrombin.over","high.bili","treatment","sex","age","ascites","hepatomegaly","spiders","edema","platelets","cholesterol","bilirubin","albumin","alkaline","prothrombin")
My.stepwise.coxph(Time = "time", Status = "delta", variable.list = my.variable.list,
                 data = Mayo_Clinic_PBC_trial_data_restricted)

My.stepwise.coxph(Time = "time", Status = "delta", variable.list = my.variable.list, in.variable = c("treatment"),
                  data = Mayo_Clinic_PBC_trial_data_restricted)


#The above was run to get the final model:

library(MASS)
cox.model = coxph(formula = Surv(time, delta) ~ + high.bili + albumin + prothrombin + bilirubin + age + edema, data = Mayo_Clinic_PBC_trial_data_restricted, method = "efron")

cox.model.trt = coxph(formula = Surv(time, delta) ~ treatment + high.bili + albumin + prothrombin + bilirubin + age + edema, data = Mayo_Clinic_PBC_trial_data_restricted, method = "efron")

summary(cox.model)

anova(cox.model,cox.model.trt)

#Other analysis

risk.factors <- cbind(high.bili,albumin,prothrombin,bilirubin,age,edema)
risk.score = rep(0,312)
for(i in 1:312){
  risk.score[i] = sum(cox.model$coefficients*risk.factors[i,])
}
risk.score

sv.highrisk = Surv(time[risk.score>3],delta[risk.score>3])
sv.lowrisk = Surv(time[risk.score<=3],delta[risk.score<=3])


fithigh = survfit(surv.highrisk~1)
fitlow = survfit(surv.lowrisk~1)
band.1 = km.ci(fithigh,method="logep")
plot(band.1,main="Survival Function: High Risk")
band.2 = km.ci(fitlow,method="logep")
plot(band.2,main="Survival Function: Low Risk")


testhigh = survdiff(sv.highrisk~treatment[risk.score>3])
testlow = survdiff(sv.lowrisk~treatment[risk.score<=3])
