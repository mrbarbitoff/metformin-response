setwd("/media/barbitoff/DATA/Working issues/WES/Metformin_GWAS")

library(ggplot2)
library(reshape2)
library(stringr)
library(randomForest)
library(caret)
library(glmnet)
suppressMessages(library(rattle))
library(pROC)
library(ggplotify)
library(lattice)
library(colorRamps)
library(RColorBrewer)
set.seed(33)

load(".RData")
#### Revamp

gt_data = read.table('metformin_genotypes_new.csv', 
                     sep='\t', header=T, 
                     stringsAsFactors = F, 
                     na.strings = list('-', '', 'NA'))
head(gt_data)
str(gt_data)
#gt_data = na.omit(gt_data)
gt_data$metformin_response = as.factor(gt_data$response)
gt_data$therapy_type = as.factor(gt_data$therapy_type)
gt_data$mono_vs_no = as.factor(ifelse(as.character(gt_data$metformin_response) == 'yes' &
                                        as.character(gt_data$therapy_type) == 'monotherapy', 'yes', 
                                      ifelse(as.character(gt_data$metformin_response) == 'no', 'no', NA)))

gt_data$status = as.factor(gt_data$status)
gt_data$ATM_rs11212617 = str_count(gt_data$ATM_rs11212617, 'C')
gt_data$SLC22A1_rs628031 = str_count(gt_data$SLC22A1_rs628031, 'G')
gt_data$SLC22A1_rs12208357 = str_count(gt_data$SLC22A1_rs12208357, 'C')
gt_data$SLC2A2_rs8192675 = str_count(gt_data$SLC2A2_rs8192675, 'G')
gt_data$SLC47A1_rs2289669 = str_count(gt_data$SLC47A1_rs2289669, 'G')

head(gt_data)

### Case vs. control GLMs

fit_case_control = glm(status ~ ATM_rs11212617 +
                         SLC22A1_rs628031 + SLC22A1_rs12208357 +
                         SLC47A1_rs2289669 + SLC2A2_rs8192675,
                       gt_data, family = "binomial")
summary(fit_case_control)

fit_case_control_covar = glm(status ~ ATM_rs11212617 +
                         SLC22A1_rs628031 + SLC22A1_rs12208357 +
                         SLC47A1_rs2289669 + SLC2A2_rs8192675 +
                         age + sex + BMI + diabetic_inheritance,
                       gt_data, family = "binomial")
summary(fit_case_control_covar)


### Setting up ML settings with caret

train_control <- trainControl(method="cv", 
                              number=3,
                              classProbs=T,
                              summaryFunction=twoClassSummary, 
                              savePredictions = T)
train_data = downSample(x = gt_data,
                        y = gt_data$metformin_response)


model_onesnp = train(metformin_response ~ SLC22A1_rs12208357,
                data=train_data, 
                trControl=train_control, 
                method="glm", na.action = 'na.omit')
print(model_onesnp)
print(model_onesnp$results)


model_1 = train(metformin_response ~ ATM_rs11212617 +
                  SLC22A1_rs628031 + SLC22A1_rs12208357 +
                  SLC47A1_rs2289669 + SLC2A2_rs8192675,
                data=train_data, 
                trControl=train_control, 
                method="glm", na.action = 'na.omit')
print(model_1)
print(model_1$results)

model_2 = train(metformin_response ~ ATM_rs11212617 +
                  SLC22A1_rs628031 + SLC22A1_rs12208357 +
                  SLC47A1_rs2289669 + SLC2A2_rs8192675,
                data=train_data, 
                trControl=train_control, 
                method="rf", na.action = 'na.omit')
print(model_2)
print(model_2$results)


model_3 = train(metformin_response ~ ATM_rs11212617 +
                  SLC22A1_rs628031 + SLC22A1_rs12208357 +
                  SLC47A1_rs2289669 + SLC2A2_rs8192675 +
                  age + sex + BMI + diabetic_inheritance,
                data=train_data, 
                trControl=train_control, 
                method="rf", na.action = 'na.omit')
print(model_3)
print(model_3$results)

selectedIndices <- model_3$pred$mtry == 9

pROC_obj <- roc(model_3$pred$obs[selectedIndices], 
                model_3$pred$yes[selectedIndices],
                smoothed = TRUE,
                # arguments for ci
                ci=TRUE, ci.alpha=0.9, stratified=FALSE,
                # arguments for plot
                plot=TRUE, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE,
                print.auc=TRUE, show.thres=TRUE)

sens.ci <- ci.se(pROC_obj)
plot(sens.ci, type="shape", col="lightblue")

# Validation on multiple data sub-samples


aucs_validation <- data.frame(replicate=1:1000, snp=rep(NA, 1000), 
    all=rep(NA, 1000), random_5=rep(NA, 1000), 
    random_9=rep(NA, 1000), glm_slc22a1=rep(NA, 1000))
varimps <- data.frame(feature=rownames(varImp(model_3)$importance))
varimps_snp <- data.frame(feature=rownames(varImp(model_2)$importance))
for (i in 1:1000){
  train_data_reshuf = downSample(x = gt_data,
                                 y = gt_data$metformin_response)
  for (j in 1:9) {
    train_data_reshuf[, paste0('rv_', j)] = 
      rnorm(nrow(train_data_reshuf))
  }
  model_x = train(metformin_response ~ ATM_rs11212617 +
                    SLC22A1_rs628031 + SLC22A1_rs12208357 +
                    SLC47A1_rs2289669 + SLC2A2_rs8192675,
                  data=train_data_reshuf, 
                  trControl=train_control, 
                  method="glm", na.action = 'na.omit')
  
  model_y = train(metformin_response ~ ATM_rs11212617 +
                    SLC22A1_rs628031 + SLC22A1_rs12208357 +
                    SLC47A1_rs2289669 + SLC2A2_rs8192675 +
                    age + sex + BMI + diabetic_inheritance,
                  data=train_data_reshuf, 
                  trControl=train_control, 
                  method="glm", na.action = 'na.omit')

  model_w = train(metformin_response ~ rv_1 + rv_2 +
                    rv_3 + rv_4 + rv_5,
                  data=train_data_reshuf, 
                  trControl=train_control, 
                  method="glm", na.action = 'na.omit')
  
  model_z = train(metformin_response ~ rv_1 + rv_2 +
                    rv_3 + rv_4 + rv_5 + rv_6 + 
                    rv_7 + rv_8 + rv_9,
                  data=train_data_reshuf, 
                  trControl=train_control, 
                  method="glm", na.action = 'na.omit')
  
  model_t <- train(metformin_response ~ SLC22A1_rs12208357,
                  data=train_data_reshuf, 
                  trControl=train_control, 
                  method="glm", na.action = 'na.omit')
  
  aucs_validation[i, ] = c(i, model_x$results$ROC[3], 
                      model_y$results$ROC[3],
                      model_w$results$ROC[3],
                      model_z$results$ROC[3],
                      model_t$results$ROC)
  varimps[, paste0('rep_', i)] = varImp(model_y)$importance$Overall
  varimps_snp[, paste0('rep_', i)] = 
    varImp(model_x)$importance$Overall
}

apply(aucs_validation[, 2:5], 2, function(elt) table(elt > 0.5))
colMeans(aucs_validation)

val_data <- melt(aucs_validation, id.vars=c('replicate'),
                 measure.vars = c('snp', 'all', 'random_9',
                                  'glm_slc22a1'))
val_data$variable = rep(c('SNPs', 'All', 'Random', 'SLC22A1'), 
                        each=1000)
ggplot(val_data, aes(x=variable, y=value, fill=variable)) +
  geom_violin(scale='width') + 
  geom_boxplot(width=0.2, fill='white', outlier.shape=NA) +
  theme_bw() + scale_fill_brewer(palette = "YlOrBr") +
  guides(fill=F) + xlab('Predictor group') +
  ylab('Area under curve (AUC)')

rownames(varimps) = varimps$feature

varimps_df = data.frame(feature = varimps$feature, 
        mda_mean = rowMeans(varimps[, 2:1001]),
        mda_sd = apply(varimps[, 2:1001], 1, sd))
varimps_df$lower = varimps_df$mda_mean - varimps_df$mda_sd
varimps_df$upper = varimps_df$mda_mean + varimps_df$mda_sd
varimps_df$lower[varimps_df$lower < 0] = 0
varimps_df$upper[varimps_df$upper > 100] = 100
varimps_df$feature = factor(varimps_df$feature, 
        levels=varimps_df$feature[order(varimps_df$mda_mean)])


ggplot(varimps_df, aes(x=feature, y=mda_mean, fill=feature)) +
  geom_bar(col='black', stat='identity') + 
  geom_errorbar(aes(ymin=lower, ymax=upper),
                width=0.5) +
  theme_bw() + guides(fill=F) +
  theme(axis.text.x=element_text(angle=45, hjust=1)) + 
  scale_fill_brewer(palette = 'PuRd') +
  xlab('Variable') + ylab('Importance')


varimps_good <- varimps[, which(aucs_validation$all > 0.65) + 1]

varimps_g_df = data.frame(feature = varimps$feature, 
                        mda_mean = rowMeans(varimps_good),
                        mda_sd = apply(varimps_good, 1, sd))
varimps_g_df$lower = varimps_g_df$mda_mean - varimps_g_df$mda_sd
varimps_g_df$upper = varimps_g_df$mda_mean + varimps_g_df$mda_sd
varimps_g_df$lower[varimps_g_df$lower < 0] = 0
varimps_g_df$upper[varimps_g_df$upper > 100] = 100
varimps_g_df$feature = factor(varimps_g_df$feature, 
        levels=varimps_g_df$feature[order(varimps_g_df$mda_mean)])


ggplot(varimps_g_df, aes(x=feature, y=mda_mean, fill=feature)) +
  geom_bar(col='black', stat='identity') + 
  geom_errorbar(aes(ymin=lower, ymax=upper),
                width=0.5) +
  theme_bw() + guides(fill=F) +
  theme(axis.text.x=element_text(angle=45, hjust=1)) + 
  scale_fill_brewer(palette = 'PuRd') +
  xlab('Variable') + ylab('Importance')


varimps_snp_df = data.frame(feature = varimps_snp$feature, 
                  mda_mean = rowMeans(varimps_snp[, 2:1001]),
                  mda_sd = apply(varimps_snp[, 2:1001], 1, sd))
varimps_snp_df$lower = varimps_snp_df$mda_mean - varimps_snp_df$mda_sd
varimps_snp_df$upper = varimps_snp_df$mda_mean + varimps_snp_df$mda_sd
varimps_snp_df$lower[varimps_snp_df$lower < 0] = 0
varimps_snp_df$upper[varimps_snp_df$upper > 100] = 100
varimps_snp_df$feature = factor(varimps_snp_df$feature, 
    levels=varimps_snp_df$feature[order(varimps_snp_df$mda_mean)])


ggplot(varimps_snp_df, aes(x=feature, y=mda_mean, fill=feature)) +
  geom_bar(col='black', stat='identity') + 
  geom_errorbar(aes(ymin=lower, ymax=upper),
                width=0.5) +
  theme_bw() + guides(fill=F) +
  theme(axis.text.x=element_text(angle=45, hjust=1)) + 
  scale_fill_brewer(palette = 'PuRd') +
  xlab('Variable') + ylab('Importance')

# Variable correlation
cm <- cor(gt_data[, 5:13], use = "pairwise.complete.obs")


drawCoolHM = function(df){
  e = round(df, digits=2)
  myPanel_a <- function(x, y, z, ...) {
    panel.levelplot(x,y,z,...)
    panel.text(x, y,  e[cbind(x,y)]) ## use handy matrix indexing
  }
  return(levelplot(df, col.regions=jet, 
                   at=seq(-0.1, 1, length.out=100), 
                   aspect='fill', colorkey=list(width=3, labels=list(cex=1.0)),
                   scales=list(x=list(rot=45)), xlab=list(label=''), 
                   ylab=list(label=''), panel=myPanel_a))
}

mypal = colorRampPalette(c('white', '#ecad2f'))
mypal = colorRampPalette(c('white', '#d50f0dff'))
jet = mypal(100)
drawCoolHM(cm)



######################################################
################ More phenotypes #####################
######################################################

gt_data = read.table('metformin_with_more.csv', 
                     sep='\t', header=T, 
                     stringsAsFactors = F, 
                     na.strings = list('-', '', 'NA'))
head(gt_data)
str(gt_data)
#gt_data = na.omit(gt_data)
gt_data$response = as.factor(gt_data$response)
gt_data$therapy_type = as.factor(gt_data$therapy_type)

gt_data$ATM_rs11212617 = str_count(gt_data$ATM_rs11212617, 'C')
gt_data$SLC22A1_rs628031 = str_count(gt_data$SLC22A1_rs628031, 'G')
gt_data$SLC22A1_rs12208357 = str_count(gt_data$SLC22A1_rs12208357, 'C')
gt_data$SLC2A2_rs8192675 = str_count(gt_data$SLC2A2_rs8192675, 'G')
gt_data$SLC47A1_rs2289669 = str_count(gt_data$SLC47A1_rs2289669, 'G')

head(gt_data)

apply(gt_data, 2, function(x) table(is.na(x)))
gt_data$WHR[gt_data$WHR == 15] = 1.5


train_control <- trainControl(method="cv", 
                              number=4,
                              classProbs=T,
                              summaryFunction=twoClassSummary, 
                              savePredictions = T)
train_data = downSample(x = gt_data,
                        y = gt_data$response)
train_data = train_data[!(is.na(train_data$response)), ]


model_onesnp = train(response ~ SLC22A1_rs12208357,
                     data=train_data, 
                     trControl=train_control, 
                     method="glm", na.action = 'na.omit')
print(model_onesnp)
print(model_onesnp$results)


model_all = train(response ~ ATM_rs11212617 + SLC22A1_rs628031 +
                    SLC22A1_rs12208357 + SLC47A1_rs2289669 +
                    SLC2A2_rs8192675 + age + sex + BMI +
                    WHR + creatinine + HbA1C + fasting_glucose +
                    diabetic_inheritance,
                     data=train_data, 
                     trControl=train_control, 
                     method="glm", na.action = 'na.omit',)
print(model_all)
print(model_all$results)

varImp(model_all)

model_no_bc = train(response ~ ATM_rs11212617 + SLC22A1_rs628031 +
                    SLC22A1_rs12208357 + SLC47A1_rs2289669 +
                    SLC2A2_rs8192675 + age + sex + BMI +
                    WHR + diabetic_inheritance + creatinine,
                  data=train_data, 
                  trControl=train_control, 
                  method="glm", na.action = 'na.omit',)
print(model_no_bc)
print(model_no_bc$results)


aucs_validation <- data.frame(replicate=1:1000, snp=rep(NA, 1000), 
                              all=rep(NA, 1000))
varimps <- data.frame(feature=rownames(varImp(model_all)$importance))

for (i in 1:1000){
  train_data_reshuf = downSample(x = gt_data,
                                 y = gt_data$response)

  model_x = train(response ~ ATM_rs11212617 + SLC22A1_rs628031 +
                      SLC22A1_rs12208357 + SLC47A1_rs2289669 +
                      SLC2A2_rs8192675 + age + sex + BMI +
                      WHR + creatinine + HbA1C + fasting_glucose +
                      diabetic_inheritance,
                    data=train_data_reshuf, 
                    trControl=train_control, 
                    method="glm", na.action = 'na.omit',)
  
  model_t <- train(response ~ SLC22A1_rs12208357,
                   data=train_data_reshuf, 
                   trControl=train_control, 
                   method="glm", na.action = 'na.omit')
  
  aucs_validation[i, ] = c(i, model_x$results$ROC,
                           model_t$results$ROC)
  varimps[, paste0('rep_', i)] = varImp(model_x)$importance$Overall
}

val_data <- melt(aucs_validation, id.vars=c('replicate'),
                 measure.vars = c('snp', 'all'))
val_data$variable = rep(c('All', 'rs12208357'), 
                        each=1000)
ggplot(val_data, aes(x=variable, y=value, fill=variable)) +
  geom_violin(scale='width') + 
  geom_boxplot(width=0.2, fill='white', outlier.shape=NA) +
  theme_bw() + scale_fill_brewer(palette = "YlOrBr") +
  guides(fill=F) + xlab('Predictor group') +
  ylab('Area under curve (AUC)')


varimps_df = data.frame(feature = varimps$feature, 
                        mda_mean = rowMeans(varimps[, 2:1001]),
                        mda_sd = apply(varimps[, 2:1001], 1, sd))
varimps_df$lower = varimps_df$mda_mean - varimps_df$mda_sd
varimps_df$upper = varimps_df$mda_mean + varimps_df$mda_sd
varimps_df$lower[varimps_df$lower < 0] = 0
varimps_df$upper[varimps_df$upper > 100] = 100
varimps_df$feature = factor(varimps_df$feature, 
                            levels=varimps_df$feature[order(varimps_df$mda_mean)])


brewfun <- colorRampPalette(brewer.pal(9,"RdYlGn"))
importance_pal <- brewfun(13)

ggplot(varimps_df, aes(x=feature, y=mda_mean, fill=feature)) +
  geom_bar(col='black', stat='identity') + 
  geom_errorbar(aes(ymin=lower, ymax=upper),
                width=0.5) +
  theme_bw() + guides(fill=F) +
  theme(axis.text.x=element_text(angle=45, hjust=1)) + 
  scale_fill_manual(values=importance_pal) +
  xlab('Variable') + ylab('Importance')



cm <- cor(gt_data[, 5:17], use = "pairwise.complete.obs")


drawCoolHM = function(df){
  e = round(df, digits=2)
  myPanel_a <- function(x, y, z, ...) {
    panel.levelplot(x,y,z,...)
    panel.text(x, y,  e[cbind(x,y)]) ## use handy matrix indexing
  }
  return(levelplot(df, col.regions=jet, 
                   at=seq(-0.1, 1, length.out=100), 
                   aspect='fill', colorkey=list(width=3, labels=list(cex=1.0)),
                   scales=list(x=list(rot=45)), xlab=list(label=''), 
                   ylab=list(label=''), panel=myPanel_a))
}

mypal = colorRampPalette(c('white', '#ecad2f'))
mypal = colorRampPalette(c('white', '#d50f0dff'))
jet = mypal(100)
drawCoolHM(cm)


# Trying lasso

attach(na.omit(train_data))
mm <- model.matrix(response ~ ATM_rs11212617 + SLC22A1_rs12208357 +
        SLC22A1_rs628031 + SLC2A2_rs8192675 +
        SLC47A1_rs2289669 + sex + diabetic_inheritance)[, -1]

x <- as.matrix(data.frame(age, BMI, WHR, creatinine, mm))

glmmod <- glmnet(x, y=response, 
                 alpha=1, family="binomial")

plot(glmmod, xvar="lambda")

cv.glmmod <- cv.glmnet(x, y=as.numeric(response), alpha=1)
plot(cv.glmmod)

coef(glmmod)[, 2]
detach(na.omit(train_data))

train_data$predict_glmnet = 0.3068210774 + 
  (0.0748807518 * train_data$SLC22A1_rs12208357) +
  (0.0008889572 * train_data$age)

pROC_obj <- roc(train_data$response, 
                train_data$predict_glmnet,
                smoothed = TRUE,
                # arguments for ci
                ci=TRUE, ci.alpha=0.9, stratified=FALSE,
                # arguments for plot
                plot=TRUE, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE,
                print.auc=TRUE, show.thres=TRUE)

sens.ci <- ci.se(pROC_obj)
plot(sens.ci, type="shape", col="lightblue")

lasso_auc <- c()
for (i in 1:1000) {
  train_data_reshuf = downSample(x = gt_data,
                                 y = gt_data$response)
  train_data_reshuf = train_data_reshuf[
    !(is.na(train_data_reshuf$response)), ]
  train_data_reshuf$predict_glmnet = 0.3068210774 + 
    (0.0748807518 * train_data_reshuf$SLC22A1_rs12208357) +
    (0.0008889572 * train_data_reshuf$age)
  pROC_obj <- roc(train_data_reshuf$response, 
                  train_data_reshuf$predict_glmnet,
                  smoothed = TRUE,
                  # arguments for ci
                  ci=TRUE, ci.alpha=0.9, stratified=FALSE,
                  # arguments for plot
                  plot=TRUE, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE,
                  print.auc=TRUE, show.thres=TRUE)
  lasso_auc <- c(lasso_auc, as.numeric(pROC_obj$auc))
}

aucs_validation$lasso <- lasso_auc
val_data <- melt(aucs_validation, id.vars=c('replicate'),
                 measure.vars = c('snp', 'all', 'lasso'))
val_data$variable = rep(c('rs12208357', 'All', 'Lasso'), 
                        each=1000)
ggplot(val_data, aes(x=variable, y=value, fill=variable)) +
  geom_violin(scale='width') + 
  geom_boxplot(width=0.2, fill='white', outlier.shape=NA) +
  theme_bw() + scale_fill_brewer(palette = "YlOrBr") +
  guides(fill=F) + xlab('Predictor group') +
  ylab('Area under curve (AUC)')


# Comparing GLMs and Lasso using CV, full-dataset AUC & 
# a reshuffled validation dataset

aucs_validation <- data.frame(replicate=1:1000, 
                              snp_cv=rep(NA, 1000),
                              gt_cv=rep(NA, 1000),
                              all_cv=rep(NA, 1000), 
                              subset_cv=rep(NA, 1000),
                              lasso_cv=rep(NA, 1000),
                              four_cv=rep(NA, 1000),
                              snp_total=rep(NA, 1000),
                              gt_total=rep(NA, 1000),
                              all_total=rep(NA, 1000), 
                              subset_total=rep(NA, 1000),
                              lasso_total=rep(NA, 1000),
                              four_total=rep(NA, 1000),
                              snp_shuf=rep(NA, 1000),
                              gt_shuf=rep(NA, 1000),
                              all_shuf=rep(NA, 1000), 
                              subset_shuf=rep(NA, 1000),
                              lasso_shuf=rep(NA, 1000),
                              four_shuf=rep(NA, 1000))
varimps <- data.frame(feature=rownames(varImp(model_no_bc)$importance))
lasso_non_zero <- c()
roc_curves = data.frame()

get_auc_proc <- function(predicted, observed){
  pROC_obj <- roc(observed, 
                  predicted,
                  smoothed = TRUE,
                  # arguments for ci
                  ci=TRUE, ci.alpha=0.9, stratified=FALSE,
                  # arguments for plot
                  plot=FALSE)
  return(pROC_obj)
}

for (i in 1:1000){
  print(i)
  train_data_reshuf = downSample(x = gt_data,
                                 y = gt_data$response)
  train_data_reshuf = train_data_reshuf[
    !(is.na(train_data_reshuf$response)), ]
#  train_data_reshuf <- na.omit(train_data_reshuf)
  
  responsive = rownames(train_data_reshuf)[
    train_data_reshuf$response == 'yes']
  test_data = gt_data[!(rownames(gt_data) %in% responsive), ] 
  train_data_validate = downSample(x = test_data,
                                 y = test_data$response)
  train_data_validate = train_data_validate[
    !(is.na(train_data_validate$response)), ]
#  train_data_validate <- na.omit(train_data_validate)
  
  model_x = train(response ~ ATM_rs11212617 + SLC22A1_rs628031 +
                    SLC22A1_rs12208357 + SLC47A1_rs2289669 +
                    SLC2A2_rs8192675 + age + sex + BMI +
                    WHR + creatinine + HbA1C + fasting_glucose +
                    diabetic_inheritance,
                  data=train_data_reshuf, 
                  trControl=train_control, 
                  method="glm", na.action = 'na.omit')
  train_data_reshuf$glm_all = predict(model_x, 
                                        train_data_reshuf, 
                                        type = "prob",
                                        na.action='na.pass')[, 2]
  glm_all_auc_total <- get_auc_proc(train_data_reshuf$glm_all,
                                   train_data_reshuf$response)$auc
  train_data_validate$glm_all = predict(model_x, 
                                      train_data_validate, 
                                      type = "prob",
                                      na.action='na.pass')[, 2]
  glm_all_ROC_shuf <- get_auc_proc(train_data_validate$glm_all,
                              train_data_validate$response)
  glm_all_auc_shuf <- glm_all_ROC_shuf$auc
  glm_all_curve <- data.frame(sensitivity=glm_all_ROC_shuf$sensitivities,
                              specificity=glm_all_ROC_shuf$specificities,
                              model_type='all')
  
  glm_all_auc_cv = as.numeric(model_x$results$ROC)
  
  model_y = train(response ~ ATM_rs11212617 + SLC22A1_rs628031 +
                    SLC22A1_rs12208357 + SLC47A1_rs2289669 +
                    SLC2A2_rs8192675 + age + sex + BMI +
                    WHR + creatinine + diabetic_inheritance,
                  data=train_data_reshuf, 
                  trControl=train_control, 
                  method="glm", na.action = 'na.omit')
  
  train_data_reshuf$glm_subset <- predict(model_y, 
                                            train_data_reshuf, 
                                            type = "prob",
                                          na.action='na.pass')[, 2]
  glm_subset_auc_total <- get_auc_proc(train_data_reshuf$glm_subset,
                                      train_data_reshuf$response)$auc
  train_data_validate$glm_subset <- predict(model_y, 
                                      train_data_validate, 
                                      type = "prob",
                                      na.action='na.pass')[, 2]

  glm_subset_ROC_shuf <- get_auc_proc(train_data_validate$glm_subset,
                                   train_data_validate$response)
  glm_subset_auc_shuf <- glm_subset_ROC_shuf$auc
  glm_subset_curve <- data.frame(sensitivity=glm_subset_ROC_shuf$sensitivities,
                              specificity=glm_subset_ROC_shuf$specificities,
                              model_type='subset')
  
  glm_subset_auc_cv <- as.numeric(model_y$results$ROC)
#  print(c(as.numeric(glm_subset_auc), 
#          as.numeric(model_y$results$ROC)))
  
  model_z = train(response ~ ATM_rs11212617 + SLC22A1_rs628031 +
                    SLC22A1_rs12208357 + SLC47A1_rs2289669 +
                    SLC2A2_rs8192675,
                  data=train_data_reshuf, 
                  trControl=train_control, 
                  method="glm", na.action = 'na.omit')
  
  train_data_reshuf$glm_gt <- predict(model_z, 
                                          train_data_reshuf, 
                                          type = "prob",
                                          na.action='na.pass')[, 2]
  glm_gt_auc_total <- get_auc_proc(train_data_reshuf$glm_gt,
                                       train_data_reshuf$response)$auc
  train_data_validate$glm_gt <- predict(model_z, 
                                            train_data_validate, 
                                            type = "prob",
                                            na.action='na.pass')[, 2]
  glm_gt_ROC_shuf <- get_auc_proc(train_data_validate$glm_gt,
                                   train_data_validate$response)
  glm_gt_auc_shuf <- glm_gt_ROC_shuf$auc
  glm_gt_curve <- data.frame(sensitivity=glm_gt_ROC_shuf$sensitivities,
                              specificity=glm_gt_ROC_shuf$specificities,
                              model_type='gt')
  
  glm_gt_auc_cv <- as.numeric(model_z$results$ROC)
  
  
  model_t <- train(response ~ SLC22A1_rs12208357,
                   data=train_data_reshuf, 
                   trControl=train_control, 
                   method="glm", na.action = 'na.omit')
  
  train_data_reshuf$glm_rs12208357 <- predict(model_t, 
                                                train_data_reshuf, 
                                                type = "prob",
                                                na.action='na.pass')[, 2]
  glm_rs12208357_auc_total <- get_auc_proc(train_data_reshuf$glm_rs12208357,
                                          train_data_reshuf$response)$auc
  train_data_validate$glm_rs12208357 <- predict(model_t, 
                                          train_data_validate, 
                                          type = "prob",
                                          na.action='na.pass')[, 2]
  glm_rs12208357_ROC_shuf <- get_auc_proc(train_data_validate$glm_rs12208357,
                                   train_data_validate$response)
  glm_rs12208357_auc_shuf <- glm_rs12208357_ROC_shuf$auc
  glm_rs12208357_curve <- data.frame(sensitivity=glm_rs12208357_ROC_shuf$sensitivities,
                              specificity=glm_rs12208357_ROC_shuf$specificities,
                              model_type='rs12208357')
  glm_rs12208357_auc_cv <- as.numeric(model_t$results$ROC)
  
  # attach(train_data_reshuf)
  # mm <- model.matrix(response ~ ATM_rs11212617 + SLC22A1_rs12208357 +
  #                      SLC22A1_rs628031 + SLC2A2_rs8192675 +
  #                      SLC47A1_rs2289669 + sex + diabetic_inheritance)[, -1]
  # x <- as.matrix(data.frame(age, BMI, WHR, creatinine, mm))
  # 
  # glmmod <- glmnet(x, y=response, 
  #                  alpha=1, family="binomial")
  # cv.glmmod <- cv.glmnet(x, y=response, alpha=1,
  #                        family='binomial')
  # mf_mod <- cv.glmmod$nzero > 1
  # detach(train_data_reshuf)
  # 
  # attach(train_data_validate)
  # mm <- model.matrix(response ~ ATM_rs11212617 + SLC22A1_rs12208357 +
  #                      SLC22A1_rs628031 + SLC2A2_rs8192675 +
  #                      SLC47A1_rs2289669 + sex + diabetic_inheritance)[, -1]
  # x <- as.matrix(data.frame(age, BMI, WHR, creatinine, mm))
  # 
  # train_data_validate$glmnet <- predict(glmmod, x,
  #                                     type='response')[
  #           , which(cv.glmmod$lambda == 
  #                     cv.glmmod$lambda[mf_mod][
  #                       which.min(cv.glmmod$cvm[mf_mod])])]
  # 
  # detach(train_data_validate)
  #  glmnet_auc <- get_auc_proc(train_data_validate$glmnet,
  #                             train_data_validate$response)
  
  lassoGrid <- expand.grid(lambda=glmmod$lambda,
                           alpha=c(1))
  model_l = train(response ~ ATM_rs11212617 + SLC22A1_rs628031 +
                    SLC22A1_rs12208357 + SLC47A1_rs2289669 +
                    SLC2A2_rs8192675 + age + sex + BMI +
                    WHR + creatinine + diabetic_inheritance,
                  data=train_data_reshuf, 
                  trControl=train_control, 
                  method="glmnet", na.action = 'na.omit',
                  tuneGrid=lassoGrid)
  
  non_zero_count <- sum(as.logical(coef(model_l$finalModel, 
                      model_l$bestTune$lambda) > 0)) - 1
  lasso_non_zero <- c(lasso_non_zero, non_zero_count)
  
  train_data_reshuf$glmnet <- predict(model_l, 
                                              train_data_reshuf, 
                                              type = "prob",
                                              na.action='na.pass')[, 2]
  glmnet_auc_total <- get_auc_proc(train_data_reshuf$glmnet,
                                           train_data_reshuf$response)$auc
  train_data_validate$glmnet <- predict(model_l, 
                                              train_data_validate, 
                                              type = "prob",
                                              na.action='na.pass')[, 2]
  glmnet_ROC_shuf <- get_auc_proc(train_data_validate$glmnet,
                                          train_data_validate$response)
  glmnet_auc_shuf <- glmnet_ROC_shuf$auc
  glmnet_curve <- data.frame(sensitivity=glmnet_ROC_shuf$sensitivities,
                                     specificity=glmnet_ROC_shuf$specificities,
                                     model_type='lasso')
  
  glmnet_auc_cv<- max(model_l$results$ROC)
  
  ### Cheating - taking 4 best variables from Lasso
  model_w = train(response ~ SLC22A1_rs12208357 + sex +
                    WHR + diabetic_inheritance,
                  data=train_data_reshuf, 
                  trControl=train_control, 
                  method="glm", na.action = 'na.omit')
  
  train_data_reshuf$glm_four <- predict(model_w, 
                                          train_data_reshuf, 
                                          type = "prob",
                                          na.action='na.pass')[, 2]
  glm_four_auc_total <- get_auc_proc(train_data_reshuf$glm_four,
                                       train_data_reshuf$response)$auc
  train_data_validate$glm_four <- predict(model_w, 
                                            train_data_validate, 
                                            type = "prob",
                                            na.action='na.pass')[, 2]
  glm_four_ROC_shuf <- get_auc_proc(train_data_validate$glm_four,
                                          train_data_validate$response)
  glm_four_auc_shuf <- glm_four_ROC_shuf$auc
  glm_four_curve <- data.frame(sensitivity=glm_four_ROC_shuf$sensitivities,
                                     specificity=glm_four_ROC_shuf$specificities,
                                     model_type='four')
  glm_four_auc_cv <- as.numeric(model_w$results$ROC)
  
  
  
  aucs_validation[i, ] = c(i, 
                           glm_rs12208357_auc_cv,
                           glm_gt_auc_cv,
                           glm_all_auc_cv, 
                           glm_subset_auc_cv,
                           glmnet_auc_cv,
                           glm_four_auc_cv,
                           glm_rs12208357_auc_total,
                           glm_gt_auc_total,
                           glm_all_auc_total,
                           glm_subset_auc_total,
                           glmnet_auc_total,
                           glm_four_auc_total,
                           glm_rs12208357_auc_shuf,
                           glm_gt_auc_shuf,
                           glm_all_auc_shuf,
                           glm_subset_auc_shuf,
                           glmnet_auc_shuf,
                           glm_four_auc_shuf)
  varimps[, paste0('rep_', i)] = varImp(model_l)$importance$Overall
  roc_curves_iter <- as.data.frame(rbind(glm_all_curve,
                                         glm_gt_curve,
                                         glm_rs12208357_curve,
                                         glm_subset_curve,
                                         glmnet_curve,
                                         glm_four_curve))
  roc_curves = as.data.frame(rbind(roc_curves, roc_curves_iter))
}

mean(lasso_non_zero)
median(lasso_non_zero)

val_data <- melt(aucs_validation, id.vars=c('replicate'))

ylbrfun <- colorRampPalette(brewer.pal(9,"YlOrBr"))
auc_pal <- ylbrfun(18)


ggplot(val_data, aes(x=variable, y=value, fill=variable)) +
  geom_violin(scale='width') + 
  geom_boxplot(width=0.2, fill='white', outlier.shape=NA) +
  theme_bw() + scale_fill_manual(values=auc_pal) +
  guides(fill=F) + xlab('Predictor group') +
  ylab('Area under curve (AUC)') +
  theme(axis.text.x=element_text(angle=45, hjust=1))

val_data$score_type = sapply(as.character(val_data$variable),
                    function(x) strsplit(x, '_')[[1]][2])
val_data$model_type = sapply(as.character(val_data$variable),
                    function(x) strsplit(x, '_')[[1]][1])

ggplot(val_data, aes(x=model_type, y=value, fill=model_type)) +
  geom_violin(scale='width') + 
  geom_boxplot(width=0.2, fill='white', outlier.shape=NA) +
  theme_bw() + scale_fill_brewer(palette = "YlOrBr") +
  guides(fill=F) + xlab('Predictor group') +
  ylab('Area under curve (AUC)') +
  theme(axis.text.x=element_text(angle=45, hjust=1)) +
  facet_wrap(~score_type, ncol=3)


aggregate(value~variable, val_data, median)
wilcox.test(aucs_validation$snp_cv, aucs_validation$lasso_cv,
            alternative = "less")
wilcox.test(aucs_validation$subset_shuf, 
            aucs_validation$lasso_shuf,
            alternative = "less")

importance_pal <- brewfun(11)

varimps_df = data.frame(feature = varimps$feature, 
                        mda_mean = rowMeans(varimps[, 2:1001],
                                            na.rm = T),
                        mda_sd = apply(varimps[, 2:1001], 1, 
                                       sd, na.rm=T))
varimps_df$mda_se = varimps_df$mda_sd/sqrt(1000)
t_val = abs(distributions3::quantile(distributions3::StudentsT(999), 0.05))
varimps_df$lower = varimps_df$mda_mean - t_val * varimps_df$mda_se
varimps_df$upper = varimps_df$mda_mean + t_val * varimps_df$mda_se
varimps_df$lower[varimps_df$lower < 0] = 0
varimps_df$upper[varimps_df$upper > 100] = 100
varimps_df$feature = factor(varimps_df$feature, 
                            levels=varimps_df$feature[order(varimps_df$mda_mean)])


ggplot(varimps_df, aes(x=feature, y=mda_mean, fill=feature)) +
  geom_bar(col='black', stat='identity') + 
  geom_errorbar(aes(ymin=lower, ymax=upper),
                width=0.5) +
  theme_bw() + guides(fill=F) +
  theme(axis.text.x=element_text(angle=45, hjust=1)) + 
  scale_fill_manual(values=importance_pal) +
  xlab('Variable') + ylab('Importance')

di_in <- as.matrix(aggregate(diabetic_inheritance~response, 
                     gt_data, table)$diabetic_inheritance)
chisq.test(di_in)

wilcox.test(WHR~response, gt_data)
aggregate(WHR~response, gt_data, median)

ggplot(gt_data, aes(x=response, y=WHR)) + geom_violin()

roc_curves$spec_bin = round(roc_curves$specificity/0.05, 0)/20
roc_agg = aggregate(sensitivity~spec_bin+model_type, roc_curves,
                    function(x) round(mean(x), 2))
roc_agg$sens_sd = aggregate(sensitivity~spec_bin+model_type, 
                            roc_curves,
                       function(x) round(sd(x), 3))$sensitivity

endpoints = roc_agg[roc_agg$spec_bin == 1, ]
endpoints$sensitivity = 0
endpoints$sens_sd = 0

roc_agg[roc_agg$spec_bin == 1, 'spec_bin'] = 0.97
roc_agg = rbind(roc_agg, endpoints)

roc_agg$min_env <- roc_agg$sensitivity - roc_agg$sens_sd
roc_agg$max_env <- roc_agg$sensitivity + roc_agg$sens_sd


ggplot(roc_agg, aes(x=1-spec_bin, y=sensitivity,
                    ymin=min_env, ymax=max_env,
                    col=model_type, fill=model_type)) +
  geom_line(lwd=1) + theme_bw() + 
  geom_ribbon(alpha=0.25) + 
  geom_abline(slope=1, lty=2)  +
#  theme(legend.position = c(0.87, 0.4)) + 
  facet_wrap(~model_type, nrow=3, ncol=2) + 
  guides(fill=F, col=F)

# For a set of lambda values, save the coefficients and evaluate
# Area Under Curve (AUCs) as above. Look at perforemcance across
# data subsets


# A final validation of test statistics by building GLM on
# purely random noise

get_noice_aucs <- function(pred_n) {
  noise_validation <- data.frame(replicate=1:1000, 
                                random_cv=rep(NA, 1000),
                                random_shuf=rep(NA, 1000),
                                random_total=rep(NA, 1000))
  
  for (i in 1:1000){
    print(i)
    train_data_reshuf = downSample(x = gt_data,
                                   y = gt_data$response)
    train_data_reshuf = train_data_reshuf[
      !(is.na(train_data_reshuf$response)), ]
    #  train_data_reshuf <- na.omit(train_data_reshuf)
    
    responsive = rownames(train_data_reshuf)[
      train_data_reshuf$response == 'yes']
    test_data = gt_data[!(rownames(gt_data) %in% responsive), ] 
    train_data_validate = downSample(x = test_data,
                                     y = test_data$response)
    train_data_validate = train_data_validate[
      !(is.na(train_data_validate$response)), ]
    #  train_data_validate <- na.omit(train_data_validate)
    
    for (j in 1:pred_n){
      train_data_reshuf[, paste0('noise_', j)] = 
        rnorm(nrow(train_data_reshuf))
      train_data_validate[, paste0('noise_', j)] = 
        rnorm(nrow(train_data_validate))
    }
    
    train_data_reshuf = train_data_reshuf[, 
                      c(2, 19:ncol(train_data_reshuf))]
    train_data_validate = train_data_validate[, 
                      c(2, 19:ncol(train_data_validate))]
    model_n = train(response ~ .,
                    data=train_data_reshuf, 
                    trControl=train_control, 
                    method="glm", na.action = 'na.omit')
    train_data_reshuf$glm_noise = predict(model_n, 
                                        train_data_reshuf, 
                                        type = "prob",
                                        na.action='na.pass')[, 2]
    glm_noise_auc_total <- get_auc_proc(train_data_reshuf$glm_noise,
                                      train_data_reshuf$response)$auc
    train_data_validate$glm_noise = predict(model_n, 
                                          train_data_validate, 
                                          type = "prob",
                                          na.action='na.pass')[, 2]
    glm_noise_ROC_shuf <- get_auc_proc(train_data_validate$glm_noise,
                                     train_data_validate$response)
    glm_noise_auc_shuf <- glm_noise_ROC_shuf$auc
    
    glm_noise_auc_cv = as.numeric(model_n$results$ROC)
    
    noise_validation[i, ] = c(i, glm_noise_auc_cv,
                              glm_noise_auc_shuf,
                              glm_noise_auc_total)
  }
  return(c(colMeans(noise_validation[, 2:4], na.rm=T),
           apply(noise_validation[, 2:4], 2, sd)))
}

noise_means <- data.frame(pred_n = 1:13,
                          random_cv = rep(NA, 13),
                          random_shuf = rep(NA, 13),
                          random_total = rep(NA, 13),
                          random_cv_sd=rep(NA, 13),
                          random_shuf_sd=rep(NA, 13),
                          random_total_sd=rep(NA, 13))
for (n in 1:13){
  noise_means[n, ] = c(n, get_noice_aucs(n))
}

noise_data <- cbind(melt(noise_means[, 1:4], id.vars='pred_n'),
          melt(noise_means[, c(1, 5:7)], id.vars='pred_n')$value)
colnames(noise_data) = c('pred_n', 'score_type', 'AUC',
                         'AUC_sd')
noise_data$ymin = noise_data$AUC - noise_data$AUC_sd
noise_data$ymax = noise_data$AUC + noise_data$AUC_sd

ggplot(noise_data, aes(x=pred_n, y=AUC, col=score_type,
                       ymin=ymin, ymax=ymax, group=score_type)) +
  geom_point(size=2) + geom_line(lwd=1) +
  geom_errorbar(width=0.4, lwd=0.5) +
  theme_bw() + xlab('Number of noise variables') +
  ylab('Area under curve') +
  facet_wrap(~score_type, ncol=3)

#######################################################
############## Old data and old functions #############
#######################################################


# fit_mono_vs_none = glm(mono_vs_no ~ ATM_rs11212617 +
#                          SLC22A1_rs628031 + SLC22A1_rs122083571 +
#                          SLC47A1_rs2289669 + SLC2A2_rs8192675,
#                        gt_data, family = "binomial")
# summary(fit_mono_vs_none)


model_1 = train(metformin_response ~ ATM_rs11212617 +
                  SLC22A1_rs628031 + SLC22A1_rs122083571 +
                  SLC47A1_rs2289669 + SLC2A2_rs8192675,
                data=train_data, 
                trControl=train_control, 
                method="rf", na.action = 'na.omit')
print(model_1)

model_1_glm = train(metformin_response ~ ATM_rs11212617 +
                  SLC22A1_rs628031 + SLC22A1_rs122083571 +
                  SLC47A1_rs2289669 + SLC2A2_rs8192675,
                data=train_data, 
                trControl=train_control, 
                method="glm", na.action = 'na.omit')
print(model_1_glm)

model_1_dt = train(metformin_response ~ ATM_rs11212617 +
                      SLC22A1_rs628031 + SLC22A1_rs122083571 +
                      SLC47A1_rs2289669 + SLC2A2_rs8192675,
                    data=train_data, 
                    trControl=train_control, 
                    tuneLength=10, 
                    method="rpart", na.action = 'na.omit')
print(model_1_dt)
fancyRpartPlot(model_1_dt$finalModel)

model_1_dt_covar = train(metformin_response ~ ATM_rs11212617 +
                     SLC22A1_rs628031 + SLC22A1_rs122083571 +
                     SLC47A1_rs2289669 + SLC2A2_rs8192675 +
                     age + sex + BMI + diabetic_inheritance,
                   data=train_data, 
                   trControl=train_control, 
                   tuneLength=10, 
                   method="rpart", na.action = 'na.omit')
print(model_1_dt_covar)
fancyRpartPlot(model_1_dt_covar$finalModel)


model_1_rf_covar = train(metformin_response ~ ATM_rs11212617 +
                           SLC22A1_rs628031 + SLC22A1_rs122083571 +
                           SLC47A1_rs2289669 + SLC2A2_rs8192675 +
                           age + sex + BMI + diabetic_inheritance,
                         data=train_data, 
                         trControl=train_control, 
                         tuneLength=10, 
                         method="rf", na.action = 'na.omit')
print(model_1_rf_covar)
print(model_1_rf_covar$results)

selectedIndices <- model_1_rf_covar$pred$mtry == 6

plot.roc(model_1_rf_covar$pred$obs[selectedIndices],
         model_1_rf_covar$pred$yes[selectedIndices])

plot.roc(model_1_rf_covar$pred$obs[selectedIndices],
         model_1_rf_covar$pred$yes[selectedIndices])


# model_2 = train(mono_vs_no ~ ATM_rs11212617 +
#                   SLC22A1_rs628031 + SLC22A1_rs122083571 +
#                   SLC47A1_rs2289669 + SLC2A2_rs8192675,
#                 data=gt_data, 
#                 trControl=train_control, 
#                 method="rf", na.action = 'na.omit')
# print(model_2)


fit_metformin_covar = glm(metformin_response ~ ATM_rs11212617 +
                            SLC22A1_rs628031 + SLC22A1_rs122083571 +
                            SLC47A1_rs2289669 + SLC2A2_rs8192675 +
                            age + sex + BMI,
                          gt_data, family = "binomial")
summary(fit_metformin_covar)

# fit_mono_vs_none_covar = glm(mono_vs_no ~ ATM_rs11212617 +
#                                SLC22A1_rs628031 + SLC22A1_rs122083571 +
#                                SLC47A1_rs2289669 + SLC2A2_rs8192675 +
#                                age + sex + BMI,
#                              gt_data, family = "binomial")
# summary(fit_mono_vs_none_covar)


model_1_covar = train(metformin_response ~ ATM_rs11212617 +
                        SLC22A1_rs628031 + SLC22A1_rs122083571 +
                        SLC47A1_rs2289669 + SLC2A2_rs8192675 +
                        age + sex + BMI,
                      data=gt_data, 
                      trControl=train_control, 
                      method="rf", na.action = 'na.omit')
print(model_1_covar)

# model_2_covar = train(mono_vs_no ~ ATM_rs11212617 +
#                         SLC22A1_rs628031 + SLC22A1_rs122083571 +
#                         SLC47A1_rs2289669 + SLC2A2_rs8192675 +
#                         age + sex + BMI,
#                       data=gt_data, 
#                       trControl=train_control, 
#                       method="rf", na.action = 'na.omit')
# print(model_2_covar)
# varImp(model_2_covar)

gt_data$predicted = predict(model_2_covar, gt_data, type='vote')




gt_recessive = gt_data
for (i in 9:13) {
  gt_recessive[, i] = as.numeric(as.numeric(gt_data[, i]) == 2)
}

fit_metformin = glm(metformin_response ~ ATM_rs11212617 +
                      SLC22A1_rs628031 + SLC22A1_rs122083571 +
                      SLC47A1_rs2289669 + SLC2A2_rs8192675,
                    gt_recessive, family = "binomial")
summary(fit_metformin)

model_1_covar = train(metformin_response ~ ATM_rs11212617 +
                        SLC22A1_rs628031 + SLC22A1_rs122083571 +
                        SLC47A1_rs2289669 + SLC2A2_rs8192675 +
                        age + sex + BMI,
                      data=gt_recessive, 
                      trControl=train_control, 
                      method="rf", na.action = 'na.omit')
print(model_1_covar)

model_2_covar = train(mono_vs_no ~ ATM_rs11212617 +
                        SLC22A1_rs628031 + SLC22A1_rs122083571 +
                        SLC47A1_rs2289669 + SLC2A2_rs8192675 +
                        age + sex + BMI,
                      data=gt_recessive, 
                      trControl=train_control, 
                      method="rf", na.action = 'na.omit')
print(model_2_covar)
varImp(model_2_covar)




gt_dominant = gt_data
for (i in 9:13) {
  gt_dominant[, i] = as.numeric(as.numeric(gt_data[, i]) > 1)
}

fit_metformin = glm(metformin_response ~ ATM_rs11212617 +
                      SLC22A1_rs628031 + SLC22A1_rs122083571 +
                      SLC47A1_rs2289669 + SLC2A2_rs8192675,
                    gt_dominant, family = "binomial")
summary(fit_metformin)

model_1_covar = train(metformin_response ~ ATM_rs11212617 +
                        SLC22A1_rs628031 + SLC22A1_rs122083571 +
                        SLC47A1_rs2289669 + SLC2A2_rs8192675 +
                        age + sex + BMI,
                      data=gt_dominant, 
                      trControl=train_control, 
                      method="rf", na.action = 'na.omit')
print(model_1_covar)

model_2_covar = train(mono_vs_no ~ ATM_rs11212617 +
                        SLC22A1_rs628031 + SLC22A1_rs122083571 +
                        SLC47A1_rs2289669 + SLC2A2_rs8192675 +
                        age + sex + BMI,
                      data=gt_dominant, 
                      trControl=train_control, 
                      method="rf", na.action = 'na.omit')
print(model_2_covar)
varImp(model_2_covar)

