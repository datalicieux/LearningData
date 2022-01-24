#######################################
##### PREDICTING VIOLENCE IN MALI #####
#######################################


### Load Packages and Data 
library(here)
library(tidyverse)
library(caret)
library(glmnet)
library(xgboost)
library(doParallel)


load(here("newdata/mali_osv.RData"))


##### DATA PREPARATION ##### 

# Outcome as factor
predosv$violence = factor(predosv$violence, levels = c(1, 0), labels = c("violence", "none"))

# DF sorted by month
predosv = predosv[order(predosv$month), ]

# Number of cells per month
n_cells = length(unique(predosv$gid))

### Split sample by month
trainIndex = as.vector(predosv$month)
trainIndex = trainIndex %in% sort(unique(predosv$month))[1:33]

# Share of training data 
sum(trainIndex)/nrow(predosv)

train = predosv[trainIndex, ] %>% select(-c(gid, month))
test = predosv[!trainIndex, ]



# Training Outcome vector 
y_train = train$violence

# Training Predictor matrix (december as baseline)
x_train = train %>% 
  select(-c(violence, december))


## Setting for Forward chaining (validation):  
# Start with 8 months for training
# Move forward by 6 months 
# Use next 3 months for validation
fitControl = trainControl(method = "timeslice", 
                          initialWindow = n_cells*8, 
                          horizon = n_cells*3, 
                          fixedWindow = FALSE, 
                          skip = (n_cells*3-1), 
                          classProbs = TRUE, 
                          summaryFunction=twoClassSummary)


### PARALLELIZE COMPUTATIONS 
cl <- parallel::makeCluster(6)
registerDoParallel(cl)


### Baseline : Logistic Regression 
set.seed(53)
lr_base = train(x = x_train, y = y_train,
               family = "binomial", 
               metric = "ROC", 
               method = "glm", 
               trControl = fitControl)
lr_base
summary(lr_base)



### Model 1: Lasso 

lasso_tune = expand.grid(alpha = 1, lambda = seq(0.0005, 0.01, by = 0.0005))

cl <- parallel::makeCluster(6)
registerDoParallel(cl)

set.seed(53)
lasso = train(x = x_train, y = y_train,
               metric = "ROC", 
               method = "glmnet",
               family = "binomial",
               trControl = fitControl,
               tuneGrid = lasso_tune,
               preProcess = c("center","scale"))

# Training performance of best tune 
perf_lasso = lasso$results %>% 
  filter(lasso$results$lambda == 
           lasso$bestTune$lambda) %>% 
  mutate(model = "lasso")

ggplot(lasso)


### Model 2: Elastic Net 

elnet_tune = expand.grid(alpha = seq(0.01, 1, by = 0.02), lambda = seq(0.001, 0.1, by = 0.05))

cl <- parallel::makeCluster(6)
registerDoParallel(cl)

set.seed(53)
elnet = train(x = x_train, y = y_train,
              metric = "ROC", 
              method = "glmnet",
              family = "binomial",
              trControl = fitControl,
              tuneGrid = elnet_tune,
              preProcess = c("center","scale"))

# Training performance of best tune 
perf_elnet = elnet$results %>% 
  filter(elnet$results$lambda == 
           elnet$bestTune$lambda & 
           elnet$results$alpha == 
           elnet$bestTune$alpha) %>% 
  mutate(model = "elnet")



### Model 3: XG Boost

xgb_tune = expand.grid(nrounds = c(50, 500, 1000),
                       max_depth = 3,
                       eta = c(0.01, 0.1, 0.2),
                       gamma = 0,
                       colsample_bytree = c(0.5, 0.7),
                       min_child_weight = 5,
                       subsample = c(0.5, 0.8))

cl <- parallel::makeCluster(6)
registerDoParallel(cl)

set.seed(53)
xgb <- train(x = x_train, y = y_train,
                    method = "xgbTree",
                    trControl = fitControl,
                    tuneGrid = xgb_tune, 
                    metric = "ROC")


# Training performance of best tune 
perf_xgb = xgb$results %>% 
  filter(xgb$results$nrounds == 
           xgb$bestTune$nrounds & 
           xgb$results$max_depth == 
           xgb$bestTune$max_depth & 
           xgb$results$eta == 
           xgb$bestTune$eta &
           xgb$results$colsample_bytree == 
           xgb$bestTune$colsample_bytree &
           xgb$results$min_child_weight == 
           xgb$bestTune$min_child_weight &
           xgb$results$subsample == 
           xgb$bestTune$subsample) %>% 
  mutate(model = "xgb")


### Plot Variable importance of XGB
imp_xgb = matrix(varImp(xgb))
imp_xgb = data.frame(imp_xgb[[1]])
imp_xgb$name = rownames(imp_xgb)
imp_xgb = imp_xgb[1:20, ]
  
plot_varimp = ggplot(data = imp_xgb) + 
  geom_bar(aes(x = reorder(name, -rev(Overall)), y = Overall), stat = "identity", fill = "darkblue") + 
  labs(title = "XG Boost Variable Importance (Top 20)", 
       x = NULL, y = "Overall Importance") + 
  theme_minimal() + 
  coord_flip()

ggsave(plot_varimp, filename = here("xgb_varimp.png"))



########## PREDICTION ##########

# Cut-off at 0.5
pred_y = predict(xgb, test, type = "prob") %>% 
  mutate(pred_viol = factor(ifelse(violence >= 0.5, "violence", "none")))

confusionMatrix(pred_y$pred_viol, test$violence)



### ROC Curve (holdout)

# Save predicted probabilities
pred_probs <- data.frame(base = predict(lr_base, test, type = "prob"),
                        xgb = predict(xgb, test, type = "prob"))

# Setup loop 
sensitiv1 <- rep(NA, length(seq(0,1,by =0.01)))
sensitiv2 <- rep(NA, length(seq(0,1,by =0.01)))
names(sensitiv1) <- paste0("prob", seq(0,1,by =0.01))
names(sensitiv2) <- paste0("prob", seq(0,1,by =0.01))

fallout1 <- rep(NA, length(seq(0,1,by =0.01)))
fallout2 <- rep(NA, length(seq(0,1,by =0.01)))
names(fallout1) <- paste0("prob", seq(0,1,by =0.01))
names(fallout2) <- paste0("prob", seq(0,1,by =0.01))

j = 1

# Run loop 
for(prob in seq(0,1,by =0.01)){
  pred1 <- factor(ifelse(pred_probs$base.violence >= prob, "violence", "none"))
  pred2 <- factor(ifelse(pred_probs$xgb.violence >= prob, "violence", "none"))
  
  sensitiv1[j] <- sensitivity(data = pred1, reference = test$violence, positive = "violence")
  sensitiv2[j] <- sensitivity(data = pred2, reference = test$violence, positive = "violence")
  
  fallout1[j] <- 1 - specificity(data = pred1, reference = test$violence, negative = "none")
  fallout2[j] <- 1 - specificity(data = pred2, reference = test$violence, negative = "none")
  
  j = j+1
}

# Save in dataframe
df_roc <- bind_cols(sensitiv1 = sensitiv1, fallout1 = fallout1,
                    sensitiv2 = sensitiv2, fallout2 = fallout2)

# Plot ROC curve
plot_roc = ggplot(df_roc) + 
  geom_line(aes(fallout1, sensitiv1, color = "Base (LR)")) + 
  geom_line(aes(fallout2, sensitiv2, color = "XGB")) + 
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") + 
  labs(x = "Fallout", y = "Sensitivity", title = "ROC Curve") + 
  scale_color_manual(name = "Model",
                     values = c("Base (LR)" = "darkgreen", "XGB" = "orange"))


ggsave(plot_roc, filename = here("prediction_plot_1218.png"))



##### Map prediction performance 

# Bring data to right format (incl. geocodes)

load(here("newdata/mali_grid.RData"))

plot_pred = pred_y %>% rename(prob_viol = violence) %>% 
  mutate(prob_viol = round(prob_viol, digits = 2))

mali_grid = as.data.frame(mali_grid)
perform = cbind(test, plot_pred) %>% select(gid, violence, month, prob_viol) %>% 
  left_join(mali_grid, by = "gid")

perform = st_as_sf(perform)



# Plot for December 2018
perform_12 = perform %>% filter(month == "2018_12")

plot_perform_12 = ggplot(data = perform_12) + 
  geom_sf(aes(fill = prob_viol)) + 
  scale_fill_continuous(low = "white", high = "red", name = "Pred. Violence") + 
  geom_point(aes(x = xcoord, y = ycoord, shape = violence), 
             alpha = .7, size = 1) + 
  scale_shape_manual(values = c(8, 0), labels = c("Violence", "None"), name = "Real Violence") + 
  labs(y = NULL, x = NULL, title = "Figure 1. Predicted and Real Violence in Mali, December 2018.") + 
  theme_minimal() + 
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        legend.key.size = unit(0.5, 'cm'), #change legend key size
        legend.key.height = unit(0.5, 'cm'), #change legend key height
        legend.key.width = unit(0.5, 'cm'), #change legend key width
        legend.title = element_text(size = 8), #change legend title font size
        legend.text = element_text(size = 7),
        plot.title = element_text(size=12), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

plot_perform_12

### Save December Plot 
ggsave(plot_perform_12, filename = here("prediction_plot_1218.png"), dpi = 1000)

