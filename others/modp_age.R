
library(tidyverse)
library(keras)
use_condaenv("r-reticulate", required = TRUE)

load("/Users/pol/Dropbox/gmv_project/data/22Q2/multiomics_data_processed_all.RData")
load("/Users/pol/Dropbox/modp/data/multiomics_train.RData")

yvar_train <- sample_info %>%
  dplyr::filter(DepMap_ID %in% expression_train$DepMap_ID) %>% 
  dplyr::select(DepMap_ID, cell_name, age, lineage, lineage_subtype) %>% 
  dplyr::mutate(age = as.numeric(age)) %>% 
  dplyr::filter(!is.na(age))

####

mirna_train <- mirna_train %>% 
  dplyr::filter(DepMap_ID %in% yvar_train$DepMap_ID) %>% 
  tibble::column_to_rownames("DepMap_ID") %>%
  as.matrix()

metabolomics_train <- metabolomics_train %>%
  dplyr::filter(DepMap_ID %in% yvar_train$DepMap_ID) %>% 
  tibble::column_to_rownames("DepMap_ID") %>%
  as.matrix()

methylation_train <- methylation_train %>% 
  dplyr::filter(DepMap_ID %in% yvar_train$DepMap_ID) %>% 
  tibble::column_to_rownames("DepMap_ID") %>%
  as.matrix()

expression_train <- expression_train %>% 
  dplyr::filter(DepMap_ID %in% yvar_train$DepMap_ID) %>% 
  tibble::column_to_rownames("DepMap_ID") %>%
  as.matrix()

achilles_train <- achilles_train %>% 
  dplyr::filter(DepMap_ID %in% yvar_train$DepMap_ID) %>% 
  tibble::column_to_rownames("DepMap_ID") %>%
  as.matrix()

# TEST --------------------------
load("/Users/pol/Dropbox/modp/data/multiomics_test.RData")

yvar_test <- sample_info %>%
  dplyr::filter(DepMap_ID %in% expression_test$DepMap_ID) %>% 
  dplyr::select(DepMap_ID, cell_name, age, lineage, lineage_subtype) %>% 
  dplyr::mutate(age = as.numeric(age)) %>% 
  dplyr::filter(!is.na(age))

mirna_test <- mirna_test %>% 
  dplyr::filter(DepMap_ID %in% yvar_test$DepMap_ID) %>% 
  tibble::column_to_rownames("DepMap_ID") %>%
  as.matrix()

metabolomics_test <- metabolomics_test %>%
  dplyr::filter(DepMap_ID %in% yvar_test$DepMap_ID) %>% 
  tibble::column_to_rownames("DepMap_ID") %>%
  as.matrix()

methylation_test <- methylation_test %>% 
  dplyr::filter(DepMap_ID %in% yvar_test$DepMap_ID) %>% 
  tibble::column_to_rownames("DepMap_ID") %>%
  as.matrix()

expression_test <- expression_test %>% 
  dplyr::filter(DepMap_ID %in% yvar_test$DepMap_ID) %>% 
  tibble::column_to_rownames("DepMap_ID") %>%
  as.matrix()

achilles_test <- achilles_test %>% 
  dplyr::filter(DepMap_ID %in% yvar_test$DepMap_ID) %>% 
  tibble::column_to_rownames("DepMap_ID") %>%
  as.matrix()

##

omics <- list(t(expression_train), t(methylation_train), t(mirna_train), t(metabolomics_train)) 
fctrs <- ubmi(omics)

# mod <- xgboost::xgboost(data = cbind(expression_train, methylation_train, mirna_train, metabolomics_train),
#                         label = yvar_train$age,
#                         lambda = 1, eta = 0.3, gamma = 100, max_depth = 10, subsample = 0.95,
#                         nrounds = 10,
#                         verbose = FALSE,
#                         nthread = parallel::detectCores() - 2,
#                         early_stopping_rounds = 8)

mod <- xgboost::xgboost(data = as.matrix(fctrs@factors[,1:2]),
                        label = yvar_train$age,
                        lambda = 1, eta = 0.3, gamma = 100, max_depth = 10, subsample = 0.95,
                        nrounds = 10,
                        verbose = FALSE,
                        nthread = parallel::detectCores() - 2,
                        early_stopping_rounds = 8)

# shap_values <- SHAPforxgboost::shap.values(xgb_model = mod, X_train = x)
# 
# shap_contrib <- shap_values$shap_score
# ranked_col <- names(colMeans(abs(shap_contrib))[order(colMeans(abs(shap_contrib)), decreasing = TRUE)])
# 
# shap_contrib_df <- data.frame(shap_contrib)
# shap_contrib_nonzero <- shap_contrib_df[, apply(shap_contrib_df, 2, function(x) !all(x == 0, na.rm = TRUE))]

omics_test <- list(t(expression_test), t(methylation_test), t(mirna_test), t(metabolomics_test)) 
fctrs_test <- ubmi(omics_test)

pred <- predict(mod, as.matrix(fctrs_test@factors[,1:2]))

cor(yvar_test$age, pred)
plot(yvar_test$age, pred)

# MODEL --------------------------
## Expression
expression_input <- 
  layer_input(shape = ncol(expression_train), name = "main_input"
  )

expression_output <- expression_input %>% 
  layer_dense(units = 256,
              activation = "relu",
              kernel_initializer = "random_normal"
  ) %>%
  layer_dense(units = 128,
              activation = "relu",
              kernel_initializer = "random_normal"
  ) %>%
  layer_dense(units = 64,
              activation = "relu",
              kernel_initializer = "random_normal"
  )

auxiliary_output <- expression_output %>% 
  layer_dense(units = 1, 
              activation = "linear",
              kernel_initializer = "random_normal",
              name = "auxiliary_output")

## Metabolomics
metabolomics_input <-
  layer_input(shape = ncol(metabolomics_train), name = "metabolomics_input"
  )

metabolomics_output <- metabolomics_input %>% 
  layer_dense(units = 64,
              activation = "relu",
              kernel_initializer = "random_normal"
  )

## miRNA
mirna_input <- 
  layer_input(shape = ncol(mirna_train), name = "mirna_input"
  )

mirna_output <- mirna_input %>% 
  layer_dense(units = 64,
              activation = "relu",
              kernel_initializer = "random_normal"
  )

## Methylation
methylation_input <- 
  layer_input(shape = ncol(methylation_train), name = "methylation_input"
  )

methylation_output <- methylation_input %>% 
  layer_dense(units = 256,
              activation = "relu",
              kernel_initializer = "random_normal"
  ) %>%
  layer_dense(units = 128,
              activation = "relu",
              kernel_initializer = "random_normal"
  ) %>%
  layer_dense(units = 64,
              activation = "relu",
              kernel_initializer = "random_normal"
  )

## Output layer
main_output <-
  layer_concatenate(c(expression_output, 
                      mirna_output,
                      methylation_output,
                      metabolomics_output #,
                      # metadata_output
  ), 
  axis = -1
  ) %>% 
  layer_dense(units = 64, 
              activation = "relu", 
              kernel_initializer = "random_normal"
  ) %>% 
  layer_dense(units = 64, 
              activation = "relu", 
              kernel_initializer = "random_normal"
  ) %>% 
  layer_dropout(0.7) %>% # 0.1
  layer_dense(units = 1, 
              activation = "linear",
              kernel_initializer = "random_normal",
              name = "main_output")

## Final model
model <- keras_model(inputs = c(expression_input,
                                mirna_input,
                                methylation_input,
                                metabolomics_input #,
                                # metadata_input
), 
outputs = c(main_output,
            auxiliary_output
))

## Compile model
summary(model)
# deepviz::plot_model(model)

model <- model %>%
  compile(
    optimizer = optimizer_rmsprop(), # "adam"
    loss = "mse",
    loss_weights = list(main_output = 1.0, auxiliary_output = 0.2),
    metrics = list("mean_absolute_error")
  )

history <- model %>%
  fit(
    x = list(main_input = expression_train, 
             mirna_input = mirna_train,
             methylation_input = methylation_train,
             metabolomics_input = metabolomics_train#,
             # metadata_input = cell_metadata_train
    ),
    y = list(main_output = yvar_train$age, 
             auxiliary_output = yvar_train$age
    ),
    validation_split = 0.2,
    epochs = 50)

# Save model
# save_model_hdf5(model, "multiomics/multiomics_model.h5")

# PREDICTION --------------------------
test_pred <- model %>%
  predict(list(expression_test,
               mirna_test,
               methylation_test,
               metabolomics_test#,
               # cell_metadata_test
  )
  )

plot(test_pred[[1]], yvar_test$age)


