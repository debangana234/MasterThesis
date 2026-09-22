source("package_helpers.R")

run_predictive_models <- function(model_data,
                                  outcome_col = "disease_state",
                                  sample_col = "Sample_Name",
                                  cpg_columns = grep("^cg", colnames(model_data), value = TRUE),
                                  cell_type_columns = c("CD8T", "CD4T", "NK", "Bcell", "Mono", "Gran"),
                                  number_of_cpgs = 200,
                                  train_fraction = 0.8,
                                  seed = 42,
                                  output_dir = NULL) {
  install_if_missing(packages = c("caret", "e1071", "pROC", "ggplot2"))
  library(caret)
  library(e1071)
  library(pROC)
  library(ggplot2)

  set.seed(seed)

  model_data <- model_data[, !colnames(model_data) %in% c("", "X")]
  model_data[[outcome_col]][model_data[[outcome_col]] == "control"] <- "Control"
  model_data[[outcome_col]] <- factor(model_data[[outcome_col]], levels = c("Control", "Parkinson"))

  cpg_variances <- apply(model_data[, cpg_columns], 2, var, na.rm = TRUE)
  selected_cpgs <- names(sort(cpg_variances, decreasing = TRUE)[seq_len(number_of_cpgs)])
  selected_features <- c(selected_cpgs, cell_type_columns)

  train_index <- createDataPartition(model_data[[outcome_col]], p = train_fraction, list = FALSE)
  train_data <- model_data[train_index, ]
  test_data <- model_data[-train_index, ]

  x_train <- train_data[, selected_features]
  x_test <- test_data[, selected_features]
  y_train <- train_data[[outcome_col]]
  y_test <- test_data[[outcome_col]]

  preprocess_steps <- preProcess(x_train, method = c("medianImpute", "center", "scale"))
  x_train_scaled <- predict(preprocess_steps, x_train)
  x_test_scaled <- predict(preprocess_steps, x_test)

  glm_data <- data.frame(disease_state = y_train, x_train_scaled)
  logistic_model <- glm(disease_state ~ ., data = glm_data, family = "binomial")
  logistic_prob <- predict(logistic_model, newdata = as.data.frame(x_test_scaled), type = "response")
  logistic_pred <- ifelse(logistic_prob >= 0.5, "Parkinson", "Control")
  logistic_pred <- factor(logistic_pred, levels = c("Control", "Parkinson"))

  svm_linear <- svm(x_train_scaled, y_train, kernel = "linear", probability = TRUE)
  svm_linear_pred <- predict(svm_linear, x_test_scaled, probability = TRUE)
  svm_linear_prob <- attr(svm_linear_pred, "probabilities")[, "Parkinson"]

  svm_radial <- svm(x_train_scaled, y_train, kernel = "radial", probability = TRUE)
  svm_radial_pred <- predict(svm_radial, x_test_scaled, probability = TRUE)
  svm_radial_prob <- attr(svm_radial_pred, "probabilities")[, "Parkinson"]

  logistic_metrics <- confusionMatrix(logistic_pred, y_test, positive = "Parkinson")
  svm_linear_metrics <- confusionMatrix(svm_linear_pred, y_test, positive = "Parkinson")
  svm_radial_metrics <- confusionMatrix(svm_radial_pred, y_test, positive = "Parkinson")

  logistic_roc <- roc(y_test, logistic_prob, levels = c("Control", "Parkinson"), direction = "<", quiet = TRUE)
  svm_linear_roc <- roc(y_test, svm_linear_prob, levels = c("Control", "Parkinson"), direction = "<", quiet = TRUE)
  svm_radial_roc <- roc(y_test, svm_radial_prob, levels = c("Control", "Parkinson"), direction = "<", quiet = TRUE)

  metrics_table <- data.frame(
    model = c("logistic", "svm_linear", "svm_radial"),
    accuracy = c(logistic_metrics$overall["Accuracy"],
      svm_linear_metrics$overall["Accuracy"],
      svm_radial_metrics$overall["Accuracy"]),
    sensitivity = c(logistic_metrics$byClass["Sensitivity"],
      svm_linear_metrics$byClass["Sensitivity"],
      svm_radial_metrics$byClass["Sensitivity"]),
    specificity = c(logistic_metrics$byClass["Specificity"],
      svm_linear_metrics$byClass["Specificity"],
      svm_radial_metrics$byClass["Specificity"]),
    AUC = c(as.numeric(auc(logistic_roc)),
      as.numeric(auc(svm_linear_roc)),
      as.numeric(auc(svm_radial_roc))))

  predictions <- data.frame(
    Sample_Name = test_data[[sample_col]],
    actual = y_test,
    logistic_probability = logistic_prob,
    svm_linear_probability = svm_linear_prob,
    svm_radial_probability = svm_radial_prob,
    logistic_prediction = logistic_pred,
    svm_linear_prediction = svm_linear_pred,
    svm_radial_prediction = svm_radial_pred
  )

  logistic_confusion <- as.data.frame(logistic_metrics$table)
  logistic_confusion$model <- "logistic"

  svm_linear_confusion <- as.data.frame(svm_linear_metrics$table)
  svm_linear_confusion$model <- "svm_linear"

  svm_radial_confusion <- as.data.frame(svm_radial_metrics$table)
  svm_radial_confusion$model <- "svm_radial"

  confusion_table <- rbind(logistic_confusion, svm_linear_confusion, svm_radial_confusion)

  logistic_roc_df <- data.frame(false_positive_rate = 1 - logistic_roc$specificities,
    true_positive_rate = logistic_roc$sensitivities,
    model = "logistic")

  svm_linear_roc_df <- data.frame(false_positive_rate = 1 - svm_linear_roc$specificities,
    true_positive_rate = svm_linear_roc$sensitivities,
    model = "svm_linear")

  svm_radial_roc_df <- data.frame(false_positive_rate = 1 - svm_radial_roc$specificities,
    true_positive_rate = svm_radial_roc$sensitivities,
    model = "svm_radial")

  roc_table <- rbind(logistic_roc_df, svm_linear_roc_df, svm_radial_roc_df)

  auc_labels <- metrics_table
  auc_labels$label <- paste0(auc_labels$model, " AUC = ", round(auc_labels$AUC, 3))

  roc_plot <- ggplot(roc_table, aes(x = false_positive_rate, y = true_positive_rate, color = model)) +
    geom_line(linewidth = 1) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey55") +
    theme_minimal() +
    labs(title = "ROC Curves for Predictive Models",
      x = "False positive rate",
      y = "True positive rate",
      color = "Model") +
    annotate("text",
      x = 0.65,
      y = c(0.25, 0.18, 0.11),
      label = auc_labels$label,
      hjust = 0,
      size = 3.5)

  confusion_plot <- ggplot(confusion_table, aes(x = Reference, y = Prediction, fill = Freq)) +
    geom_tile(color = "white") +
    geom_text(aes(label = Freq), size = 5, fontface = "bold") +
    facet_wrap(~model) +
    scale_fill_gradient(low = "white", high = "#2c7fb8") +
    theme_minimal() +
    labs(title = "Confusion Matrices for Predictive Models",
      x = "Actual class",
      y = "Predicted class",
      fill = "Count")

  if (!is.null(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    write.csv(metrics_table, file.path(output_dir, "predictive_model_metrics.csv"), row.names = FALSE)
    write.csv(predictions, file.path(output_dir, "predictive_model_predictions.csv"), row.names = FALSE)
    write.csv(confusion_table, file.path(output_dir, "predictive_model_confusion_matrices.csv"), row.names = FALSE)
    write.csv(roc_table, file.path(output_dir, "predictive_model_roc_curve_data.csv"), row.names = FALSE)
    write.csv(data.frame(CpG = selected_cpgs), file.path(output_dir, "selected_cpgs.csv"), row.names = FALSE)
    ggsave(file.path(output_dir, "predictive_model_roc_curves.png"), roc_plot, width = 7, height = 5, dpi = 300)
    ggsave(file.path(output_dir, "predictive_model_roc_curves.pdf"), roc_plot, width = 7, height = 5)
    ggsave(file.path(output_dir, "predictive_model_confusion_matrices.png"), confusion_plot, width = 9, height = 4.5, dpi = 300)
    ggsave(file.path(output_dir, "predictive_model_confusion_matrices.pdf"), confusion_plot, width = 9, height = 4.5)
  }

  list(selected_cpgs = selected_cpgs,
    preprocess = preprocess_steps,
    logistic_model = logistic_model,
    svm_linear = svm_linear,
    svm_radial = svm_radial,
    metrics = metrics_table,
    predictions = predictions,
    confusion_table = confusion_table,
    roc_table = roc_table,
    roc_plot = roc_plot,
    confusion_plot = confusion_plot)}
