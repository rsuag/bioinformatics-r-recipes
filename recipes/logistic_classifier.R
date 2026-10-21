# Task: fit and evaluate a binary logistic classifier on a held-out test set.
# Input: data.frame with a 0/1 outcome and finite numeric predictor columns.
# Output: model, disjoint train/test rows, predictions, and test-set metrics.
# For real studies, split by patient/site/batch where observations are not independent.

evaluate_logistic <- function(data, outcome, features, train_fraction = 0.7, seed = 1L) {
  if (!is.data.frame(data) || !is.character(outcome) || length(outcome) != 1L ||
      !outcome %in% names(data) || !is.character(features) || !length(features) ||
      anyDuplicated(features) || any(!features %in% names(data)) || outcome %in% features) {
    stop("supply a data.frame, one outcome column, and distinct predictor columns")
  }
  y <- data[[outcome]]
  if (!is.numeric(y) || anyNA(y) || !all(y %in% c(0, 1)) ||
      any(tabulate(y + 1L, nbins = 2L) < 2L) ||
      any(!vapply(data[features], is.numeric, logical(1))) ||
      any(!vapply(data[features], function(x) all(is.finite(x)), logical(1)))) {
    stop("outcome must be 0/1 with at least two rows per class; predictors must be finite numeric")
  }
  if (length(train_fraction) != 1L || !is.numeric(train_fraction) || is.na(train_fraction) ||
      train_fraction <= 0 || train_fraction >= 1 || length(seed) != 1L ||
      !is.numeric(seed) || is.na(seed) || seed != floor(seed)) {
    stop("train_fraction must be in (0, 1) and seed must be one integer")
  }
  set.seed(seed)
  train <- unlist(lapply(0:1, function(class) {
    rows <- which(y == class)
    sample(rows, max(1L, min(length(rows) - 1L, floor(length(rows) * train_fraction))))
  }), use.names = FALSE)
  test <- setdiff(seq_len(nrow(data)), train)
  model <- glm(reformulate(features, response = outcome), data = data[train, , drop = FALSE],
               family = binomial())
  probability <- unname(predict(model, newdata = data[test, , drop = FALSE], type = "response"))
  predicted <- as.integer(probability >= 0.5)
  actual <- y[test]
  confusion <- table(actual = factor(actual, levels = 0:1),
                     predicted = factor(predicted, levels = 0:1))
  safe_ratio <- function(numerator, denominator) if (denominator) numerator / denominator else NA_real_
  metrics <- c(accuracy = mean(predicted == actual),
               sensitivity = safe_ratio(confusion["1", "1"], sum(confusion["1", ])),
               specificity = safe_ratio(confusion["0", "0"], sum(confusion["0", ])))
  list(model = model, train_rows = train, test_rows = test,
       predictions = data.frame(row = test, actual = actual, probability = probability,
                                predicted = predicted), confusion = confusion, metrics = metrics)
}

if (sys.nframe() == 0L) {
  set.seed(24)
  demo <- data.frame(signal_a = rnorm(60), signal_b = rnorm(60))
  demo$case <- rbinom(60, 1, plogis(-0.2 + demo$signal_a - 0.5 * demo$signal_b))
  result <- evaluate_logistic(demo, "case", c("signal_a", "signal_b"), seed = 24)
  print(result$confusion)
  print(round(result$metrics, 3))
}
