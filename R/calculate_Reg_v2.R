calculate_Reg <- function(formula, df, target) {
  model <- lm(formula, data = df)
  preds <- predict(model, df)
  
  rmse <- sqrt(mean((preds - df[[target]])^2, na.rm = TRUE))
  adj_r2 <- summary(model)$adj.r.squared
  aic <- AIC(model)
  n <- nobs(model)
  
  return(c(
    RMSE = rmse,
    AdjR2 = adj_r2,
    AIC = aic,
    n = n
  ))
}
