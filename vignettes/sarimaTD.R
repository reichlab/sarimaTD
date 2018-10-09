## ------------------------------------------------------------------------
library(sarimaTD)
library(dplyr)
library(ggplot2)

ggplot(data = san_juan_dengue,
    mapping = aes(x = week_start_date, y = total_cases)) +
  geom_line()

## ------------------------------------------------------------------------
train_indices <- san_juan_dengue$season %in% paste0(1990:2004, "/", 1991:2005)

sarima_fit_bc_transform <- fit_sarima(
  y = san_juan_dengue$total_cases[train_indices],
  ts_frequency = 52,
  transformation = "box-cox",
  seasonal_difference = TRUE)

#sarima_fit_log_transform <- fit_sarima(
#  y = san_juan_dengue$total_cases[train_indices] + 1, # there are some observations of 0 cases
#  ts_frequency = 52,
#  transformation = "log",
#  seasonal_difference = TRUE)

#sarima_fit_no_transform <- fit_sarima(
#  y = san_juan_dengue$total_cases[train_indices],
#  ts_frequency = 52,
#  transformation = "none",
#  seasonal_difference = TRUE)

## ------------------------------------------------------------------------
eval_indices <- san_juan_dengue$season %in% paste0(1990:2005, "/", 1991:2006)

sampled_trajectories_bc_transform <-
  simulate(
    object = sarima_fit_bc_transform,
    nsim = 1000,
    seed = 1,
    newdata = san_juan_dengue$total_cases[eval_indices],
    h = 52
  )

plot_indices <- san_juan_dengue$season %in% paste0(2003:2006, "/", 2004:2007)

preds_df <- san_juan_dengue %>%
  filter(season == "2006/2007") %>%
  mutate(
    pred_total_cases = apply(sampled_trajectories_bc_transform, 2, mean),
    pred_95_lb = apply(sampled_trajectories_bc_transform, 2, quantile, probs = 0.025),
    pred_95_ub = apply(sampled_trajectories_bc_transform, 2, quantile, probs = 0.975),
    pred_80_lb = apply(sampled_trajectories_bc_transform, 2, quantile, probs = 0.05),
    pred_80_ub = apply(sampled_trajectories_bc_transform, 2, quantile, probs = 0.95),
    pred_50_lb = apply(sampled_trajectories_bc_transform, 2, quantile, probs = 0.25),
    pred_50_ub = apply(sampled_trajectories_bc_transform, 2, quantile, probs = 0.75),
  )

ggplot() +
  geom_ribbon(
    mapping = aes(x = week_start_date, ymin = pred_95_lb, ymax = pred_95_ub),
    fill = "cornflowerblue",
    alpha = 0.2,
    data = preds_df) +
  geom_ribbon(
    mapping = aes(x = week_start_date, ymin = pred_80_lb, ymax = pred_80_ub),
    fill = "cornflowerblue",
    alpha = 0.2,
    data = preds_df) +
  geom_ribbon(
    mapping = aes(x = week_start_date, ymin = pred_50_lb, ymax = pred_50_ub),
    fill = "cornflowerblue",
    alpha = 0.2,
    data = preds_df) +
  geom_line(
    mapping = aes(x = week_start_date, y = pred_total_cases),
    color = "cornflowerblue",
    data = preds_df) +
  geom_line(mapping = aes(x = week_start_date, y = total_cases),
    data = san_juan_dengue[plot_indices, , drop = FALSE])

