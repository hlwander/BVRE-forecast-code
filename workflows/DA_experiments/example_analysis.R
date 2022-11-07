score_dir <- arrow::SubTreeFileSystem$create(file.path(lake_directory,"scores/da_study"))
scores <- arrow::open_dataset(score_dir) |> collect()


rmse <- scores |> 
  filter(!is.na(observation),
         variable == "temperature",
         horizon > 0) |> 
  mutate(sq_error = (mean - observation)^2) |> 
  group_by(horizon, model_id) |> 
  summarize(rmse = sqrt(mean(sq_error, na.rm = TRUE)), .groups = "drop")


ggplot(rmse, aes(x = horizon, y = rmse, color = model_id)) +
  geom_line()


scores |> 
  filter(variable == "lw_factor",
         horizon <= 0) |>
  ggplot(aes(x = datetime, y = mean, color = model_id)) +
  geom_line()