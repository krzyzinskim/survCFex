#saveRDS(experiment2_results, "experiments/results/experiment2_results.rds")


experiment2_results <- readRDS("experiments/results/experiment2_results.rds")
experiment2_results$moc_results[[1]]$target_envelope
analyze(experiment2_results$moc_results[[1]])


sapply(experiment2_results$moc_results, function(x) max(x$history$generation))
