## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(funkycells)

## ----define_functions, echo=FALSE---------------------------------------------
full_run <- FALSE # This runs only computations that are feasibly quick

## ----simulate_cells-----------------------------------------------------------
set.seed(123)
cell_data <- simulatePP(
  agentVarData =
    data.frame(
      "outcome" = c(0, 1),
      "A" = c(0, 0),
      "B" = c(1 / 100, 1 / 500),
      "C" = c(1 / 500, 1 / 250),
      "D" = c(1 / 100, 1 / 100)
    ),
  agentKappaData =
    data.frame(
      "agent" = c("A", "B", "C", "D"),
      "clusterAgent" = c(NA, "A", NA, NA),
      "kappa" = c(20, 5, 15, 15)
    ),
  unitsPerOutcome = 40,
  replicatesPerUnit = 1,
  silent = FALSE
)

## ----plot_data----------------------------------------------------------------
plotPP(cell_data[cell_data$replicate == 1, c("x", "y", "type")],
  dropAxes = TRUE,
  colorGuide = "none",
  xlim = c(0, 1),
  ylim = c(0, 1)
)

## ----compute_k----------------------------------------------------------------
k_data <- getKFunction(
  data = cell_data[cell_data$unit == "u1", -1],
  agents = c("A", "A"),
  unit = "unit",
  replicate = "replicate",
  rCheckVals = seq(0, 0.25, 0.01),
  xRange = c(0, 1),
  yRange = c(0, 1)
)

## ----compute_pca--------------------------------------------------------------
cells <- unique(cell_data$type)
cells_interactions <- rbind(
  data.frame(t(combn(cells, 2))),
  data.frame("X1" = cells, "X2" = cells)
)

pca_data <- getKsPCAData(
  data = cell_data,
  outcome = "outcome",
  unit = "unit",
  replicate = "replicate",
  rCheckVals = seq(0, 0.25, 0.01),
  agents_df = cells_interactions,
  xRange = c(0, 1), yRange = c(0, 1),
  nPCs = 3
)

## ----add_meta-----------------------------------------------------------------
set.seed(123)
pcaMeta <- simulateMeta(pca_data,
  outcome = "outcome",
  metaInfo = data.frame(
    "var" = c("gender", "age"),
    "rdist" = c("rbinom", "rnorm"),
    "outcome_0" = c("0.5", "30"),
    "outcome_1" = c("0.5", "31")
  )
)

## ----funkyModel, eval=full_run------------------------------------------------
#  set.seed(123)
#  model_fc <- funkyModel(
#    data = pcaMeta,
#    outcome = "outcome",
#    unit = "unit",
#    metaNames = c("gender", "age")
#  )

## ----variableimportance, eval=full_run----------------------------------------
#  model_fc$viPlot

## ----predict_data, eval=full_run----------------------------------------------
#  set.seed(12345)
#  cell_data_pred <- simulatePP(
#    agentVarData =
#      data.frame(
#        "outcome" = c(0, 1),
#        "A" = c(0, 0),
#        "B" = c(1 / 100, 1 / 500),
#        "C" = c(1 / 500, 1 / 250),
#        "D" = c(1 / 100, 1 / 100)
#      ),
#    agentKappaData =
#      data.frame(
#        "agent" = c("A", "B", "C", "D"),
#        "clusterAgent" = c(NA, "A", NA, NA),
#        "kappa" = c(20, 5, 15, 15)
#      ),
#    unitsPerOutcome = 2,
#    replicatesPerUnit = 2,
#    silent = FALSE
#  )
#  
#  pca_data_pred <- getKsPCAData(
#    data = cell_data_pred,
#    outcome = "outcome",
#    unit = "unit",
#    replicate = "replicate",
#    rCheckVals = seq(0, 0.25, 0.01),
#    agents_df = cells_interactions,
#    xRange = c(0, 1), yRange = c(0, 1),
#    nPCs = 3
#  )
#  
#  set.seed(12345)
#  pcaMeta_pred <- simulateMeta(pca_data_pred,
#    outcome = "outcome",
#    metaInfo = data.frame(
#      "var" = c("gender", "age"),
#      "rdist" = c("rbinom", "rnorm"),
#      "outcome_0" = c("0.5", "30"),
#      "outcome_1" = c("0.5", "31")
#    )
#  )
#  
#  predictions <- predict_funkyForest(model = model_fc$model, data_pred = pcaMeta_pred[-1])

## ----compare_k_functions------------------------------------------------------
ab_outcome0 <- getKFunction(cell_data[cell_data$outcome == 0, -1],
  agents = c("A", "B"), unit = "unit",
  replicate = "replicate",
  rCheckVals = seq(0, 0.25, 0.01)
)
ab_outcome1 <- getKFunction(cell_data[cell_data$outcome == 1, -1],
  agents = c("A", "B"), unit = "unit",
  replicate = "replicate",
  rCheckVals = seq(0, 0.25, 0.01)
)

ab_outcome0_long <- tidyr::pivot_longer(data = ab_outcome0, cols = K1:K15)
ab_outcome1_long <- tidyr::pivot_longer(data = ab_outcome1, cols = K1:K15)

data_k_plot <- rbind(
  data.frame(
    "r" = ab_outcome0_long$r,
    "K" = ab_outcome0_long$value,
    "unit" = ab_outcome0_long$name,
    "outcome" = "0"
  ),
  data.frame(
    "r" = ab_outcome1_long$r,
    "K" = ab_outcome1_long$value,
    "unit" = paste0(ab_outcome1_long$name, "_1"),
    "outcome" = "1"
  )
)

plot_K_functions(data_k_plot)

## ----compute_roc, eval=full_run-----------------------------------------------
#  pred_roc <- predict_funkyForest(
#    model = model_fc$model,
#    data_pred = pcaMeta[-2],
#    data = pcaMeta[-2]
#  )
#  computePseudoROCCurves(pcaMeta$outcome, pred_roc$PredPerc[-1])

