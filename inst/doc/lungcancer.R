## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(funkycells)

## ----define_vars, echo=FALSE--------------------------------------------------
full_run <- FALSE # This runs only computations that are feasibly quick

## ----load_data, eval=FALSE----------------------------------------------------
#  dat <- lungcancer

## ----define_functions, eval=FALSE---------------------------------------------
#  getInteractions <- function(cells,
#                              dropCells = c(
#                                "Epithelialcells",
#                                "Stromalcells", "Other"
#                              )) {
#    cells <- cells[!(cells %in% dropCells)]
#    return_df <- data.frame("c1" = cells[1], "c2" = cells[1])
#    if (length(cells) > 1) {
#      return_df <- rbind(
#        data.frame(t(combn(cells, 2))),
#        data.frame("X1" = cells, "X2" = cells)
#      )
#    }
#  
#    return_df
#  }
#  
#  clean_figures <- function(model_rf) {
#    specify_decimal <- .specify_decimal
#    # Get Vars
#    viData <- model_rf$VariableImportance
#    accData <- model_rf$AccuracyEstimate
#    NoiseCutoff <- model_rf$NoiseCutoff
#    InterpolationCutoff <- model_rf$InterpolationCutoff
#    subsetPlotSize <- model_rf$AdditionalParams$subsetPlotSize
#  
#    # Plot Full
#    maxVal <- max(InterpolationCutoff, NoiseCutoff, viData$est)
#    plot_full <- ggplot2::ggplot(
#      data = viData,
#      mapping = ggplot2::aes(
#        x = factor(stats::reorder(var, est)),
#        y = ifelse(est / maxVal > 1, 1,
#          ifelse(est / maxVal < 0, 0,
#            est / maxVal
#          )
#        ),
#        group = 1
#      )
#    ) +
#      ggplot2::geom_errorbar(
#        ggplot2::aes(
#          ymin = ifelse((est - sd) / maxVal < 0, 0, (est - sd) / maxVal),
#          ymax = ifelse((est + sd) / maxVal > 1, 1, (est + sd) / maxVal)
#        ),
#        color = "black", width = 0.2
#      ) +
#      ggplot2::geom_point(color = "black", size = 3) +
#      ggplot2::geom_hline(
#        ggplot2::aes(yintercept = max(0, min(1, NoiseCutoff / maxVal))),
#        color = "red", linetype = "dotted", linewidth = 2
#      ) +
#      ggplot2::coord_flip(ylim = c(0, 1)) +
#      ggplot2::xlab(NULL) +
#      ggplot2::ylim(c(0, 1)) +
#      ggplot2::ylab(NULL) +
#      ggplot2::theme_bw() +
#      ggplot2::theme(
#        axis.text = ggplot2::element_text(size = 18),
#        axis.text.y = ggplot2::element_blank(),
#        axis.ticks.y = ggplot2::element_blank(),
#        panel.grid.major.y = ggplot2::element_blank(),
#        axis.title = ggplot2::element_text(size = 20)
#      ) +
#      ggplot2::geom_line(ggplot2::aes(
#        x = ordered(viData[order(-est), "var"]),
#        y = InterpolationCutoff / maxVal
#      ), color = "orange", linetype = "dashed", linewidth = 2) +
#      ggplot2::ylab(paste0(
#        "OOB (",
#        specify_decimal(accData$OOB, 2),
#        "), Guess (",
#        specify_decimal(accData$guess, 2),
#        "), Bias (",
#        specify_decimal(accData$bias, 2), ")"
#      ))
#  
#    # Plot Subset
#    InterpolationCutoff1 <- InterpolationCutoff[1:subsetPlotSize]
#    viData1 <- viData[order(-viData$est), ]
#    viData1 <- viData1[1:subsetPlotSize, ]
#  
#    maxVal <- max(InterpolationCutoff1, NoiseCutoff, viData1$est)
#    plot_top25 <- ggplot2::ggplot(
#      data = viData1,
#      mapping = ggplot2::aes(
#        x = factor(stats::reorder(var, est)),
#        y = ifelse(est / maxVal > 1, 1,
#          ifelse(est / maxVal < 0, 0,
#            est / maxVal
#          )
#        ),
#        group = 1
#      )
#    ) +
#      ggplot2::geom_errorbar(
#        ggplot2::aes(
#          ymin = ifelse((est - sd) / maxVal < 0, 0, (est - sd) / maxVal),
#          ymax = ifelse((est + sd) / maxVal > 1, 1, (est + sd) / maxVal)
#        ),
#        color = "black", width = 0.2
#      ) +
#      ggplot2::geom_point(color = "black", size = 3) +
#      ggplot2::geom_hline(
#        ggplot2::aes(yintercept = max(0, min(1, NoiseCutoff / maxVal))),
#        color = "red", linetype = "dotted", linewidth = 2
#      ) +
#      ggplot2::coord_flip(ylim = c(0, 1)) +
#      ggplot2::xlab(NULL) +
#      ggplot2::ylim(c(0, 1)) +
#      ggplot2::ylab(NULL) +
#      ggplot2::theme_bw() +
#      ggplot2::theme(
#        axis.text.x = ggplot2::element_text(size = 17),
#        axis.text.y = ggplot2::element_text(size = 14),
#        axis.ticks.y = ggplot2::element_blank(),
#        panel.grid.major.y = ggplot2::element_blank(),
#        axis.title = ggplot2::element_text(size = 17)
#      ) +
#      ggplot2::geom_line(ggplot2::aes(
#        x = ordered(viData1[order(-est), "var"]),
#        y = InterpolationCutoff1 / maxVal
#      ), color = "orange", linetype = "dashed", linewidth = 2) +
#      ggplot2::ylab(paste0(
#        "OOB (",
#        specify_decimal(accData$OOB, 2),
#        "), Guess (",
#        specify_decimal(accData$guess, 2),
#        "), Bias (",
#        specify_decimal(accData$bias, 2), ")"
#      ))
#  
#    list(
#      "full" = plot_full,
#      "top" = plot_top25
#    )
#  }

## ----pca_data, eval=full_run--------------------------------------------------
#  pcaData <- getKsPCAData(
#    data = use_data[!(colnames(use_data) %in% c("SmokingStatusLong"))],
#    outcome = "Cancer",
#    unit = "Patient",
#    replicate = "Image",
#    rCheckVals = seq(0, 50, 1),
#    nPCs = 3,
#    agents_df = ints,
#    displayTVE = FALSE
#  )
#  
#  data_rf <- merge(pcaData,
#    unique(use_data[, c("Patient", "SmokingStatusLong")]),
#    keep.x = TRUE, by.x = "Patient", by.y = "Patient"
#  )

## ----fit_model, eval=full_run-------------------------------------------------
#  set.seed(12345)
#  model_fm <- funkyModel(
#    data = data_rf[!is.na(data_rf[["Cancer"]]), ],
#    outcome = "Cancer", unit = "Patient",
#    metaNames = "SmokingStatusLong"
#  )

## ----plot_data, eval=full_run-------------------------------------------------
#  figs <- clean_figures(model_fm)
#  figs$full
#  fits$top

