#' Compute a Modified Random Forest Model
#'
#' This function creates a modified random forest model for principal component
#'  and meta-data. This can be useful to get a final model, but we recommend
#'  use of randomForest_CVPC in general, which includes the final model.
#'
#' @param data Data.frame of outcome and predictors. The predictors include
#'  groups of variables which are finite projections of a higher dimensional
#'  variables as well as single meta-variables.
#'
#'  Any replicate data, i.e. repeated observations, should already be handled.
#'    The unit column is needed just to drop data (so pre-removing and giving
#'    NULL works). Typically use the results from getKsPCAData, potentially with
#'    meta-variables attached.
#' @param outcome (Optional) String indicating the outcome column name in data.
#'   Default is the first column of data.
#' @param unit (Optional) String indicating the unit column name in data.
#'   Default is the second column of data.
#' @param nTrees (Optional) Numeric indicating the number of trees to use in the
#'     random forest model. Default is 500.
#' @param varImpPlot (Optional) Boolean indicating if variable importance plots
#'     should also be returned with the model. Default is TRUE.
#' @param metaNames (Optional) Vector with the column names of data that
#'     correspond to metavariables. Default is NULL.
#' @param keepModels (Optional) Boolean indicating if the individual models
#'     should be kept. Can get large in size. Default is TRUE as it is needed
#'     for predictions.
#' @param varSelPercent (Optional) Numeric in (0,1) indicating (approx)
#'  percentage of variables to keep for each tree. Default is 0.8.
#' @param method (Optional) Method for rpart tree to build random forest.
#'  Default is "class". Currently this is the only tested method. This will
#'  be expanded in future releases.
#'
#' @return A list with  entries
#'     \enumerate{
#'         \item varImportanceData: Data.frame for variable importance
#'                                  information.
#'         \item (Optional) model: List of CART that builds the random forest model.
#'         \item (Optional) varImportancePlot: Variable importance plots.
#'     }
#' @export
#'
#' @examples
#' ff <- funkyForest(
#'   data = TNBC[, c(1:8, ncol(TNBC))],
#'   outcome = "Class", unit = "Person",
#'   metaNames = c("Age")
#'   )
funkyForest <- function(data, outcome = colnames(data)[1],
                        unit = colnames(data)[2],
                        nTrees = 500, varImpPlot = TRUE, metaNames = NULL,
                        keepModels = TRUE, varSelPercent = 0.8,
                        method = "class") {
  # Ensure this is worthwhile
  if (length(unique(data[, outcome])) == 1) {
    stop("Error: Only 1 outcome in data, cannot do RF")
  }

  # Setup
  underlyingDataAlignedFunctions <- .getUnderlyingVariable(colnames(data),
    returnUnique = FALSE
  )
  underlyingVars <- unique(
    underlyingDataAlignedFunctions[!(underlyingDataAlignedFunctions %in%
      c(outcome, unit))]
  )

  data_result <- data.frame(
    "var" = underlyingVars,
    "VI" = 0
  )
  if (keepModels) RF <- list()
  # To do CART
  for (i in 1:nTrees) {
    # Bootstrap observations
    data_rf <- data[sample(1:nrow(data), nrow(data), replace = TRUE), ]
    # Ensure there are multiple responses in bootstrapped sample.
    while (length(unique(data_rf[, outcome])) == 1) {
      data_rf <- data[sample(1:nrow(data), nrow(data), replace = TRUE), ]
    }

    # Variable selection
    underlyingVarsSelected <- unique(sample(underlyingVars,
      varSelPercent * length(underlyingVars),
      replace = FALSE
    ))
    # Build RF data. Note estimators are scrambled to avoid computational bias
    data_rf_vars <-
      data_rf[underlyingDataAlignedFunctions %in% underlyingVarsSelected]
    data_rf <- cbind(data_rf[outcome], data_rf_vars[, sample(1:ncol(data_rf_vars))])

    # Fit CART
    model <- rpart::rpart(paste0(outcome, " ~ .-", outcome),
      data = data_rf, method = method,
      control = rpart::rpart.control(minsplit = 1, minbucket = 1, cp = 0)
    )

    # Get Var importance
    data_result <- .computeTotalVarImportance(model, data_result)

    if (keepModels) RF[[i]] <- model # Save model
  }

  # Get mean results for variable importance
  #   TODO:: See if this is even necessary
  data_result$avgVI <- data_result$VI / nTrees

  returnResults <- list("varImportanceData" = data_result)

  if (keepModels) {
    returnResults <- append(returnResults, list("model" = RF))
  }

  # Get Variable importance
  if (varImpPlot) {
    returnResults <- append(
      returnResults,
      list(
        "varImportancePlot" =
          .plotVariableImportance(varImportanceData = data_result)
      )
    )
  }

  returnResults
}


#' Compute Total Variable Importances
#'
#' This (internal) function computes the variable importance for a tree model
#'     and combines it with totVarImportance data.frame.
#'
#' See use in funkyForest.
#'
#' @param tree Model from rpart for fitting a CART tree.
#' @param totVarImportance Data.frame indicating the total (cumulative) variable
#'     importance.
#'
#' @return Data.frame which is the updated totVarImportance adding in the new
#'     tree information.
#' @noRd
.computeTotalVarImportance <- function(tree, totVarImportance) {
  tmp1 <- as.data.frame(tree$variable.importance)
  tmp2 <- data.frame(
    "var" = .getUnderlyingVariable(rownames(tmp1), returnUnique = FALSE),
    "VI" = tmp1[, 1]
  )

  # Combine Data
  # if (is.null(totVarImportance)) {
  #   totVarImportance <- data.frame("var" = unique(tmp2$var), "VI" = 0)
  # }

  for (v in unique(tmp2$var)) {
    totVarImportance[totVarImportance$var == v, "VI"] <-
      totVarImportance[totVarImportance$var == v, "VI"] + sum(tmp2[tmp2$var == v, "VI"])
  }

  totVarImportance
}


#' Plot Variable Importance (Ensure Data)
#'
#' This (internal) function organizes the data to plot the variable importance
#'     for a funkyForest model. One of the inputs is necessary.
#'
#' See use in funkyForest.
#'
#' Upcoming: See handle model without Var importance
#'
#' @param varImportanceData (Optional) Data.frame for the variable importance
#'     information.
#' @param model (Optional) Random forest model from funkyForest.
#'
#' @return grid.arrange containing two ggplots
#' @noRd
.plotVariableImportance <- function(varImportanceData = NULL, model = NULL) {
  if ((is.null(varImportanceData) && is.null(model)) ||
    (!is.null(varImportanceData) && !is.null(model))) {
    stop("Error: Give one varImportanceData or model")
  }

  if (is.null(varImportanceData)) {
    stop("Give Var Importance")
    # cols_preds never defined
    # varImportanceData <- data.frame('var'=cols_preds,
    #                                 # 'splits'=0, 'giniDec'=0,
    #                                 # 'varImp'=0, 'varImpCt'=0,
    #                                 'VI'=0)
    #
    # for(i in 1:length(model)){
    #   varImportanceData <- .computeTotalVarImportance(model[[i]], varImportanceData)
    # }

    # Get mean results for variable importance
    #     TEST: Use average over all or just when appears in model
    # varImportanceData$avgVI <- varImportanceData$VI / nTrees#data_result$varImpCt
  }

  .plotImportance(varImportanceData)
}


#' Plot Variable Importance (Plot Data)
#'
#' This (internal) function plots the variable importance for a
#' funkyForest model.
#'
#' See use in .plotVariableImportance.
#'
#' @param varImportanceData Data.frame for the variable importance information.
#'
#' @return grid.arrange containing ggplot
#' @noRd
.plotImportance <- function(varImportanceData) {
  # Add this to prevent NOTEs when building package
  var <- avgVI <- NULL

  varImportance <- ggplot2::ggplot() +
    ggplot2::geom_point(
      ggplot2::aes(
        x = stats::reorder(var, avgVI),
        y = avgVI / max(avgVI)
      ),
      data = varImportanceData
    ) +
    ggplot2::coord_flip() +
    ggplot2::xlab(NULL) +
    ggplot2::ylim(c(0, 1)) +
    ggplot2::ylab("VarImp") +
    ggplot2::theme_bw()

  varImportance
}


#' Predict a funkyForest
#'
#' This function gets the predicted value from a funkyForest model.
#'
#' @param model funkyForest model. See funkyForest. A list of
#'     CART models from rpart. Additionally this is given in funkyModel.
#' @param data_pred data.frame of the data to be predicted.
#' @param type (Optional) String indicating type of analysis. Options are pred
#'     or all. The choice changes the return to best fit intended use.
#' @param data (Optional) Data.frame of full data. The data used to fit the
#'     model will be extracted (by row name).
#'
#' @return The returned data depends on type:
#'     \itemize{
#'         \item type='pred': returns a vector of the predictions
#'         \item type='all': returns a vector of the predictions
#'     }
#' @export
#'
#' @examples
#' data_pp <- simulatePP(
#'   agentVarData =
#'     data.frame(
#'       "outcome" = c(0, 1),
#'       "A" = c(0, 0),
#'       "B" = c(1 / 50, 1 / 50)
#'     ),
#'   agentKappaData = data.frame(
#'     "agent" = c("A", "B"),
#'     "clusterAgent" = c(NA, "A"),
#'     "kappa" = c(10, 5)
#'   ),
#'   unitsPerOutcome = 5,
#'   replicatesPerUnit = 1,
#'   silent = FALSE
#' )
#' pcaData <- getKsPCAData(data_pp,
#'   replicate = "replicate",
#'   xRange = c(0, 1), yRange = c(0, 1), silent = FALSE
#' )
#' RF <- funkyForest(data = pcaData[-2], nTrees = 5) #
#' pred <- predict_funkyForest(
#'   model = RF$model, type = "all",
#'   data_pred = pcaData[-2],
#'   data = pcaData[-2]
#' )
predict_funkyForest <- function(model, data_pred, type = "all", data = NULL) {
  # Verification
  #   This checks the data to see if there are any columns that are characters
  #   or factors. This will be used later to manage missing data. We obviously
  #   can't predict a group we've never seen.
  if (!is.null(data)) {
    checkcols <- c()
    for (i in 1:ncol(data)) {
      checkcols <- c(
        checkcols,
        ifelse(is.character(data_pred[, i]) ||
          is.factor(data_pred[, i]), i, NA)
      )
    }
    checkcols <- checkcols[!is.na(checkcols)]
    if (length(checkcols) == 0) {
      data <- NULL
    } # No need for data if nothing to check
  }

  # Setup
  #   predictions: each row is an observation prediction and columns are models
  predictions <- data.frame(matrix(
    nrow = nrow(data_pred),
    ncol = length(model)
  ))
  data_pred_use <- data_pred # This will be used to manage missing data
  drop_pred_rows <- rep(0, nrow(data_pred)) # Rows of data_pred_use

  # Get prediction for each tree
  for (i in 1:length(model)) {
    # If original data given, this will allow verification
    #   Will ensure that all characters/factors are present in the fit model or
    #   set the data to missing. Missing data is managed by model, but not
    #   previously unseen data.
    if (!is.null(data)) {
      # Reset
      data_pred_use <- data_pred

      # Get rows of data used
      data_row <- names(unlist(model[[i]]$where))
      unique_rows <- unique(gsub("\\..*", "", data_row)) # Since bootstrap

      # Check each variable to see if contained in fit data
      #   Only check chars and factors (nums won't have this problem)
      drop_spot <- data.frame(matrix(0,
        ncol = length(checkcols),
        nrow = nrow(data_pred_use)
      )) # 0 keep, 1 drop
      for (j in checkcols) {
        # This flags any NA or data factors/characters the model hasn't seen
        #   If any are seen, set their values to NA. This will allow the model
        #   to use secondary splits.
        if (sum(!(data_pred_use[, j] %in% data[unique_rows, j])) > 0) {
          for (k in 1:nrow(data_pred_use)) {
            if (!(data_pred_use[k, j] %in% data[unique_rows, j])) {
              data_pred_use[k, j] <- NA
            }
          }
        }
      }
    }

    ## Predict the outcome using the model
    pred_mat <- stats::predict(model[[i]], data_pred_use)

    sols <- colnames(pred_mat)
    for (j in 1:nrow(data_pred)) {
      max_pred <- which.max(pred_mat[j, ])
      if (length(max_pred) != 0) {
        predictions[j, i] <- sols[max_pred]
      }
    }
  }


  # Combine predictions
  if (is.numeric(predictions)) {
    warning(paste(
      "Sorry: modelPercents incomplete. Only Preds will be returned"
    ))
    type <- "pred"
    # Take Mean
    modelPredictions <- rowMeans(predictions)
    # Perhaps return L2 error instead?
  } else {
    # Vote (Randomly selects on ties)
    modelPredictions <- apply(predictions,
      MARGIN = 1,
      function(x) {
        tmp <- .getMode(x)
        tmp[sample(1:length(tmp), 1)]
      }
    )
    # Percent correctly classified
    modelPercents <- as.data.frame(t(apply(predictions,
      MARGIN = 1,
      function(x, options) {
        sapply(options, function(option, x) {
          sum(x == option) / length(x)
        },
        x = x
        )
      }, options = unique(data_pred[, 1])
    )))
    colnames(modelPercents) <- unique(data_pred[, 1])
  }

  # Return Results
  if (type == "pred") {
    # This returns only the predictions
    return(modelPredictions)
  } else if (type == "all") {
    return(list(
      "PredPerc" = cbind(modelPredictions, modelPercents),
      "Acc" = sum(modelPredictions == data_pred[, 1]) / nrow(data_pred)
    ))
  }
}


#' Get Underlying Function
#'
#' This (internal) function is used to get the underlying function (i.e. no
#'     _PC\# in it). This is used to centralize called from funkyForest.R and
#'     RandomForest_CVPC.R. If there is no PC (i.e. in meta-variables), the
#'     meta-variable is returned. See use in funkyForest.
#'
#' @param names Vector of names to get underlying function on each
#'
#' @return Vector of underlying functions for each in name.
#' @noRd
.getUnderlyingVariable <- function(names, returnUnique = TRUE) {
  if (returnUnique) {
    return(unique(stringr::str_remove(names, "(_PC)[0-9]+$")))
  }
  return(stringr::str_remove(names, "(_PC)[0-9]+$"))
}
