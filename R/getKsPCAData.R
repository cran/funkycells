#' Get K Functions and Compute Principal Components
#'
#' This function computes K functions from point process data then converts it
#'  into PCs. Note, if there are replicates, i.e. multiple observations per
#'  unit, the K functions will be a weighted average based on the number of
#'  the first agents.
#'
#' @inheritParams getKFunction
#' @param data Data.frame with column titles for at least outcome, x, y, agents,
#'     and unit. For consistency (and avoiding errors), use that order.
#'     Additionally, replicate can be added.
#' @param outcome (Optional) String of the column name in data indicating the
#'     outcome or response. Default is the 1st column.
#' @param unit (Optional) String of the column name in data indicating a unit or
#'     base thing. Note this unit may have replicates. Default is the 4th
#'     column.
#' @param replicate (Optional) String of the column name in data indicating the
#'  replicate id. Default is NULL.
#' @param rCheckVals (Optional) Numeric vector indicating the radius to check.
#'    Note, if not specified, this could take a lot of memory, particularly
#'    with many units and replicates.
#' @param nPCs (Optional) Numeric indicating the number of principal components.
#' @param agents_df (Optional) Two-column data.frame. The first for agent 1
#'     and the second for agent 2. Both should be in data agents column. This
#'     determines which K functions to compute. Default is to compute all, but
#'     may be misspecified if the data is in a different order.
#' @param xRange,yRange (Optional) Two value numeric vector indicating the min
#'  and max x/y values. Note this is re-used for all replicates. The default
#'  just takes the min and max x from each replicate. This allows different
#'  sized images, but the edges are defined by some agent location.
#' @param nbasis (Optional) Numeric indicating number of basis functions to fit
#'     K functions in order to compute PCA. Current implementation uses a
#'     b-spline basis.
#' @param silent (Optional) Boolean indicating if progress should be printed.
#' @param displayTVE (Optional) Boolean to  indicate if total variance explained
#'   (TVE) should be displayed. Default is FALSE.
#'
#' @return Data.frame with the outcome, unit and principle components of
#'     computed K functions.
#' @export
#'
#' @examples
#' dataPCA_pheno <- getKsPCAData(
#'   data = TNBC_pheno, unit = "Person",
#'   agents_df = data.frame(rep("B", 2), c("Tumor", "Fake")),
#'   nPCs = 3,
#'   rCheckVals = seq(0, 50, 1),
#'   displayTVE = TRUE
#' )
getKsPCAData <- function(data, outcome = colnames(data)[1],
                         unit = colnames(data)[5],
                         replicate = NULL,
                         rCheckVals = NULL, nPCs = 3,
                         agents_df = as.data.frame(expand.grid(
                           unique(data[, 4]),
                           unique(data[, 4])
                         )),
                         xRange = NULL, yRange = NULL,
                         edgeCorrection = "isotropic", nbasis = 21,
                         silent = FALSE, displayTVE = FALSE) {
  # Define pcaData
  pcaData <- unique(data[, c(outcome, unit)])

  ## Compute PCA for each agent-agent K-function
  if (!silent) cat("PCA Pairs (", nrow(agents_df), "): ", sep = "")
  pcaData_list <- apply(
    cbind(
      1:nrow(agents_df),
      as.data.frame(agents_df)
    ),
    MARGIN = 1,
    FUN = function(agents, nPCs, rCheckVals, data,
                   outcome, unit,
                   replicate,
                   xRange, yRange,
                   edgeCorrection, nbasis,
                   maxIters) {
      if (!silent) cat(trimws(agents[1]), sep = "")

      ## Compute K-Function for each unit
      #     Reminder, replicates -> weighted averages
      evaled_fd_K <- getKFunction(
        agents = agents[-1],
        unit = unit,
        replicate = replicate,
        data = data[, colnames(data) != outcome],
        rCheckVals = rCheckVals,
        xRange = xRange, yRange = yRange,
        edgeCorrection = edgeCorrection
      )

      ## Get PCA Scores
      K_pca_scores <- .getPCs(
        rKData = evaled_fd_K,
        agents = agents[-1], nPCs = nPCs,
        nbasis = nbasis, silent = !displayTVE
      )

      if (!silent && agents[1] < maxIters) cat(", ", sep = "")


      # Set up data (add counts as desired)
      retData <- cbind(
        "Unit" = unique(data[, unit]),
        as.data.frame(K_pca_scores)
      )
      colnames(retData) <- c(unit, colnames(K_pca_scores))
      retData
    },
    nPCs = nPCs, rCheckVals = rCheckVals,
    data = data, outcome = outcome, unit = unit,
    replicate = replicate,
    xRange = xRange, yRange = yRange,
    edgeCorrection = edgeCorrection,
    nbasis = nbasis, maxIters = nrow(agents_df)
  )
  if (!silent) cat("\n")

  # Potentially convert to do.call("rbind"), but unsure if order always same
  .mergeListsToDF(
    df = pcaData, lists = pcaData_list,
    dfCol = unit, listsDFCol = unit
  )
}


#' Get Principal Components
#'
#' This (internal) function converts the radius and K functions into principal
#'  components and returns the scores.
#'
#' See getKsPCAData for use.
#'
#' @inheritParams getKsPCAData
#' @param rKData data.frame with the first column being the checked radius and the
#'     rest relating to the K function for each unit at those points. NA columns
#'     for any K functions that could not be computed will be handled.
#' @param agents Two value vector indicating the two agents used for the K
#'     function, the first to the second.
#' @param nPCs Numeric indicating the number of principal components to compute
#' @param nbasis Numeric indicating the number of basis functions to use to fit
#'    the data to a bspline basis.
#'
#' @return Data.frame with the outcomes, units, then principal component scores.
#' @noRd
.getPCs <- function(rKData, agents, nPCs, nbasis = 21, silent = FALSE) {
  # Setup Data
  KData <- rKData[, -1]
  evalPts <- rKData[[1]]

  # Check if any are missing
  dropIdx <- which(colSums(is.na(KData)) != 0)
  if (length(dropIdx) > 0) {
    KData <- KData[, -dropIdx]
  }

  # Skip if not atleast 1 nonNA col
  if (ncol(KData) <= 1 || nrow(KData) == 0) {
    # No multiple people have K-function so no PCA
    K_pca_scores <- matrix(nrow = ncol(KData), ncol = 1)
    colnames(K_pca_scores) <- paste0(agents[1], "_", agents[2], "_PC")
    return(K_pca_scores)
  }

  K_func <- fda::Data2fd(
    argvals = evalPts,
    y = as.matrix(KData),
    basisobj =
      fda::create.bspline.basis(
        rangeval = range(evalPts),
        nbasis = nbasis
      )
  )
  K_pca <- fda::pca.fd(K_func, nharm = nPCs)
  if (!silent) cat(paste0(" (TVE: ", .specify_decimal(sum(K_pca$varprop), 3), ")"))

  if (length(dropIdx) > 0) {
    K_pca_scores <- .insertMissingRows(K_pca$scores, dropIdx)
  } else {
    K_pca_scores <- K_pca$scores
  }
  colnames(K_pca_scores) <- paste0(agents[1], "_", agents[2], "_PC", 1:nPCs)

  K_pca_scores
}


#' Insert Missing Rows Back into Data
#'
#' This (internal) function inserting missing observations back into PCA data.
#'     The data was dropped in order to allow PCA, but the NA need to be
#'     re-inserted in order to correctly analyze the data.
#'
#' See .getPCs for use.
#'
#' @param data_add A data.frame with the PCs.
#' @param insertRows Numeric vector indicating the original rows that were
#'     dropped.
#'
#' @return A data.frame with the PCs, now including the NA rows.
#' @noRd
.insertMissingRows <- function(data_add, insertRows) {
  data_return <- matrix(
    ncol = ncol(data_add), # nPCs
    nrow = (nrow(data_add) + length(insertRows))
  )
  currentAdd <- 1 # Current row I take data from

  for (i in 1:nrow(data_return)) {
    if (!(i %in% insertRows)) { # If this row shouldn't be NA
      data_return[i, ] <- as.numeric(data_add[currentAdd, ])
      currentAdd <- currentAdd + 1
    }
  }

  data_return
}
