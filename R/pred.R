validatePredInputs <- function( X, y, vTest )
{
    ## Argument verification
    stopifnot( is.matrix(X) )
    stopifnot( all(names(y) == rownames(X)) )
    stopifnot( all(vTest %in% rownames(X)) )
    stopifnot( length(setdiff(XY$Label, c("neg","pos", as.character(seq(0, 6))))) == 0 )
}

## Train-test for a single pair using logistic regression
## X - matrix of expression values; rownames are sample IDs; colnames are genes
## y - character vector of labels, sampled from {"neg","pos"}
## vTest - IDs to withhold for testing
## Returns a length(vTest)-by-3 data frame containing test IDs, true Labels and predictions
liblinear_lgr <- function( X, y, vTest )
{
    validatePredInputs( X, y, vTest )

    ## Split the data into train and test
    vTrain <- setdiff( rownames(X), vTest )
    Xte <- X[vTest,]
    Xtr <- X[vTrain,]
    ytr <- y[vTrain]

    ## Train a model and apply it to test data
    m <- LiblineaR::LiblineaR( Xtr, ytr, type=0 )
    ypred <- predict( m, Xte, proba=TRUE )$probabilities[,"pos"]
    tibble::enframe( y[vTest], "ID", "Label" ) %>% dplyr::mutate( Pred = ypred )
}

## Train-test for a single pair using support vector machines
## X - matrix of expression values; rownames are sample IDs; colnames are genes
## y - character vector of labels, sampled from {"neg","pos"}
## vTest - IDs to withhold for testing
## Returns a length(vTest)-by-3 data frame containing test IDs, true Labels and predictions
liblinear_svm <- function( X, y, vTest )
{
    validatePredInputs( X, y, vTest )

    ## Split the data into train and test
    vTrain <- setdiff( rownames(X), vTest )
    Xte <- X[vTest,]
    Xtr <- X[vTrain,]
    ytr <- y[vTrain]

    ## Train a model and apply it to test data
    m <- LiblineaR::LiblineaR( Xtr, ytr, type=2 )
    p <- predict( m, Xte, decisionValues=TRUE )
    ypred <- `if`( !identical(p$decisionValues[,"pos"], c(0,0)),
                  p$decisionValues[,"pos"],
                  -p$decisionValues[,"neg"] )
    tibble::enframe( y[vTest], "ID", "Label" ) %>% dplyr::mutate( Pred = ypred )
}

xgboost <- function( X, y, vTest )
{
    validatePredInputs( X, y, vTest )

    ## Convert response to 0,1
    y01 <- ifelse( y == "pos", 1, 0 )

    ## Split the data into train and test
    vTrain <- setdiff( rownames(X), vTest )
    Xte <- X[vTest,]
    Xtr <- X[vTrain,]
    ytr <- y01[vTrain]

    mdl <- xgboost::xgboost( Xtr, ytr, nrounds=20, verbose=0 )

    ## Train a model and apply it to test data
    ypred <- predict( mdl, Xte )
    tibble::enframe( y[vTest], "ID", "Label" ) %>% dplyr::mutate( Pred = ypred )
}

nnet <- function( X, y, vTest )
{
    validatePredInputs( X, y, vTest )

    ## Convert response to 0,1
    y01 <- ifelse( y == "pos", 1, 0 )

    ## Split the data into train and test
    vTrain <- setdiff( rownames(X), vTest )
    Xte <- X[vTest,]
    Xtr <- X[vTrain,]
    ytr <- y01[vTrain]

    mdl <- purrr::quietly(nnet::nnet)( Xtr, ytr, size=2, MaxNWts=10000 )$result

    ## Train a model and apply it to test data
    ypred <- predict( mdl, Xte ) %>% as.data.frame %>%
        tibble::rownames_to_column("ID") %>% dplyr::rename( Pred = V1 )
    tibble::enframe( y[vTest], "ID", "Label" ) %>% dplyr::inner_join(ypred, by="ID")
}

oforest <- function( X, y, vTest )
{
    validatePredInputs( X, y, vTest )

    ## Convert response to ordered factor
    y01 <- factor( y, levels = as.character(seq(0, 6)), ordered = TRUE )

    ## Split the data into train and test
    vTrain <- setdiff( rownames(X), vTest )
    dfs <- X %>%
        as.data.frame() %>%
        dplyr::mutate( Label = y01 ) %>%
        split( if_else(rownames(.) %in% vTest, "test", "train") )

    mdl <- ordinalForest::ordfor( depvar = "Label", data = dfs$train )

    preds <- predict(mdl, newdata = dfs$test)

    dfs$test %>%
        tibble::as_tibble( rownames = "ID" ) %>%
        transmute( ID, Label, Pred = preds$ypred )
}

onet <- function( X, y, vTest )
{
    validatePredInputs( X, y, vTest )

    ## Convert response to ordered factor
    y01 <- factor( y, levels = as.character(seq(0, 6)), ordered = TRUE )

    XX <- X %>% t() %>% cov()

    ## Split the data into train and test
    vTrain <- setdiff( rownames(XX), vTest )
    Xte <- XX[vTest, vTrain]
    Xtr <- XX[vTrain, vTrain]
    ytr <- y01[vTrain]

    ## Train a model and apply it to test data
    mdl <- ordinalNet::ordinalNet(
        Xtr, ytr, alpha = 0, threshIn = 1e-5, threshOut = 1e-5,
        keepTrainingData = FALSE
    )
    preds <- predict(mdl, Xte, type = "response")

    # Encoding order of pair so that correct ordering is c(1, 2)
    # and incorrect is c(2, 1). Ties are c(1, 1)
    ytrue <- match(y01[vTest], sort(y01[vTest]))
    ypred <- match(preds[, "P[Y=7]"], sort(preds[, "P[Y=7]"]))

    tibble(
        ID = vTest, Label = ytrue, Pred = ypred
    )
}
