## Data directory
if (!dir.exists("../data"))
    dir.create(path = "../data")

## Download dataset
if (!file.exists("../data/Griner_et_al.zip")) {
    utils::download.file(
        url = "https://ars.els-cdn.com/content/image/1-s2.0-S2001037015000422-mmc4.zip",
        method = "wget",
        destfile = "../data/Griner_et_al.zip"
    )
}

if (!all(file.exists(
    c("../data/metadata.csv",
      "../data/responses.csv")
    ))) {
    utils::unzip(
        zipfile = "../data/Griner_et_al.zip",
        files = c("metadata.csv", "responses.csv"), ## extract dataset files only
        exdir = "../data"
    )
}

## install required packages
required_pkgs <- c(
    "drc",
    "caTools",
    "kriging",
    "lattice",
    "reshape2",
    "plotrix",
    "compiler",
    "truncnorm",
    "data.table"
)
if (!require(required_pkgs, quietly = TRUE))
    install.packages(setdiff(required_pkgs, rownames(installed.packages())))

## original paper result
source(file.path(getwd(), "Delta_calculation_original.R"))

## result after changing logistic model to log-logistic
## cannot define parameter boundaries due to many negative responses
source(file.path(getwd(), "Delta_calculation.R"))

library(data.table)

delta_score_published <- fread("../data/Delta_score_original.csv")
delta_score_new <- fread("../data/Delta_score.csv")

## Experiments with different delta scores
delta_score_diff <- delta[
    which(
        sign(delta_score_published$Delta) !=
        sign(delta_score_new$Delta)
    )
] 
