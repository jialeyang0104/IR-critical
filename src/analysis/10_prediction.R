#!/usr/bin/env Rscript
root <- Sys.getenv("IR_REVISION_ROOT")
if (!nzchar(root)) stop("Set IR_REVISION_ROOT")
source(file.path(root,"src/shared/04_prediction.R"))
