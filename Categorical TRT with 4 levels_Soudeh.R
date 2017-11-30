setwd("C:/Users/soudi/Desktop/Multi Treatment Causal Modeling")

# ----------------------------------------------------------------------
# Get Libraries
# ----------------------------------------------------------------------
library(stremr)
library(data.table)
library(magrittr)
library(h2o)
options(stremr.verbose=TRUE)

# ----------------------------------------------------------------------
# Read Data 
# ----------------------------------------------------------------------
AD = readRDS('sampleAD.rds')

# ----------------------------------------------------------------------
# counterfactual regimen of interest are identifyied
# by indicator varaibles PATH1- PATH5 to show if pateint is on 
# dynamic regimen of interst or not
# ----------------------------------------------------------------------
AD[, ("TRT1.set")  := 1L]
AD[, ("TRT2.set")  := 2L]
AD[, ("TRT3.set")  := 3L]
AD[, ("TRT4.set")  := 4L]


# ----------------------------------------------------------------------
# Import Data
# ----------------------------------------------------------------------
OData.2  <-  importData(AD, ID = "ID", t_name = "SEQ", 
                        covars = c("CAT_VAR1","CAT_VAR2","CONT_VAR1"),           
                        CENS = c("CNS","ADM_CNS"), 
                        TRT = "TRTN",
                        MONITOR = NULL, OUTCOME = "STATUS",
                        weights = NULL, remove_extra_rows = TRUE,
                        verbose = getOption("stremr.verbose"))

# ----------------------------------------------------------------------
# Look at the input data object
# ----------------------------------------------------------------------
print(OData.2)

# ----------------------------------------------------------------------
# Access the input data
# ----------------------------------------------------------------------
get_data(OData.2)

# ----------------------------------------------------------------------
# Model the Right Censroing and Adminstrative Censoring
# ----------------------------------------------------------------------
gform_CENS <- "CNS + ADM_CNS ~ CAT_VAR1 + CONT_VAR1"

# ----------------------------------------------------------------------
# Estimate Propensity Scores
# fitPRopensity score with all defult option has error
# tried modeing treatmtnet with Gradient Boosting machines same error
# ----------------------------------------------------------------------
gform_TRT = "TRTN ~ CAT_VAR1 + CAT_VAR2 + CONT_VAR1"

OData.2 <- fitPropensity(OData.2, gform_CENS = gform_CENS,ngform_TRT = gform_TRT )

# ----------------------------------------------------------------------
# Fitting treatment model with Gradient Boosting machines:
# ----------------------------------------------------------------------
h2o::h2o.init(nthreads = -1)
models_TRT <- sl3::Lrnr_h2o_grid$new(algorithm = "gbm")
OData.2 <- fitPropensity(OData.2, gform_CENS = gform_CENS,
                         gform_TRT = gform_TRT,
                         models_TRT = models_TRT)


# Use `H2O-3` distributed implementation of GLM for treatment model estimator:
models_TRT <- sl3::Lrnr_h2o_glm$new(family = "multinomial")
OData.2 <- fitPropensity(OData.2, gform_CENS = gform_CENS,
                           gform_TRT = gform_TRT,
                           models_TRT = models_TRT)


# ----------------------------------------------------------------------
# Dynamic Treatment Patern of interset are defined as Dummy PATH1 - PATH5
# Always Treated with TRT1  then PATH1=1
# Always Treated with TRT2  then PATH2=1
# Always Treated with TRT3  then PATH3=1
# start with TRT2 switch to TRT1 PATH4=1
# start with TRT3 switch to TRT1 PATH5=1
#
#  how should the useonly_t_TRT & intervened_TRT be defined?
# ----------------------------------------------------------------------

# ----------------------------------------------------------------------
# Defining weights for two dynamic regimens
#how should the useonly_t_TRT & intervened_TRT be defined?
# with multiple Binary exposure and with categorical exposure
#  Defining weights for two dynamic regimens 
# ----------------------------------------------------------------------
wts.DT.1 <- getIPWeights(OData = OData.2)

wts.DT.1 <- getIPWeights(OData = OData.2, intervened_TRT="TRT1.set",
                         rule_name = "Always TRT1")

wts.DT.1 <- getIPWeights(OData = OData.2, intervened_TRT="TRT1.set",
                         useonly_t_TRT="PATH1==1", rule_name = "Always TRT1")

