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
# Creating Treatment Nodes to use in intervened_TRT
# ----------------------------------------------------------------------
AD[, ("zero.set") := 0L]
AD[, ("TRT1.set")  := 1L] 


# ----------------------------------------------------------------------
# Import Data
# ----------------------------------------------------------------------
OData.1  <-  importData(AD, ID = "ID", t_name = "SEQ", 
                        covars = c("CAT_VAR1","CAT_VAR2","CONT_VAR1"),           
                        CENS = c("CNS","ADM_CNS"), 
                        TRT = c("TRT1","TRT2","TRT3","TRT4"),
                        MONITOR = NULL, OUTCOME = "STATUS",
                        weights = NULL, remove_extra_rows = TRUE,
                        verbose = getOption("stremr.verbose"))

# ----------------------------------------------------------------------
# Look at the input data object
# ----------------------------------------------------------------------
print(OData.1)

# ----------------------------------------------------------------------
# Access the input data
# ----------------------------------------------------------------------
get_data(OData.1)

# ----------------------------------------------------------------------
# Model the Right Censroing and Adminstrative Censoring
# ----------------------------------------------------------------------
gform_CENS <- "CNS + ADM_CNS ~ CAT_VAR1 + CONT_VAR1"

# ----------------------------------------------------------------------
# Fit Propensity Scores
# ----------------------------------------------------------------------
gform_TRT = "TRT1+TRT2+TRT3+TRT4 ~ CAT_VAR1 + CAT_VAR2 + CONT_VAR1"
OData.1 <- fitPropensity(OData.1, gform_CENS = gform_CENS,ngform_TRT = gform_TRT )

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
# IPW-MSM for hazard :how should the useonly_t_TRT & intervened_TRT be defined?
# with multiple Binary exposure and with categorical exposure
# ----------------------------------------------------------------------
wts.DT.1 <- getIPWeights(OData = OData.1)

wts.DT.1 <- getIPWeights(OData = OData.1, intervened_TRT="TRT1.set",
                         useonly_t_TRT="PATH1==1", rule_name = "Always TRT1")

wts.DT.1 <- getIPWeights(OData = OData.1, 
                         useonly_t_TRT="PATH1==1", rule_name = "Always TRT1")
wts.DT.1

wts.DT.2 <- getIPWeights(OData = OData.1, 
                         useonly_t_TRT="PATH2==1", rule_name = "Always TRT2")
wts.DT.2



# ----------------------------------------------------------------------
# Fitting treatment model with Gradient Boosting machines:
# ----------------------------------------------------------------------
h2o::h2o.init(nthreads = -1)
models_TRT <- sl3::Lrnr_h2o_grid$new(algorithm = "gbm")
OData.1 <- fitPropensity(OData.1, gform_CENS = gform_CENS,
                       gform_TRT = gform_TRT,
                       models_TRT = models_TRT)


# Use `H2O-3` distributed implementation of GLM for treatment model estimator:
models_TRT <- sl3::Lrnr_h2o_glm$new(family = "binomial")
OData.1 <- fitPropensity(OData.1, gform_CENS = gform_CENS,
                       gform_TRT = gform_TRT,
                       models_TRT = models_TRT)

# Use Deep Neural Nets:
models_TRT <- sl3::Lrnr_h2o_grid$new(algorithm = "deeplearning", family = "AUTO")
OData.1 <- fitPropensity(OData.1, gform_CENS = gform_CENS,
                       gform_TRT = gform_TRT,
                       models_TRT = models_TRT)

# ----------------------------------------------------------------------
# Fitting different models with different algorithms
# Fine tuning modeling with optional tuning parameters.
# ----------------------------------------------------------------------
## Not run: 
params_TRT <- sl3::Lrnr_h2o_grid$new(algorithm = "gbm",
                                     ntrees = 50,
                                     learn_rate = 0.05,
                                     sample_rate = 0.8,
                                     col_sample_rate = 0.8,
                                     balance_classes = TRUE)
params_CENS <- sl3::Lrnr_glm_fast$new()

OData.1 <- fitPropensity(OData.1,
                       gform_CENS = gform_CENS,  params_CENS = params_CENS,
                       gform_TRT = gform_TRT, params_TRT = params_TRT)
