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
# Creating counterfactual node values to use in intervened_TRT
# how should we set up these nodes to get the Hazard for Trt regmines of
# interest? 
# example: 
#         Patient Stays only on TRT1 (PATH1==1):  
#               (TRT1==1 & TRT2==0 & TRT3==0 & TRT4==0) 
#         Patient Start with TRT2 then switch to TRT1, PATH3==1:
#              (TRT1==1 & TRT2==1 & TRT3==0 & TRT4==0) 
# ----------------------------------------------------------------------
AD[, ("zero.set") := 0L]
AD[, ("TRT1.set")  := 1L] 


# ----------------------------------------------------------------------
# Import Data: using 4 dummy variables for exposure
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
# Dynamic Treatment Pattern of interest are defined as Dummy PATH1 - PATH5
# Always Treated with TRT1  then PATH1=1
# Always Treated with TRT2  then PATH2=1
# Always Treated with TRT3  then PATH3=1
# start with TRT2 switch to TRT1 PATH4=1
# start with TRT3 switch to TRT1 PATH5=1
#
#  how should the useonly_t_TRT & intervened_TRT be defined?
#
# ----------------------------------------------------------------------

# ----------------------------------------------------------------------
# IPW-MSM for hazard :how should the useonly_t_TRT & intervened_TRT be defined?
# with multiple Binary exposure so we dont get an error about 
# ----------------------------------------------------------------------
wts.DT.1 <- getIPWeights(OData = OData.1)

#Error: length(intervened_NODE) not equal to length(NodeNames)
wts.DT.1 <- getIPWeights(OData = OData.1, intervened_TRT="TRT1.set",
                         useonly_t_TRT="PATH1==1", rule_name = "Always TRT1")


#Error in modelfit.g$getPsAsW.models()[[i]] : subscript out of bounds
wts.DT.1 <- getIPWeights(OData = OData.1, 
                         intervened_TRT=c("TRT1.set","zero.set","zero.set","zero.set"),
                         useonly_t_TRT="PATH1==1", rule_name = "Always TRT1")

# result dosent change with PATH1==1 or PATH2==1
wts.DT.1 <- getIPWeights(OData = OData.1, 
                         useonly_t_TRT="PATH1==1", rule_name = "Always TRT1")
wts.DT.1

wts.DT.2 <- getIPWeights(OData = OData.1, 
                         useonly_t_TRT="PATH2==1", rule_name = "Always TRT2")
wts.DT.2




