# percentile-residuals

Main code: percentile_residuals.r
Example usage for the gamma distribution is given in the code gamma.R

The percentile_residuals.r currently supports normal (redundant, but just there for sanity checks), gamma (covers exponential as a special case) and mixture_normal (no longer important as of now)

The normal one can be easily modified to get the log_normal (using the exp transformation)

Test options include both, left and right sided alternatives, and range of the type I error over which we will plot the ROC curve. 
