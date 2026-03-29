# DEPRECATED: R implementation moved to C++ in src/resampling_entropy_rcpp.cpp
# This function is no longer used. All jackknife operations now use the
# optimized C++/Rcpp implementation for better performance.
#
# The C++ implementation is called via .jackknife_resampling_hybrid()
# See R/jackknife_rcpp_wrapper.R for integration details.

# DEPRECATED: Batch jackknife also moved to C++ implementation
# Use .jackknife_resampling_hybrid() for all jackknife operations
