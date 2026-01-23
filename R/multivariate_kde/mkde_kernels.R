gaussian_kernel <- function (normalized_observation){
  1/sqrt(2*pi)*exp((-normalized_observation^2)/2)
}