constrain_weight <- function(w, w_LB = -0.5, w_UB = 1.5) {

  pmin(pmax(w, w_LB), w_UB)

}
