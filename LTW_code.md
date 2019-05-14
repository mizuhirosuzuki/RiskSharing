Replication of Figure 1 in Ligon, Thomas & Worrall (2002)
================
Mizuhiro Suzuki
5/9/2018

Replicate Figure 1 of Ligon, Thomas & Worrall (2002)
----------------------------------------------------

In this document, my goal is to replicate Figure 1 of Ligon, Thomas & Worrall (2002). For this, I need to derive upper and lower bounds of intervals of the ratio of marginal utilities. I derive them by considering 4 cases separately; autarky case, non-overlapping case, overlapping case, and first best case.

Before diving into the derivations, there are several facts that we know (p.221 in the paper):

1.  since all penalties are zero, from Proposition 2(iv), *λ*<sub>*l**h*, *l**o**w**e**r*</sub> = *ξ*<sub>*l**h*</sub> = 1/2 and *λ*<sub>*h**l*, *u**p**p**e**r*</sub> = *ξ*<sub>*h**l*</sub> = 2, and
2.  because of the log-form utility functions, the *h**h* and *l**l* *λ*-intervals are identical. <!--3. since preferences are identical, by symmetry $\lambda_{lh, lower} = 1 / \lambda_{hl, upper}$ and $\lambda_{hh, lower} = 1 / \lambda_{hh, upper}$ (Why?). -->

Here, I define common parameters and the utility functions:

``` r
# y_low and y_high (I decide)
y_l <- 2 / 3
y_h <- 4 / 3
# p = 0.1 (given in the paper)
p <- 0.1
# utility functions are log form (and common for two households)
# Note that this utility functions cannot be changed in this document
# since the log form is used when I derive transfer amounts and
# the log form makes the hh and ll lambda-intervals identical.
u <- function(c){
  object <- log(c)
  return(object)
}
v <- function(c){
  object <- log(c)
  return(object)
}
```

I choose *y*<sub>*l*</sub> and *y*<sub>*h*</sub> so that the mean is 1 and *y*<sub>*l*</sub>/*y*<sub>*h*</sub> = 1/2 as in the original paper. The log utility functions and the equation (11) in the paper tell us the transfer amounts are given as *λ*: $\\tau\_s = \\frac{y\_1(s) - y\_2(s) \\lambda}{1 + \\lambda}$.

### Case 1: Autarky case

Under autarky, there is no transfer and thus the intervals become degenerate. Therefore, *λ*<sub>*h**h*, *l**o**w**e**r*</sub> = *λ*<sub>*h**h*, *u**p**p**e**r*</sub> = 1, *λ*<sub>*h**l*, *l**o**w**e**r*</sub> = *λ*<sub>*h**l*, *u**p**p**e**r*</sub> = 2, and *λ*<sub>*l**h*, *l**o**w**e**r*</sub> = *λ*<sub>*l**h*, *u**p**p**e**r*</sub> = 1/2.

### Case 2: Non-overlapping case

To calculate the end-points of the intervals, evaluate profits at the following end-points:

-   *λ*<sub>*l**h*, *u**p**p**e**r*</sub>, at which the discounted surplus of 2 is 0 and the next-period ratio of marginal utilities are *λ*<sub>*l**h*, *u**p**p**e**r*</sub>, *λ*<sub>*h**h*, *l**o**w**e**r*</sub>, *λ*<sub>*l**l*, *l**o**w**e**r*</sub>, *λ*<sub>*h**l*, *l**o**w**e**r*</sub>, with probabilities *p*(1 − *p*), (1 − *p*)<sup>2</sup>, *p*<sup>2</sup>, and *p*(1 − *p*), respectively.
-   *λ*<sub>*h**h*, *l**o**w**e**r*</sub>, at which the discounted surplus of 1 is 0 and the next-period ratio of marginal utilities are *λ*<sub>*l**h*, *u**p**p**e**r*</sub>, *λ*<sub>*h**h*, *l**o**w**e**r*</sub>, *λ*<sub>*l**l*, *l**o**w**e**r*</sub>, *λ*<sub>*h**l*, *l**o**w**e**r*</sub>, with probabilities *p*(1 − *p*), (1 − *p*)<sup>2</sup>, *p*<sup>2</sup>, and *p*(1 − *p*), respectively.
-   *λ*<sub>*h**h*, *u**p**p**e**r*</sub>, at which the discounted surplus of 2 is 0 and the next-period ratio of marginal utilities are *λ*<sub>*l**h*, *u**p**p**e**r*</sub>, *λ*<sub>*h**h*, *u**p**p**e**r*</sub>, *λ*<sub>*l**l*, *u**p**p**e**r*</sub>, *λ*<sub>*h**l*, *l**o**w**e**r*</sub>, with probabilities *p*(1 − *p*), (1 − *p*)<sup>2</sup>, *p*<sup>2</sup>, and *p*(1 − *p*), respectively.
-   *λ*<sub>*h**l*, *l**o**w**e**r*</sub>, at which the discounted surplus of 1 is 0 and the next-period ratio of marginal utilities are *λ*<sub>*l**h*, *u**p**p**e**r*</sub>, *λ*<sub>*h**h*, *u**p**p**e**r*</sub>, *λ*<sub>*l**l*, *u**p**p**e**r*</sub>, *λ*<sub>*h**l*, *l**o**w**e**r*</sub>, with probabilities *p*(1 − *p*), (1 − *p*)<sup>2</sup>, *p*<sup>2</sup>, and *p*(1 − *p*), respectively.

$\\begin{aligned}  0 & = v(y\_2(h) + (y\_1(l) - y\_2(h) \\lambda\_{lh, upper}) / (1 + \\lambda\_{lh, upper})) - v(y\_2(h)) + \\\\  & \\delta ((1 - 2p + 2p^2) V\_{hh, upper} + p (1 - p) V\_{hl, upper}) \\\\  U\_{lh, upper} & = u(y\_1(l) - (y\_1(l) - y\_2(h) \\lambda\_{lh, upper}) / (1 + \\lambda\_{lh, upper})) - u(y\_1(l)) + \\\\  & \\delta p (1 - p) U\_{lh, upper} \\\\  V\_{hh, upper} & = v(y\_2(h) + (y\_1(h) - y\_2(h) \\lambda\_{hh, lower}) / (1 + \\lambda\_{hh, lower})) - v(y\_2(h)) + \\\\  & \\delta ((1 - 2p + 2p^2) V\_{hh, upper} + p (1 - p) V\_{hl, upper}) \\\\  0 & = u(y\_1(h) - (y\_1(h) - y\_2(h) \\lambda\_{hh, lower}) / (1 + \\lambda\_{hh, lower})) - u(y\_1(h)) + \\\\  & \\delta p (1 - p) U\_{lh, upper} \\\\  0 & = v(y\_2(h) + (y\_1(h) - y\_2(h) \\lambda\_{hh, upper}) / (1 + \\lambda\_{hh, upper})) - v(y\_2(h)) + \\\\  & \\delta p (1 - p) V\_{hl, upper} \\\\  U\_{hh, upper} & = u(y\_1(h) - (y\_1(h) - y\_2(h) \\lambda\_{hh, upper}) / (1 + \\lambda\_{hh, upper})) - u(y\_1(h)) + \\\\  & \\delta (p (1 - p) U\_{lh, upper} + (1 - 2p + 2p^2) U\_{hh, upper}) \\\\  V\_{hl, upper} & = v(y\_2(l) + (y\_1(h) - y\_2(l) \\lambda\_{hl, lower}) / (1 + \\lambda\_{hl, lower})) - v(y\_2(l)) + \\\\  & \\delta p (1 - p) V\_{hl, upper} \\\\  0 & = u(y\_1(h) - (y\_1(h) - y\_2(l) \\lambda\_{hl, lower}) / (1 + \\lambda\_{hl, lower})) - u(y\_1(h)) + \\\\  & \\delta (p (1 - p) U\_{lh, upper} + (1 - 2p + 2p^2) U\_{hh, upper}) \\end{aligned}$

The following function solves this system of equations:

``` r
# For the case where there is no overlap, solve end-points and surpluses at these points

no_overlap <- function(x, delta, y_l, y_h, p){
  # x <- c(lambda_{lh, upper}, lambda_{hh, lower}, lambda_{hh, upper}, lambda_{hl, lower}, 
  #       U_{lh, upper}, V_{hh, upper}, U_{hh, upper}, V_{hl, upper})
  y <- numeric(8)
  y[1] <- v(y_h + (y_l - y_h * x[1]) / (1 + x[1])) - v(y_h) + 
            delta * ((1 - 2 * p + 2 * p^2) * x[6] + p * (1 - p) * x[8])
  y[2] <- u(y_l - (y_l - y_h * x[1]) / (1 + x[1])) - u(y_l) + 
            delta * (p * (1 - p) * x[5]) - x[5]
  y[3] <- v(y_h + (y_h - y_h * x[2]) / (1 + x[2])) - v(y_h) + 
            delta * ((1 - 2 * p + 2 * p^2) * x[6] + p * (1 - p) * x[8]) - x[6]
  y[4] <- u(y_h - (y_h - y_h * x[2]) / (1 + x[2])) - u(y_h) + 
            delta * (p * (1 - p) * x[5])
  y[5] <- v(y_h + (y_h - y_h * x[3]) / (1 + x[3])) - v(y_h) + 
            delta * (p * (1 - p) * x[8])
  y[6] <- u(y_h - (y_h - y_h * x[3]) / (1 + x[3])) - u(y_h) + 
            delta * (p * (1 - p) * x[5] + (1 - 2 * p + 2 * p^2) * x[7]) - x[7]
  y[7] <- v(y_l + (y_h - y_l * x[4]) / (1 + x[4])) - v(y_l) + 
            delta * (p * (1 - p) * x[8]) - x[8]
  y[8] <- u(y_h - (y_h - y_l * x[4]) / (1 + x[4])) - u(y_h) + 
            delta * (p * (1 - p) * x[5] + (1 - 2 * p + 2 * p^2) * x[7])
  return(y)
}
```

Given this function to solve each end-point, I obtain these end-points in this region, and the threshold with other regions (*δ*<sub>1</sub> (cutoff between autarky region and no-overlapping region) and *δ*<sub>2</sub> (cutoff between non-overlapping region and overlapping region)).

``` r
delta_vals <- seq(0.80, 0.999, by = 0.0001)
lambda_lh_upper_vals_no_overlap <- numeric(length(delta_vals))
lambda_hh_lower_vals_no_overlap <- numeric(length(delta_vals))
lambda_hh_upper_vals_no_overlap <- numeric(length(delta_vals))
lambda_hl_lower_vals_no_overlap <- numeric(length(delta_vals))
same_as_autarky <- numeric(length(delta_vals))
x_start <- rep((y_l + y_h) / 2, 8)
for (delta in 1:length(delta_vals)){
  no_overlap_solutions <- nleqslv(x_start, 
                                  no_overlap, 
                                  delta = delta_vals[delta],
                                  y_l = y_l, 
                                  y_h = y_h, 
                                  p = p)$x
  lambda_lh_upper_vals_no_overlap[delta] <- no_overlap_solutions[1]
  lambda_hh_lower_vals_no_overlap[delta] <- no_overlap_solutions[2]
  lambda_hh_upper_vals_no_overlap[delta] <- no_overlap_solutions[3]
  lambda_hl_lower_vals_no_overlap[delta] <- no_overlap_solutions[4]
  same_as_autarky[delta] <- is.logical(all.equal(no_overlap_solutions[1], y_l / y_h, tolerance = 1e-3))
}

delta_1 <- delta_vals[max(which(same_as_autarky == 1))]
delta_2 <- delta_vals[min(which(lambda_lh_upper_vals_no_overlap
                                >= lambda_hh_lower_vals_no_overlap))]
```

### Case 3: Overlapping case

To calculate the end-points of the intervals, evaluate profits at the following end-points:

-   *λ*<sub>*l**h*, *u**p**p**e**r*</sub>, at which the discounted surplus of 2 is 0 and the next-period ratio of marginal utilities are *λ*<sub>*l**h*, *u**p**p**e**r*</sub> and *λ*<sub>*h**l*, *l**o**w**e**r*</sub>, with probabilities 1 − *p*(1 − *p*) and *p*(1 − *p*), respectively.
-   *λ*<sub>*h**h*, *l**o**w**e**r*</sub>, at which the discounted surplus of 1 is 0 and the next-period ratio of marginal utilities are *λ*<sub>*h**h*, *l**o**w**e**r*</sub>, *λ*<sub>*l**l*, *l**o**w**e**r*</sub>, and *λ*<sub>*h**l*, *l**o**w**e**r*</sub> with probabilities (1 − *p*), *p*<sup>2</sup>, and *p*(1 − *p*), respectively.
-   *λ*<sub>*h**h*, *u**p**p**e**r*</sub>, at which the discounted surplus of 2 is 0 and the next-period ratio of marginal utilities are *λ*<sub>*l**h*, *u**p**p**e**r*</sub>, *λ*<sub>*h**h*, *u**p**p**e**r*</sub>, and *λ*<sub>*l**l*, *u**p**p**e**r*</sub> with probabilities *p*(1 − *p*), (1 − *p*), and *p*<sup>2</sup>, respectively.
-   *λ*<sub>*h**l*, *l**o**w**e**r*</sub>, at which the discounted surplus of 1 is 0 and the next-period ratio of marginal utilities are *λ*<sub>*l**h*, *u**p**p**e**r*</sub> and *λ*<sub>*h**l*, *l**o**w**e**r*</sub>, with probabilities *p*(1 − *p*) and 1 − *p*(1 − *p*), respectively.

This is very similar to the previous case, but this case is indeed more complicated. Why?

This is because, in the previous case I only had to consider discounted surpluses at the end-points, in this overlapping case that does not work. For instance, from *λ*<sub>*l**h*, *u**p**p**e**r*</sub>, if {High, High} state realized, then the ratio of marginal utilities remains *λ*<sub>*l**h*, *u**p**p**e**r*</sub> and it is interior of the intervals \[*λ*<sub>*h**h*, *l**o**w**e**r*</sub>, *λ*<sub>*h**h*, *u**p**p**e**r*</sub>\]. Therefore, we need to consider the discounted surpluses at that point too!

Here, I define the discounted surpluses at the interior of intervals as the followings:

-   *U*<sub>*h**h*</sub><sup>\*</sup>, *V*<sub>*h**h*</sub><sup>\*</sup>: discounted surpluses when income realization is {High, High} and the ratio of marginal utilities is *λ*<sub>*l**h*, *u**p**p**e**r*</sub>,
-   *U*<sub>*l**h*</sub><sup>\*</sup>, *V*<sub>*l**h*</sub><sup>\*</sup>: discounted surpluses when income realization is {Low, High} and the ratio of marginal utilities is *λ*<sub>*h**h*, *l**o**w**e**r*</sub>,
-   *U*<sub>*h**l*</sub><sup>\*</sup>, *V*<sub>*h**l*</sub><sup>\*</sup>: discounted surpluses when income realization is {High, Low} and the ratio of marginal utilities is *λ*<sub>*h**h*, *u**p**p**e**r*</sub>, and
-   *U*<sub>*h**h*</sub><sup>\*\*</sup>, *V*<sub>*h**h*</sub><sup>\*\*</sup>: discounted surpluses when income realization is {High, High} and the ratio of marginal utilities is *λ*<sub>*h**l*, *l**o**w**e**r*</sub>.

Then, I can write down the conditions:

$\\begin{aligned} 0 & = v(y\_2(h) + (y\_1(l) - y\_2(h) \\lambda\_{lh, upper}) / (1 + \\lambda\_{lh, upper})) - v(y\_2(h)) + \\\\  & \\delta ((1 - 2p + 2p^2) V\_{hh}^\* + p (1 - p) V\_{hl, upper}) \\\\ U\_{lh, upper} & = u(y\_1(l) - (y\_1(l) - y\_2(h) \\lambda\_{lh, upper}) / (1 + \\lambda\_{lh, upper})) - u(y\_1(l)) + \\\\  & \\delta ((1 - 2p + 2p^2) U\_{hh}^\* + p (1 - p) U\_{lh, upper}) \\\\ V\_{hh}^\* & = v(y\_2(h) + (y\_1(h) - y\_2(h) \\lambda\_{lh, upper}) / (1 + \\lambda\_{lh, upper})) - v(y\_2(h)) + \\\\  & \\delta (p (1 - p) V\_{hl, upper} + (1 - 2p + 2p^2) V\_{hh}^\*) \\\\ U\_{hh}^\* & = u(y\_1(h) - (y\_1(h) - y\_2(h) \\lambda\_{lh, upper}) / (1 + \\lambda\_{lh, upper})) - u(y\_1(h)) + \\\\  & \\delta (p (1 - p) U\_{lh, upper} + (1 - 2p + 2p^2) U\_{hh}^\*) \\\\ \\\\ V\_{hh, upper} & = v(y\_2(h) + (y\_1(h) - y\_2(h) \\lambda\_{hh, lower}) / (1 + \\lambda\_{hh, lower})) - v(y\_2(h)) + \\\\  & \\delta ((1 - 2p + 2p^2) V\_{hh, upper} + p (1 - p) V\_{lh}^\* + p (1 - p) V\_{hl, upper}) \\\\ 0 & = u(y\_1(h) - (y\_1(h) - y\_2(h) \\lambda\_{hh, lower}) / (1 + \\lambda\_{hh, lower})) - u(y\_1(h)) + \\\\  & \\delta p (1 - p) U\_{lh}^\* \\\\ V\_{lh}^\* & = v(y\_2(h) + (y\_1(l) - y\_2(h) \\lambda\_{hh, lower}) / (1 + \\lambda\_{hh, lower})) - v(y\_2(h)) + \\\\  & \\delta ((1 - 2p + 2p^2) V\_{hh, upper} + p (1 - p) V\_{lh}^\* + p (1 - p) V\_{hl, upper}) \\\\ U\_{lh}^\* & = u(y\_1(l) - (y\_1(l) - y\_2(h) \\lambda\_{hh, lower}) / (1 + \\lambda\_{hh, lower})) - u(y\_1(l)) + \\\\  & \\delta p (1 - p) U\_{lh}^\* \\\\ \\\\ 0 & = v(y\_2(h) + (y\_1(h) - y\_2(h) \\lambda\_{hh, upper}) / (1 + \\lambda\_{hh, upper})) - v(y\_2(h)) + \\\\  & \\delta p (1 - p) V\_{hl}^\* \\\\ U\_{hh, upper} & = u(y\_1(h) - (y\_1(h) - y\_2(h) \\lambda\_{hh, upper}) / (1 + \\lambda\_{hh, upper})) - u(y\_1(h)) + \\\\  & \\delta ((1 - 2p + 2p^2) U\_{hh, upper} + p (1 - p) U\_{hl}^\* + p (1 - p) U\_{lh, upper}) \\\\ V\_{hl}^\* & = v(y\_2(l) + (y\_1(h) - y\_2(l) \\lambda\_{hh, upper}) / (1 + \\lambda\_{hh, upper})) - v(y\_2(l)) + \\\\  & \\delta p (1 - p) V\_{hl}^\* \\\\ U\_{hl}^\* & = u(y\_1(h) - (y\_1(h) - y\_2(l) \\lambda\_{hh, upper}) / (1 + \\lambda\_{hh, upper})) - u(y\_1(h)) + \\\\  & \\delta ((1 - 2p + 2p^2) U\_{hh, upper} + p (1 - p) U\_{hl}^\* + p (1 - p) U\_{lh, upper}) \\\\ \\\\ V\_{hl, upper} & = v(y\_2(l) + (y\_1(h) - y\_2(l) \\lambda\_{hl, lower}) / (1 + \\lambda\_{hl, lower})) - v(y\_2(l)) + \\\\  & \\delta ((1 - 2p + 2p^2) V\_{hh}^{\*\*} + p (1 - p) V\_{hl, upper}) \\\\ 0 & = u(y\_1(h) - (y\_1(h) - y\_2(l) \\lambda\_{hl, lower}) / (1 + \\lambda\_{hl, lower})) - u(y\_1(h)) + \\\\  & \\delta ((1 - 2p + 2p^2) U\_{hh}^{\*\*} + p (1 - p) U\_{lh, upper}) \\\\ V\_{hh}^{\*\*} & = v(y\_2(h) + (y\_1(h) - y\_2(h) \\lambda\_{hl, lower}) / (1 + \\lambda\_{hl, lower})) - v(y\_2(h)) + \\\\  & \\delta (p (1 - p) V\_{hl, upper} + (1 - 2p + 2p^2) V\_{hh}^{\*\*}) \\\\ U\_{hh}^{\*\*} & = u(y\_1(h) - (y\_1(h) - y\_2(h) \\lambda\_{hl, lower}) / (1 + \\lambda\_{hl, lower})) - u(y\_1(h)) + \\\\  & \\delta (p (1 - p) U\_{lh, upper} + (1 - 2p + 2p^2) U\_{hh}^{\*\*}) \\\\ \\end{aligned}$

For the discounted surpluses at the interior of intervals, these relationships can be obtained:

$\\begin{aligned} V\_{hh}^\* & = \\frac{v(y\_2(h) + (y\_1(h) - y\_2(h) \\lambda\_{lh, upper}) / (1 + \\lambda\_{lh, upper})) - v(y\_2(h)) + \\delta p (1 - p) V\_{hl, upper}}{1 - \\delta (1 - 2p + 2p^2)} \\\\ U\_{hh}^\* & = \\frac{u(y\_1(h) - (y\_1(h) - y\_2(h) \\lambda\_{lh, upper}) / (1 + \\lambda\_{lh, upper})) - u(y\_1(h)) + \\delta p (1 - p) U\_{lh, upper}}{1 - \\delta (1 - 2p + 2p^2)} \\\\ V\_{lh}^\* & = \\frac{v(y\_2(h) + (y\_1(l) - y\_2(h) \\lambda\_{hh, lower}) / (1 + \\lambda\_{hh, lower})) - v(y\_2(h)) + \\delta ((1 - 2p + 2p^2) V\_{hh, upper} + p (1 - p) V\_{hl, upper})}{1 - \\delta p (1 - p)} \\\\ U\_{lh}^\* & = \\frac{u(y\_1(l) - (y\_1(l) - y\_2(h) \\lambda\_{hh, lower}) / (1 + \\lambda\_{hh, lower})) - u(y\_1(l))}{1 - \\delta p (1 - p)} \\\\ V\_{hl}^\* & = \\frac{v(y\_2(l) + (y\_1(h) - y\_2(l) \\lambda\_{hh, upper}) / (1 + \\lambda\_{hh, upper})) - v(y\_2(l))}{1 - \\delta p ( 1 - p)} \\\\ U\_{hl}^\* & = \\frac{u(y\_1(h) - (y\_1(h) - y\_2(l) \\lambda\_{hh, upper}) / (1 + \\lambda\_{hh, upper})) - u(y\_1(h)) + \\delta ((1 - 2p + 2p^2) U\_{hh, upper} + p (1 - p) U\_{lh, upper})}{1 - \\delta p (1 - p)} \\\\ V\_{hh}^{\*\*} & = \\frac{v(y\_2(h) + (y\_1(h) - y\_2(h) \\lambda\_{hl, lower}) / (1 + \\lambda\_{hl, lower})) - v(y\_2(h)) + \\delta p (1 - p) V\_{hl, upper}}{1 - \\delta (1 - 2p + 2p^2)} \\\\ U\_{hh}^{\*\*} & = \\frac{u(y\_1(h) - (y\_1(h) - y\_2(h) \\lambda\_{hl, lower}) / (1 + \\lambda\_{hl, lower})) - u(y\_1(h)) + \\delta p (1 - p) U\_{lh, upper}}{1 - \\delta (1 - 2p + 2p^2)} \\\\ \\end{aligned}$

Oh man this is A LOT!!! Finally I can define the function to solve this system of equations:

``` r
# For the case where there are overlaps, solve end-points and surpluses at these points

overlap <- function(x, delta, y_l, y_h, p){
  # x <- c(lambda_{lh, upper}, lambda_{hh, lower}, lambda_{hh, upper}, lambda_{hl, lower}, 
  #         U_{lh, upper}, V_{hh, upper}, U_{hh, upper}, V_{hl, upper})
  # First I derive the discounted surpluses at the interior of intervals
  V_hh_star  <- (v(y_h + (y_h - y_h * x[1]) / (1 + x[1])) - v(y_h) + 
                 delta * (p - p^2) * x[8]) / (1 - delta * (1 - 2 * p + 2 * p^2))
  U_hh_star  <- (u(y_h - (y_h - y_h * x[1]) / (1 + x[1])) - u(y_h) + 
                 delta * (p - p^2) * x[5]) / (1 - delta * (1 - 2 * p + 2 * p^2))

  V_lh_star  <- (v(y_h + (y_l - y_h * x[2]) / (1 + x[2])) - v(y_h) + 
                 delta * ((1 - 2 * p + 2 * p^2) * x[6] + (p - p^2) * x[8])) / (1 - delta * (p - p^2))
  U_lh_star  <- (u(y_l - (y_l - y_h * x[2]) / (1 + x[2])) - u(y_l)) / (1 - delta * (p - p^2))

  V_hl_star  <- (v(y_l + (y_h - y_l * x[3]) / (1 + x[3])) - v(y_l)) / (1 - delta * (p - p^2))
  U_hl_star  <- (u(y_h - (y_h - y_l * x[3]) / (1 + x[3])) - u(y_h) + 
                 delta * ((1 - 2 * p + 2 * p^2) * x[7] + (p - p^2) * x[5])) / (1 - delta * (p - p^2))
  
  V_hh_2star <- (v(y_h + (y_h - y_h * x[4]) / (1 + x[4])) - v(y_h) + 
                 delta * (p - p^2) * x[8]) / (1 - delta * (1 - 2 * p + 2 * p^2))
  U_hh_2star <- (u(y_h - (y_h - y_h * x[4]) / (1 + x[4])) - u(y_h) + 
                 delta * (p - p^2) * x[5]) / (1 - delta * (1 - 2 * p + 2 * p^2))


  y <- numeric(8)
  y[1] <- v(y_h + (y_l - y_h * x[1]) / (1 + x[1])) - v(y_h) + 
            delta * ((1 - 2 * p + 2 * p^2) * V_hh_star + (p - p^2) * x[8])
  y[2] <- u(y_l - (y_l - y_h * x[1]) / (1 + x[1])) - u(y_l) + 
            delta * ((1 - 2 * p + 2 * p^2) * U_hh_star + (p - p^2) * x[5]) - x[5]

  y[3] <- v(y_h + (y_h - y_h * x[2]) / (1 + x[2])) - v(y_h) + 
            delta * ((1 - p + p^2) * x[6] + (p - p^2) * V_lh_star + (p - p^2) * x[8])     - x[6]
  y[4] <- u(y_h - (y_h - y_h * x[2]) / (1 + x[2])) - u(y_h) + 
            delta * (p - p^2) * U_lh_star

  y[5] <- v(y_h + (y_h - y_h * x[3]) / (1 + x[3])) - v(y_h) + 
            delta * (p - p^2) * V_hl_star
  y[6] <- u(y_h - (y_h - y_h * x[3]) / (1 + x[3])) - u(y_h) + 
            delta * ((1 - p + p^2) * x[7] + (p - p^2) * U_hl_star + (p - p^2) * x[5])     - x[7]

  y[7] <- v(y_l + (y_h - y_l * x[4]) / (1 + x[4])) - v(y_l) + 
            delta * ((1 - 2 * p + 2 * p^2) * V_hh_2star + (p - p^2) * x[8]) - x[8]
  y[8] <- u(y_h - (y_h - y_l * x[4]) / (1 + x[4])) - u(y_h) + 
            delta * ((1 - 2 * p + 2 * p^2) * U_hh_2star + (p - p^2) * x[5])
  return(y)
}
```

Given this function to solve each end-point, I obtain these end-points in this region, and the cutoff between overlapping region and first-best region (*δ*<sub>3</sub>).

``` r
delta_vals <- seq(0.80, 0.999, by = 0.0001)
lambda_lh_upper_vals_overlap <- numeric(length(delta_vals))
lambda_hh_lower_vals_overlap <- numeric(length(delta_vals))
lambda_hh_upper_vals_overlap <- numeric(length(delta_vals))
lambda_hl_lower_vals_overlap <- numeric(length(delta_vals))
same_as_autarky <- numeric(length(delta_vals))
x_start <- rep((y_l + y_h) / 2, 8)
for (delta in 1:length(delta_vals)){
  overlap_solutions <- nleqslv(x_start, 
                               overlap, 
                               delta = delta_vals[delta],
                               y_l = y_l, 
                               y_h = y_h, 
                               p = p)$x
  lambda_lh_upper_vals_overlap[delta] <- overlap_solutions[1]
  lambda_hh_lower_vals_overlap[delta] <- overlap_solutions[2]
  lambda_hh_upper_vals_overlap[delta] <- overlap_solutions[3]
  lambda_hl_lower_vals_overlap[delta] <- overlap_solutions[4]
}

delta_3 <- delta_vals[min(which(lambda_lh_upper_vals_overlap
                        >= lambda_hl_lower_vals_overlap))]
```

### Case 4: First-best case

To calculate the end-points of the intervals, evaluate profits at the following end-points:

-   *λ*<sub>*l**h*, *u**p**p**e**r*</sub>, at which the discounted surplus of 2 is 0 and the next-period ratio of marginal utilities are *λ*<sub>*l**h*, *u**p**p**e**r*</sub> with probabilities 1.
-   *λ*<sub>*h**h*, *l**o**w**e**r*</sub>, at which the discounted surplus of 1 is 0 and the next-period ratio of marginal utilities are *λ*<sub>*h**h*, *l**o**w**e**r*</sub>, *λ*<sub>*l**l*, *l**o**w**e**r*</sub>, and *λ*<sub>*h**l*, *l**o**w**e**r*</sub> with probabilities (1 − *p*), *p*<sup>2</sup>, and *p*(1 − *p*), respectively.
-   *λ*<sub>*h**h*, *u**p**p**e**r*</sub>, at which the discounted surplus of 2 is 0 and the next-period ratio of marginal utilities are *λ*<sub>*l**h*, *u**p**p**e**r*</sub>, *λ*<sub>*h**h*, *u**p**p**e**r*</sub>, and *λ*<sub>*l**l*, *u**p**p**e**r*</sub> with probabilities *p*(1 − *p*), (1 − *p*), and *p*<sup>2</sup>, respectively.
-   *λ*<sub>*h**l*, *l**o**w**e**r*</sub>, at which the discounted surplus of 1 is 0 and the next-period ratio of marginal utilities are *λ*<sub>*h**l*, *l**o**w**e**r*</sub> with probabilities 1.

Here, I define the discounted surpluses at the interior of intervals as the followings:

-   *U*<sub>*h**h*</sub><sup>\*</sup>, *V*<sub>*h**h*</sub><sup>\*</sup>: discounted surpluses when income realization is {High, High} and the ratio of marginal utilities is *λ*<sub>*l**h*, *u**p**p**e**r*</sub>,
-   *U*<sub>*h**l*</sub><sup>\*</sup>, *V*<sub>*h**l*</sub><sup>\*</sup>: discounted surpluses when income realization is {High, Low} and the ratio of marginal utilities is *λ*<sub>*l**h*, *u**p**p**e**r*</sub>,
-   *U*<sub>*l**h*</sub><sup>\*\*</sup>, *V*<sub>*l**h*</sub><sup>\*\*</sup>: discounted surpluses when income realization is {Low, High} and the ratio of marginal utilities is *λ*<sub>*h**h*, *l**o**w**e**r*</sub>,
-   *U*<sub>*h**l*</sub><sup>\* \* \*</sup>, *V*<sub>*h**l*</sub><sup>\* \* \*</sup>: discounted surpluses when income realization is {Low, High} and the ratio of marginal utilities is *λ*<sub>*h**h*, *u**p**p**e**r*</sub>, and
-   *U*<sub>*h**h*</sub><sup>\* \* \*\*</sup>, *V*<sub>*h**h*</sub><sup>\* \* \*\*</sup>: discounted surpluses when income realization is {High, High} and the ratio of marginal utilities is *λ*<sub>*h**l*, *l**o**w**e**r*</sub>.
-   *U*<sub>*l**h*</sub><sup>\* \* \*\*</sup>, *V*<sub>*l**h*</sub><sup>\* \* \*\*</sup>: discounted surpluses when income realization is {Low, High} and the ratio of marginal utilities is *λ*<sub>*h**l*, *l**o**w**e**r*</sub>.

Then, I can write down the conditions:

$\\begin{aligned} 0 & = v(y\_2(h) + (y\_1(l) - y\_2(h) \\lambda\_{lh, upper}) / (1 + \\lambda\_{lh, upper})) - v(y\_2(h)) + \\\\  & \\delta ((1 - 2p + 2p^2) V\_{hh}^{\*} + p (1 - p) V\_{hl}^{\*}) \\\\ U\_{lh, upper} & = u(y\_1(l) - (y\_1(l) - y\_2(h) \\lambda\_{lh, upper}) / (1 + \\lambda\_{lh, upper})) - u(y\_1(l)) + \\\\  & \\delta ((1 - 2p + 2p^2) U\_{hh}^{\*} + p (1 - p) U\_{hl}^{\*} + p (1 - p) U\_{lh, upper}) \\\\ V\_{hh}^\* & = v(y\_2(h) + (y\_1(h) - y\_2(h) \\lambda\_{lh, upper}) / (1 + \\lambda\_{lh, upper})) - v(y\_2(h)) + \\\\  & \\delta ((1 - 2p + 2p^2) V\_{hh}^{\*} + p (1 - p) V\_{hl}^{\*}) \\\\ U\_{hh}^\* & = u(y\_1(h) - (y\_1(h) - y\_2(h) \\lambda\_{lh, upper}) / (1 + \\lambda\_{lh, upper})) - u(y\_1(h)) + \\\\  & \\delta ((1 - 2p + 2p^2) U\_{hh}^{\*} + p (1 - p) U\_{hl}^{\*} + p (1 - p) U\_{lh, upper}) \\\\ V\_{hl}^\* & = v(y\_2(l) + (y\_1(h) - y\_2(l) \\lambda\_{lh, upper}) / (1 + \\lambda\_{lh, upper})) - v(y\_2(l)) + \\\\  & \\delta ((1 - 2p + 2p^2) V\_{hh}^{\*} + p (1 - p) V\_{hl}^{\*}) \\\\ U\_{hl}^\* & = u(y\_1(h) - (y\_1(h) - y\_2(l) \\lambda\_{lh, upper}) / (1 + \\lambda\_{lh, upper})) - u(y\_1(h)) + \\\\  & \\delta ((1 - 2p + 2p^2) U\_{hh}^{\*} + p (1 - p) U\_{hl}^{\*} + p (1 - p) U\_{lh, upper}) \\\\ \\\\ V\_{hh, upper} & = v(y\_2(h) + (y\_1(h) - y\_2(h) \\lambda\_{hh, lower}) / (1 + \\lambda\_{hh, lower})) - v(y\_2(h)) + \\\\  & \\delta ((1 - 2p + 2p^2) V\_{hh, upper} + p (1 - p) V\_{lh}^{\*\*} + p (1 - p) V\_{hl, upper}) \\\\ 0 & = u(y\_1(h) - (y\_1(h) - y\_2(h) \\lambda\_{hh, lower}) / (1 + \\lambda\_{hh, lower})) - u(y\_1(h)) + \\\\  & \\delta p (1 - p) U\_{lh}^{\*\*} \\\\ V\_{lh}^{\*\*} & = v(y\_2(h) + (y\_1(l) - y\_2(h) \\lambda\_{hh, lower}) / (1 + \\lambda\_{hh, lower})) - v(y\_2(h)) + \\\\  & \\delta ((1 - 2p + 2p^2) V\_{hh, upper} + p (1 - p) V\_{lh}^{\*\*} + p (1 - p) V\_{hl, upper}) \\\\ U\_{lh}^{\*\*} & = u(y\_1(l) - (y\_1(l) - y\_2(h) \\lambda\_{hh, lower}) / (1 + \\lambda\_{hh, lower})) - u(y\_1(l)) + \\\\  & \\delta p (1 - p) U\_{lh}^{\*\*} \\\\ \\\\ 0 & = v(y\_2(h) + (y\_1(h) - y\_2(h) \\lambda\_{hh, upper}) / (1 + \\lambda\_{hh, upper})) - v(y\_2(h)) + \\\\  & \\delta p (1 - p) V\_{hl}^{\*\*\*} \\\\ U\_{hh, upper} & = u(y\_1(h) - (y\_1(h) - y\_2(h) \\lambda\_{hh, upper}) / (1 + \\lambda\_{hh, upper})) - u(y\_1(h)) + \\\\  & \\delta ((1 - 2p + 2p^2) U\_{hh, upper} + p (1 - p) U\_{hl}^{\*\*\*} + p (1 - p) U\_{lh, upper}) \\\\ V\_{hl}^{\*\*\*} & = v(y\_2(l) + (y\_1(h) - y\_2(l) \\lambda\_{hh, upper}) / (1 + \\lambda\_{hh, upper})) - v(y\_2(l)) + \\\\  & \\delta p (1 - p) V\_{hl}^\* \\\\ U\_{hl}^{\*\*\*} & = u(y\_1(h) - (y\_1(h) - y\_2(l) \\lambda\_{hh, upper}) / (1 + \\lambda\_{hh, upper})) - u(y\_1(h)) + \\\\  & \\delta ((1 - 2p + 2p^2) U\_{hh, upper} + p (1 - p) U\_{hl}^\* + p (1 - p) U\_{lh, upper}) \\\\ \\\\ V\_{hl, upper} & = v(y\_2(l) + (y\_1(h) - y\_2(l) \\lambda\_{hl, lower}) / (1 + \\lambda\_{hl, lower})) - v(y\_2(l)) + \\\\  & \\delta ((1 - 2p + 2p^2) V\_{hh}^{\*\*\*\*} + p (1 - p) V\_{lh}^{\*\*\*\*} + p (1 - p) V\_{hl, upper}) \\\\ 0 & = u(y\_1(h) - (y\_1(h) - y\_2(l) \\lambda\_{hl, lower}) / (1 + \\lambda\_{hl, lower})) - u(y\_1(h)) + \\\\  & \\delta ((1 - 2p + 2p^2) U\_{hh}^{\*\*\*\*} + p (1 - p) U\_{lh}^{\*\*\*\*}) \\\\ V\_{hh}^{\*\*\*\*} & = v(y\_2(h) + (y\_1(h) - y\_2(h) \\lambda\_{hl, lower}) / (1 + \\lambda\_{hl, lower})) - v(y\_2(h)) + \\\\  & \\delta ((1 - 2p + 2p^2) V\_{hh}^{\*\*\*\*} + p (1 - p) V\_{lh}^{\*\*\*\*} + p (1 - p) V\_{hl, upper}) \\\\ U\_{hh}^{\*\*\*\*} & = u(y\_1(h) - (y\_1(h) - y\_2(h) \\lambda\_{hl, lower}) / (1 + \\lambda\_{hl, lower})) - u(y\_1(h)) + \\\\  & \\delta ((1 - 2p + 2p^2) U\_{hh}^{\*\*\*\*} + p (1 - p) U\_{lh}^{\*\*\*\*}) \\\\ V\_{lh}^{\*\*\*\*} & = v(y\_2(h) + (y\_1(l) - y\_2(h) \\lambda\_{hl, lower}) / (1 + \\lambda\_{hl, lower})) - v(y\_2(h)) + \\\\  & \\delta ((1 - 2p + 2p^2) V\_{hh}^{\*\*\*\*} + p (1 - p) V\_{lh}^{\*\*\*\*} + p (1 - p) V\_{hl, upper}) \\\\ U\_{lh}^{\*\*\*\*} & = u(y\_1(l) - (y\_1(l) - y\_2(h) \\lambda\_{hl, lower}) / (1 + \\lambda\_{hl, lower})) - u(y\_1(l)) + \\\\  & \\delta ((1 - 2p + 2p^2) U\_{hh}^{\*\*\*\*} + p (1 - p) U\_{lh}^{\*\*\*\*}) \\end{aligned}$

For the discounted surpluses at the interior of intervals, these relationships can be obtained (notice that some of the future expected surpluses take the same form due to unchanged *λ*'s):

$\\begin{aligned}  V\_{hh}^\* &= v(y\_2(h) + (y\_1(h) - y\_2(h) \\lambda\_{lh, upper}) / (1 + \\lambda\_{lh, upper})) - v(y\_2(h)) - \\\\  &(v(y\_2(h) + (y\_1(l) - y\_2(h) \\lambda\_{lh, upper}) / (1 + \\lambda\_{lh, upper})) - v(y\_2(h))) \\\\  U\_{hh}^\* &= u(y\_1(h) - (y\_1(h) - y\_2(h) \\lambda\_{lh, upper}) / (1 + \\lambda\_{lh, upper})) - u(y\_1(h)) + \\\\  &(U\_{lh, upper} - (u(y\_1(l) - (y\_1(l) - y\_2(h) \\lambda\_{lh, upper}) / (1 + \\lambda\_{lh, upper})) - u(y\_1(l)))) \\\\  V\_{hl}^\* &= v(y\_2(l) + (y\_1(h) - y\_2(l) \\lambda\_{lh, upper}) / (1 + \\lambda\_{lh, upper})) - v(y\_2(l)) - \\\\  &(v(y\_2(h) + (y\_1(l) - y\_2(h) \\lambda\_{lh, upper}) / (1 + \\lambda\_{lh, upper})) - v(y\_2(h))) \\\\  U\_{hl}^\* &= u(y\_1(h) - (y\_1(h) - y\_2(l) \\lambda\_{lh, upper}) / (1 + \\lambda\_{lh, upper})) - u(y\_1(h)) + \\\\  &(U\_{lh, upper} - (u(y\_1(l) - (y\_1(l) - y\_2(h) \\lambda\_{lh, upper}) / (1 + \\lambda\_{lh, upper})) - u(y\_1(l)))) \\\\  \\\\  V\_{lh}^{\*\*} &= v(y\_2(h) + (y\_1(l) - y\_2(h) \\lambda\_{hh, lower}) / (1 + \\lambda\_{hh, lower})) - v(y\_2(h)) + \\\\  &(V\_{hh, upper} - (v(y\_2(h) + (y\_1(h) - y\_2(h) \\lambda\_{hh, lower}) / (1 + \\lambda\_{hh, lower})) - v(y\_2(h)))) \\\\  U\_{lh}^{\*\*} &= u(y\_1(l) - (y\_1(l) - y\_2(h) \\lambda\_{hh, lower}) / (1 + \\lambda\_{hh, lower})) - u(y\_1(l)) - \\\\  &(u(y\_1(h) - (y\_1(h) - y\_2(h) \\lambda\_{hh, lower}) / (1 + \\lambda\_{hh, lower})) - u(y\_1(h))) \\\\  \\\\  V\_{hl}^{\*\*\*} &= v(y\_2(l) + (y\_1(h) - y\_2(l) \\lambda\_{hh, upper}) / (1 + \\lambda\_{hh, upper})) - v(y\_2(l)) - \\\\  &(v(y\_2(h) + (y\_1(h) - y\_2(h) \\lambda\_{hh, upper}) / (1 + \\lambda\_{hh, upper})) - v(y\_2(h))) \\\\  U\_{hl}^{\*\*\*} &= u(y\_1(h) - (y\_1(h) - y\_2(l) \\lambda\_{hh, upper}) / (1 + \\lambda\_{hh, upper})) - u(y\_1(h)) + \\\\  &(U\_{hh, upper} - (u(y\_1(h) - (y\_1(h) - y\_2(h) \\lambda\_{hh, upper}) / (1 + \\lambda\_{hh, upper})) - u(y\_1(h)))) \\\\  \\\\  V\_{hh}^{\*\*\*\*} &= v(y\_2(h) + (y\_1(h) - y\_2(h) \\lambda\_{hl, lower}) / (1 + \\lambda\_{hl, lower})) - v(y\_2(h)) + \\\\  &(V\_{hl, upper} - (v(y\_2(l) + (y\_1(h) - y\_2(l) \\lambda\_{hl, lower}) / (1 + \\lambda\_{hl, lower})) - v(y\_2(l)))) \\\\  U\_{hh}^{\*\*\*\*} &= u(y\_1(h) - (y\_1(h) - y\_2(h) \\lambda\_{hl, lower}) / (1 + \\lambda\_{hl, lower})) - u(y\_1(h)) - \\\\  &(u(y\_1(h) - (y\_1(h) - y\_2(l) \\lambda\_{hl, lower}) / (1 + \\lambda\_{hl, lower})) - u(y\_1(h))) \\\\  V\_{lh}^{\*\*\*\*} &= v(y\_2(h) + (y\_1(l) - y\_2(h) \\lambda\_{hl, lower}) / (1 + \\lambda\_{hl, lower})) - v(y\_2(h)) + \\\\  &(V\_{hl, upper} - (v(y\_2(l) + (y\_1(h) - y\_2(l) \\lambda\_{hl, lower}) / (1 + \\lambda\_{hl, lower})) - v(y\_2(l)))) \\\\  U\_{lh}^{\*\*\*\*} &= u(y\_1(l) - (y\_1(l) - y\_2(h) \\lambda\_{hl, lower}) / (1 + \\lambda\_{hl, lower})) - u(y\_1(l)) - \\\\  &(u(y\_1(h) - (y\_1(h) - y\_2(l) \\lambda\_{hl, lower}) / (1 + \\lambda\_{hl, lower})) - u(y\_1(h))) \\end{aligned}$

The following function solves this system of equations:

``` r
# For the case where first best can be achieved, solve end-points and surpluses at these points

first_best <- function(x, delta, y_l, y_h, p){
  # x <- c(lambda_{lh, upper}, lambda_{hh, lower}, lambda_{hh, upper}, lambda_{hl, lower}, 
  #       U_{lh, upper}, V_{hh, upper}, U_{hh, upper}, V_{hl, upper})
  # First I derive the discounted surpluses at the interior of intervals

  V_hh_star  = v(y_h + (y_h - y_h * x[1]) / (1 + x[1])) - v(y_h) -
                (v(y_h + (y_l - y_h * x[1]) / (1 + x[1])) - v(y_h))
  U_hh_star  = u(y_h - (y_h - y_h * x[1]) / (1 + x[1])) - u(y_h) + 
                (x[5] - (u(y_l - (y_l - y_h * x[1]) / (1 + x[1])) - u(y_l)))
  V_hl_star  = v(y_l + (y_h - y_l * x[1]) / (1 + x[1])) - v(y_l) -
                (v(y_h + (y_l - y_h * x[1]) / (1 + x[1])) - v(y_h))
  U_hl_star  = u(y_h - (y_h - y_l * x[1]) / (1 + x[1])) - u(y_h) +
                (x[5] - (u(y_l - (y_l - y_h * x[1]) / (1 + x[1])) - u(y_l)))
  V_lh_2star = v(y_h + (y_l - y_h * x[2]) / (1 + x[2])) - v(y_h) + 
                (x[6] - (v(y_h + (y_h - y_h * x[2]) / (1 + x[2])) - v(y_h))) 
  U_lh_2star = u(y_l - (y_l - y_h * x[2]) / (1 + x[2])) - u(y_l) -
                (u(y_h - (y_h - y_h * x[2]) / (1 + x[2])) - u(y_h))
  V_hl_3star = v(y_l + (y_h - y_l * x[3]) / (1 + x[3])) - v(y_l) -
                (v(y_h + (y_h - y_h * x[3]) / (1 + x[3])) - v(y_h))
  U_hl_3star = u(y_h - (y_h - y_l * x[3]) / (1 + x[3])) - u(y_h) +
                (x[7] - (u(y_h - (y_h - y_h * x[3]) / (1 + x[3])) - u(y_h)))
  V_hh_4star = v(y_h + (y_h - y_h * x[4]) / (1 + x[4])) - v(y_h) +
                (x[8] - (v(y_l + (y_h - y_l * x[4]) / (1 + x[4])) - v(y_l)))
  U_hh_4star = u(y_h - (y_h - y_h * x[4]) / (1 + x[4])) - u(y_h) -
                (u(y_h - (y_h - y_l * x[4]) / (1 + x[4])) - u(y_h))
  V_lh_4star = v(y_h + (y_l - y_h * x[4]) / (1 + x[4])) - v(y_h) +
                (x[8] - (v(y_l + (y_h - y_l * x[4]) / (1 + x[4])) - v(y_l)))
  U_lh_4star = u(y_l - (y_l - y_h * x[4]) / (1 + x[4])) - u(y_l) - 
                (u(y_h - (y_h - y_l * x[4]) / (1 + x[4])) - u(y_h))


  y <- numeric(8)
  y[1] <- v(y_h + (y_l - y_h * x[1]) / (1 + x[1])) - v(y_h) + 
            delta * ((1 - 2 * p + 2 * p^2) * V_hh_star + (p - p^2) * V_hl_star)
  y[2] <- u(y_l - (y_l - y_h * x[1]) / (1 + x[1])) - u(y_l) + 
            delta * ((1 - 2 * p + 2 * p^2) * U_hh_star + (p - p^2) * U_hl_star + 
                     (p - p^2) * x[5]) - x[5]
  y[3] <- v(y_h + (y_h - y_h * x[2]) / (1 + x[2])) - v(y_h) + 
            delta * ((1 - 2 * p + 2 * p^2) * x[6] + (p - p^2) * x[8] + 
                     (p - p^2) * V_lh_2star) - x[6]
  y[4] <- u(y_h - (y_h - y_h * x[2]) / (1 + x[2])) - u(y_h) + 
            delta * (p - p^2) * U_lh_2star
  y[5] <- v(y_h + (y_h - y_h * x[3]) / (1 + x[3])) - v(y_h) + 
            delta * (p - p^2) * V_hl_3star
  y[6] <- u(y_h - (y_h - y_h * x[3]) / (1 + x[3])) - u(y_h) + 
            delta * ((1 - 2 * p + 2 * p^2) * x[5] + (p - p^2) * x[7] + 
                     (p - p^2) * U_hl_3star) - x[7]
  y[7] <- v(y_l + (y_h - y_l * x[4]) / (1 + x[4])) - v(y_l) + 
            delta * ((1 - 2 * p + 2 * p^2) * V_hh_4star + (p - p^2) * V_lh_4star + 
                     (p - p^2) * x[8]) - x[8]
  y[8] <- u(y_h - (y_h - y_l * x[4]) / (1 + x[4])) - u(y_h) + 
            delta * ((1 - 2 * p + 2 * p^2) * U_hh_4star + (p - p^2) * U_lh_4star)
  return(y)
}
```

Given this function to solve each end-point, I obtain these end-points in this region:

``` r
delta_vals <- seq(0.80, 0.999, by = 0.0001)
lambda_lh_upper_vals_first_best <- numeric(length(delta_vals))
lambda_hh_lower_vals_first_best <- numeric(length(delta_vals))
lambda_hh_upper_vals_first_best <- numeric(length(delta_vals))
lambda_hl_lower_vals_first_best <- numeric(length(delta_vals))
same_as_autarky <- numeric(length(delta_vals))
x_start <- rep((y_l + y_h) / 2, 8)
for (delta in 1:length(delta_vals)){
  first_best_solutions <- nleqslv(x_start, 
                               first_best, 
                               delta = delta_vals[delta],
                               y_l = y_l, 
                               y_h = y_h, 
                               p = p)$x
  lambda_lh_upper_vals_first_best[delta] <- first_best_solutions[1]
  lambda_hh_lower_vals_first_best[delta] <- first_best_solutions[2]
  lambda_hh_upper_vals_first_best[delta] <- first_best_solutions[3]
  lambda_hl_lower_vals_first_best[delta] <- first_best_solutions[4]
}
```

### Replicate the figure

With the information at hand, I replicate Figure 1 in the original paper:

``` r
delta_vals <- seq(0.80, 0.999, by = 0.0001)

lambda_lh_upper_vals <- numeric(length(delta_vals))
lambda_hh_lower_vals <- numeric(length(delta_vals))
lambda_hh_upper_vals <- numeric(length(delta_vals))
lambda_hl_lower_vals <- numeric(length(delta_vals))
x_start <- rep((y_l + y_h) / 2, 8)
for (delta in 1:length(delta_vals)){
  if (delta_vals[delta] <= delta_1){
    lambda_lh_upper_vals[delta] <- y_l / y_h 
    lambda_hh_lower_vals[delta] <- 1
    lambda_hh_upper_vals[delta] <- 1
    lambda_hl_lower_vals[delta] <- y_h / y_l 
  } else if (delta_vals[delta] > delta_1 & delta_vals[delta] <= delta_2){
    lambda_lh_upper_vals[delta] <- lambda_lh_upper_vals_no_overlap[delta]
    lambda_hh_lower_vals[delta] <- lambda_hh_lower_vals_no_overlap[delta]
    lambda_hh_upper_vals[delta] <- lambda_hh_upper_vals_no_overlap[delta]
    lambda_hl_lower_vals[delta] <- lambda_hl_lower_vals_no_overlap[delta]
  } else if (delta_vals[delta] > delta_2 & delta_vals[delta] <= delta_3){
    lambda_lh_upper_vals[delta] <- lambda_lh_upper_vals_overlap[delta]
    lambda_hh_lower_vals[delta] <- lambda_hh_lower_vals_overlap[delta]
    lambda_hh_upper_vals[delta] <- lambda_hh_upper_vals_overlap[delta]
    lambda_hl_lower_vals[delta] <- lambda_hl_lower_vals_overlap[delta]
  } else {
    lambda_lh_upper_vals[delta] <- lambda_lh_upper_vals_first_best[delta]
    lambda_hh_lower_vals[delta] <- lambda_hh_lower_vals_first_best[delta]
    lambda_hh_upper_vals[delta] <- lambda_hh_upper_vals_first_best[delta]
    lambda_hl_lower_vals[delta] <- lambda_hl_lower_vals_first_best[delta]
  }
}

lambda_lh_lower_vals <- rep(y_l / y_h, length(delta_vals))
lambda_hl_upper_vals <- rep(y_h / y_l, length(delta_vals))

ggplot() +
  geom_line(aes(delta_vals, log(lambda_lh_lower_vals), color="a")) + 
  geom_line(aes(delta_vals, log(lambda_lh_upper_vals), color="b")) + 
  geom_line(aes(delta_vals, log(lambda_hh_lower_vals), color="c")) + 
  geom_line(aes(delta_vals, log(lambda_hh_upper_vals), color="d")) +
  geom_line(aes(delta_vals, log(lambda_hl_lower_vals), color="e")) +
  geom_line(aes(delta_vals, log(lambda_hl_upper_vals), color="f")) +
  coord_cartesian(xlim = c(0.8, 1.0)) +
  geom_ribbon(aes(x = delta_vals, 
                  ymin = log(lambda_lh_lower_vals),
                  ymax = log(lambda_lh_upper_vals)),
                  fill = "blue", alpha = "0.2") +
  geom_ribbon(aes(x = delta_vals,
                  ymin = log(lambda_hh_lower_vals),
                  ymax = log(lambda_hh_upper_vals)),
                  fill = "red", alpha = "0.2") +
  geom_ribbon(aes(x = delta_vals,
                  ymin = log(lambda_hl_lower_vals),
                  ymax = log(lambda_hl_upper_vals)),
                  fill = "green", alpha = "0.2") + 
  scale_color_manual(name = "End-points", 
                     values = c("blue", "purple", "brown", "red", "yellow", "green"), 
                     labels = c("lambda_lh_lower", 
                                "lambda_lh_upper",
                                "lambda_hh_lower",
                                "lambda_hh_upper",
                                "lambda_hl_lower",
                                "lambda_hl_upper")) +
  xlab("Discount factor (delta)") +
  ylab("log of the ratio of marginal utilities (lambda)")
```

![](LTW_code_files/figure-markdown_github/unnamed-chunk-8-1.png)

I could replicate the figure! Yay!!

The cutoffs of the regions are *δ*<sub>1</sub> = 0.8666, *δ*<sub>2</sub> = 0.9351, and *δ*<sub>3</sub> = 0.9645.

Fun Exercise: How does the ratio of marginal utilities change over time in response to income shocks?
-----------------------------------------------------------------------------------------------------

Now that we could replicate the figure in the original paper, I will next see how *λ* changes over time in response to income shocks. This is the function to create the figure:

``` r
lambda_change_plot <- function(lambda_init, delta, n){


  # Sequence of income shocks
  income_realization <- c(1, 2, 3, 4)
  income_realization_label <- c("High, High", "High, Low", "Low, High", "Low, Low")
  income_seq <- sample(income_realization, 
                       size = n,
                       replace = TRUE,
                       prob = c((1 - p)^2,
                                p * (1 - p),
                                p * (1 - p),
                                p^2))

  # Sequence of the ratios of marginal utilities under autarky 
  xi_seq <- numeric(n)
  income_realization_label_seq <- rep("", n)
  for (i in 1:n){
    xi_seq[i] <- 1 * (income_seq[i] == 1) + 
                  (y_h / y_l) * (income_seq[i] == 2) +
                  (y_l / y_h) * (income_seq[i] == 3) +
                  1 * (income_seq[i] == 4)
    income_realization_label_seq[i] <- income_realization_label[income_seq[i]]
  }
  
  # Given the initial lambda, derive the sequence of lambda's
  
  if (delta <= delta_1){
    lambda_lh_upper <- y_l / y_h
    lambda_hh_lower <- 1
    lambda_hh_upper <- 1
    lambda_hl_lower <- y_h / y_l
  } else if (delta > delta_1 & delta <= delta_2){
    no_overlap_solutions <- nleqslv(x_start, 
                                    no_overlap, 
                                    delta = delta,
                                    y_l = y_l, 
                                    y_h = y_h, 
                                    p = p)$x
    lambda_lh_upper <- no_overlap_solutions[1]
    lambda_hh_lower <- no_overlap_solutions[2]
    lambda_hh_upper <- no_overlap_solutions[3]
    lambda_hl_lower <- no_overlap_solutions[4]
  } else if (delta > delta_2 & delta <= delta_3) {
    overlap_solutions <- nleqslv(x_start, 
                                 overlap, 
                                 delta = delta,
                                 y_l = y_l, 
                                 y_h = y_h, 
                                 p = p)$x
    lambda_lh_upper <- overlap_solutions[1]
    lambda_hh_lower <- overlap_solutions[2]
    lambda_hh_upper <- overlap_solutions[3]
    lambda_hl_lower <- overlap_solutions[4]
  } else {
    first_best_solutions <- nleqslv(x_start, 
                                 first_best, 
                                 delta = delta,
                                 y_l = y_l, 
                                 y_h = y_h, 
                                 p = p)$x
    lambda_lh_upper <- first_best_solutions[1]
    lambda_hh_lower <- first_best_solutions[2]
    lambda_hh_upper <- first_best_solutions[3]
    lambda_hl_lower <- first_best_solutions[4]
  }
  
  lambda <- lambda_init
  lambda_seq <- numeric(n)

  if (delta <= delta_1){
    lambda_seq <- xi_seq
  } else if (delta > delta_1 & delta <= delta_2){
    for (i in 1:n){
      if (lambda >= lambda_hl_lower){
        lambda_seq[i] <- lambda * (income_seq[i] == 2) +
                        lambda_hh_upper * (income_seq[i] == 1 | income_seq[i] == 4) +
                        lambda_lh_upper * (income_seq[i] == 3)
        lambda <- lambda_seq[i]
      } else if (lambda < lambda_hl_lower & lambda >= lambda_hh_upper){
        lambda_seq[i] <- lambda_hl_lower * (income_seq[i] == 2) +
                        lambda_hh_upper * (income_seq[i] == 1 | income_seq[i] == 4) +
                        lambda_lh_upper * (income_seq[i] == 3)
        lambda <- lambda_seq[i]
      } else if (lambda < lambda_hh_upper & lambda >= lambda_hh_lower){
        lambda_seq[i] <- lambda_hl_lower * (income_seq[i] == 2) +
                        lambda * (income_seq[i] == 1 | income_seq[i] == 4) +
                        lambda_lh_upper * (income_seq[i] == 3)
        lambda <- lambda_seq[i]
      } else if (lambda < lambda_hh_lower & lambda >= lambda_lh_upper){
        lambda_seq[i] <- lambda_hl_lower * (income_seq[i] == 2) +
                        lambda_hh_lower * (income_seq[i] == 1 | income_seq[i] == 4) +
                        lambda_lh_upper * (income_seq[i] == 3)
        lambda <- lambda_seq[i]
      } else {
        lambda_seq[i] <- lambda_hl_lower * (income_seq[i] == 2) +
                        lambda_hh_lower * (income_seq[i] == 1 | income_seq[i] == 4) +
                        lambda * (income_seq[i] == 3)
        lambda <- lambda_seq[i]
      }
    }
  } else if (delta > delta_2 & delta <= delta_3) {
    for (i in 1:n){
      if (lambda >= lambda_hh_upper){
        lambda_seq[i] <- lambda * (income_seq[i] == 2) +
                        lambda_hh_upper * (income_seq[i] == 1 | income_seq[i] == 4) +
                        lambda_lh_upper * (income_seq[i] == 3)
        lambda <- lambda_seq[i]
      } else if (lambda < lambda_hh_upper & lambda >= lambda_hl_lower){
        lambda_seq[i] <- lambda * (income_seq[i] == 2) +
                        lambda * (income_seq[i] == 1 | income_seq[i] == 4) +
                        lambda_lh_upper * (income_seq[i] == 3)
        lambda <- lambda_seq[i]
      } else if (lambda < lambda_hl_lower & lambda >= lambda_lh_upper){
        lambda_seq[i] <- lambda_hl_lower * (income_seq[i] == 2) +
                        lambda * (income_seq[i] == 1 | income_seq[i] == 4) +
                        lambda_lh_upper * (income_seq[i] == 3)
        lambda <- lambda_seq[i]
      } else if (lambda < lambda_lh_upper & lambda >= lambda_hh_lower){
        lambda_seq[i] <- lambda_hl_lower * (income_seq[i] == 2) +
                        lambda * (income_seq[i] == 1 | income_seq[i] == 4) +
                        lambda * (income_seq[i] == 3)
        lambda <- lambda_seq[i]
      } else {
        lambda_seq[i] <- lambda_hl_lower * (income_seq[i] == 2) +
                        lambda_hh_lower * (income_seq[i] == 1 | income_seq[i] == 4) +
                        lambda * (income_seq[i] == 3)
        lambda <- lambda_seq[i]
      }
    }
  } else {
    for (i in 1:n){
      if (lambda >= lambda_hh_upper){
        lambda_seq[i] <- lambda * (income_seq[i] == 2) +
                        lambda_hh_upper * (income_seq[i] == 1 | income_seq[i] == 4) +
                        lambda_lh_upper * (income_seq[i] == 3)
        lambda <- lambda_seq[i]
      } else if (lambda < lambda_hh_upper & lambda >= lambda_lh_upper){
        lambda_seq[i] <- lambda * (income_seq[i] == 2) +
                        lambda * (income_seq[i] == 1 | income_seq[i] == 4) +
                        lambda_lh_upper * (income_seq[i] == 3)
        lambda <- lambda_seq[i]
      } else if (lambda < lambda_lh_upper & lambda >= lambda_hl_lower){
        lambda_seq[i] <- lambda * (income_seq[i] == 2) +
                        lambda * (income_seq[i] == 1 | income_seq[i] == 4) +
                        lambda * (income_seq[i] == 3)
        lambda <- lambda_seq[i]
      } else if (lambda < lambda_hl_lower & lambda >= lambda_hh_lower){
        lambda_seq[i] <- lambda_hl_lower * (income_seq[i] == 2) +
                        lambda * (income_seq[i] == 1 | income_seq[i] == 4) +
                        lambda * (income_seq[i] == 3)
        lambda <- lambda_seq[i]
      } else {
        lambda_seq[i] <- lambda_hl_lower * (income_seq[i] == 2) +
                        lambda_hh_lower * (income_seq[i] == 1 | income_seq[i] == 4) +
                        lambda * (income_seq[i] == 3)
        lambda <- lambda_seq[i]
      }
    }
  }
  # Sequence of net transfers (1 -> 2), income, and consumptions
  transfer_seq <- numeric(n)
  inc_seq_1 <- numeric(n)
  inc_seq_2 <- numeric(n)
  cons_seq_1 <- numeric(n)
  cons_seq_2 <- numeric(n)
  for (i in 1:n){
    if (income_seq[i] == 1){
      transfer_seq[i] <- (y_h - y_h * lambda_seq[i]) / (1 + lambda_seq[i])
      inc_seq_1[i] <- y_h
      inc_seq_2[i] <- y_h
      cons_seq_1[i] <- y_h - transfer_seq[i]
      cons_seq_2[i] <- y_h + transfer_seq[i]
    } else if (income_seq[i] == 2) {
      transfer_seq[i] <- (y_h - y_l * lambda_seq[i]) / (1 + lambda_seq[i])
      inc_seq_1[i] <- y_h
      inc_seq_2[i] <- y_l
      cons_seq_1[i] <- y_h - transfer_seq[i]
      cons_seq_2[i] <- y_l + transfer_seq[i]
    } else if (income_seq[i] == 3) {
      transfer_seq[i] <- (y_l - y_h * lambda_seq[i]) / (1 + lambda_seq[i])
      inc_seq_1[i] <- y_l
      inc_seq_2[i] <- y_h
      cons_seq_1[i] <- y_l - transfer_seq[i]
      cons_seq_2[i] <- y_h + transfer_seq[i]
    } else {
      transfer_seq[i] <- (y_l - y_l * lambda_seq[i]) / (1 + lambda_seq[i])
      inc_seq_1[i] <- y_l
      inc_seq_2[i] <- y_l
      cons_seq_1[i] <- y_l - transfer_seq[i]
      cons_seq_2[i] <- y_l + transfer_seq[i]
    }
  }
  plot_test <- ggplot() +
                geom_line(aes(delta_vals, log(lambda_lh_lower_vals), color="a")) + 
                geom_line(aes(delta_vals, log(lambda_lh_upper_vals), color="b")) + 
                geom_line(aes(delta_vals, log(lambda_hh_lower_vals), color="c")) + 
                geom_line(aes(delta_vals, log(lambda_hh_upper_vals), color="d")) +
                geom_line(aes(delta_vals, log(lambda_hl_lower_vals), color="e")) +
                geom_line(aes(delta_vals, log(lambda_hl_upper_vals), color="f")) +
                coord_cartesian(xlim = c(0.8, 1.0)) +

                geom_ribbon(aes(x = delta_vals,
                                ymin = log(lambda_lh_lower_vals),
                                ymax = log(lambda_lh_upper_vals)),
                                fill = "blue", alpha = "0.2") +
                geom_ribbon(aes(x = delta_vals,
                                ymin = log(lambda_hh_lower_vals),
                                ymax = log(lambda_hh_upper_vals)),
                                fill = "red", alpha = "0.2") +
                geom_ribbon(aes(x = delta_vals,
                                ymin = log(lambda_hl_lower_vals),
                                ymax = log(lambda_hl_upper_vals)),
                                fill = "green", alpha = "0.2") + 
                geom_vline(xintercept = delta, color = "black", size = 0.8, alpha = 0.6) +
                scale_color_manual(name = "End-points",
                                   values = c("blue", "purple", "brown", "red", "yellow", "green"), 
                                   labels = c("lambda_lh_lower",
                                              "lambda_lh_upper",
                                              "lambda_hh_lower",
                                              "lambda_hh_upper",
                                              "lambda_hl_lower",
                                              "lambda_hl_upper")) +
                xlab("Discount factor (delta)") +
                ylab("log of the ratio of marginal utilities (lambda)") +
                geom_point(aes(rep(delta, (n + 1)), c(log(lambda_init), log(lambda_seq)))) +
                geom_label_repel(
                  aes(rep(delta, (n + 1)), c(log(lambda_init), log(lambda_seq)), label = seq(0,n)),
                  box.padding = 0.35, point.padding = 0.5)
  
  result_table <- data.frame(cbind(seq(0,n), 
                                   round(c(log(lambda_init), log(lambda_seq)), 3), 
                                   c(NaN, income_realization_label_seq),
                                   round(c(NaN, transfer_seq), 3),
                                   round(c(NaN, cons_seq_1), 3),
                                   round(c(NaN, cons_seq_2), 3)
                                    ))
  colnames(result_table) <- c("Period",
                              "ln(lambda)",
                              "Income shocks",
                              "Net transfer (1 -> 2)",
                              "Consumption (1)",
                              "Consumption (2)")

  summary_table <- matrix(c(mean(cons_seq_1),
                            std(cons_seq_1),
                            mean(cons_seq_2),
                            std(cons_seq_2),
                            mean(inc_seq_1),
                            std(inc_seq_1),
                            mean(inc_seq_2),
                            std(inc_seq_2)),
                          byrow = TRUE, ncol = 2)
  summary_table <- round(summary_table, 3)
  colnames(summary_table) <- c("Mean", "Std")
  rownames(summary_table) <- c("Consumption (1)",
                               "Consumption (2)",
                               "Income (1)",
                               "Income (2)")
  return(list(plot_test, result_table, summary_table))
}
```

\newpage

### Autarky region

First, let me try the autarky region.

``` r
n <- 10
lambda_init <- exp(0)
delta <- 0.85
set.seed(1)
result <- lambda_change_plot(lambda_init, delta, n)
plot_figure <- result[1]
result_table <- result[2]
summary_table <- result[3]

kable(result_table, caption = "Autarky region")
```

<table class="kable_wrapper">
<caption>
Autarky region
</caption>
<tbody>
<tr>
<td>
| Period | ln(lambda) | Income shocks | Net transfer (1 -&gt; 2) | Consumption (1) | Consumption (2) |
|:-------|:-----------|:--------------|:-------------------------|:----------------|:----------------|
| 0      | 0          | NaN           | NaN                      | NaN             | NaN             |
| 1      | 0          | High, High    | 0                        | 1.333           | 1.333           |
| 2      | 0          | High, High    | 0                        | 1.333           | 1.333           |
| 3      | 0          | High, High    | 0                        | 1.333           | 1.333           |
| 4      | 0.693      | High, Low     | 0                        | 1.333           | 0.667           |
| 5      | 0          | High, High    | 0                        | 1.333           | 1.333           |
| 6      | -0.693     | Low, High     | 0                        | 0.667           | 1.333           |
| 7      | 0.693      | High, Low     | 0                        | 1.333           | 0.667           |
| 8      | 0          | High, High    | 0                        | 1.333           | 1.333           |
| 9      | 0          | High, High    | 0                        | 1.333           | 1.333           |
| 10     | 0          | High, High    | 0                        | 1.333           | 1.333           |

</td>
</tr>
</tbody>
</table>
``` r
plot_figure
```

    ## [[1]]

![](LTW_code_files/figure-markdown_github/unnamed-chunk-10-1.png)

There is no transfer between households and the ratio of marginal utilities fluctuates a lot in response to income shocks.

### Non-overlapping region

Now, I try and see the transition in the ratio of marginal utilities and consumptions in the non-overlapping region first.

``` r
n <- 10
lambda_init <- exp(0)
delta <- 0.9
set.seed(1)
result <- lambda_change_plot(lambda_init, delta, n)
plot_figure <- result[1]
result_table <- result[2]
summary_table <- result[3]

kable(result_table, caption = "Non-overlapping region")
```

<table class="kable_wrapper">
<caption>
Non-overlapping region
</caption>
<tbody>
<tr>
<td>
| Period | ln(lambda) | Income shocks | Net transfer (1 -&gt; 2) | Consumption (1) | Consumption (2) |
|:-------|:-----------|:--------------|:-------------------------|:----------------|:----------------|
| 0      | 0          | NaN           | NaN                      | NaN             | NaN             |
| 1      | 0          | High, High    | 0                        | 1.333           | 1.333           |
| 2      | 0          | High, High    | 0                        | 1.333           | 1.333           |
| 3      | 0          | High, High    | 0                        | 1.333           | 1.333           |
| 4      | 0.396      | High, Low     | 0.138                    | 1.195           | 0.805           |
| 5      | 0.033      | High, High    | -0.022                   | 1.355           | 1.311           |
| 6      | -0.396     | Low, High     | -0.138                   | 0.805           | 1.195           |
| 7      | 0.396      | High, Low     | 0.138                    | 1.195           | 0.805           |
| 8      | 0.033      | High, High    | -0.022                   | 1.355           | 1.311           |
| 9      | 0.033      | High, High    | -0.022                   | 1.355           | 1.311           |
| 10     | 0.033      | High, High    | -0.022                   | 1.355           | 1.311           |

</td>
</tr>
</tbody>
</table>
``` r
plot_figure
```

    ## [[1]]

![](LTW_code_files/figure-markdown_github/unnamed-chunk-11-1.png)

It is interesting to see the path dependence of transfers: after receiving a bad shock, a household keeps to \`\`repay the loan'' by transferring to the other household even when the realized incomes are the same, until the other household receives a bad shock and she herself receives a good shock simultaneously. Also, I look at the mean and standard deviations of incomes and consumptions:

``` r
kable(summary_table, caption = "Non-overlapping region (summary statistics)")
```

<table class="kable_wrapper">
<caption>
Non-overlapping region (summary statistics)
</caption>
<tbody>
<tr>
<td>
|                 |   Mean|    Std|
|-----------------|------:|------:|
| Consumption (1) |  1.262|  0.173|
| Consumption (2) |  1.205|  0.215|
| Income (1)      |  1.267|  0.211|
| Income (2)      |  1.200|  0.281|

</td>
</tr>
</tbody>
</table>
The table shows that while mean consumptions are similar to mean incomes, the fluctuation in consumptions are smaller than that of incomes. This is the benefit of risk sharing, although full insurance is not achieved since *λ* fluctuates over time.

### Overlapping region

Next, I look at how *λ* changes over time in the overlapping region:

``` r
n <- 10
lambda_init <- exp(0)
delta <- 0.95
set.seed(1)
result <- lambda_change_plot(lambda_init, delta, n)
plot_figure <- result[1]
result_table <- result[2]
summary_table <- result[3]

kable(result_table, caption = "Overlapping region")
```

<table class="kable_wrapper">
<caption>
Overlapping region
</caption>
<tbody>
<tr>
<td>
| Period | ln(lambda) | Income shocks | Net transfer (1 -&gt; 2) | Consumption (1) | Consumption (2) |
|:-------|:-----------|:--------------|:-------------------------|:----------------|:----------------|
| 0      | 0          | NaN           | NaN                      | NaN             | NaN             |
| 1      | 0          | High, High    | 0                        | 1.333           | 1.333           |
| 2      | 0          | High, High    | 0                        | 1.333           | 1.333           |
| 3      | 0          | High, High    | 0                        | 1.333           | 1.333           |
| 4      | 0.04       | High, Low     | 0.313                    | 1.02            | 0.98            |
| 5      | 0.04       | High, High    | -0.027                   | 1.36            | 1.307           |
| 6      | -0.04      | Low, High     | -0.313                   | 0.98            | 1.02            |
| 7      | 0.04       | High, Low     | 0.313                    | 1.02            | 0.98            |
| 8      | 0.04       | High, High    | -0.027                   | 1.36            | 1.307           |
| 9      | 0.04       | High, High    | -0.027                   | 1.36            | 1.307           |
| 10     | 0.04       | High, High    | -0.027                   | 1.36            | 1.307           |

</td>
</tr>
</tbody>
</table>
``` r
plot_figure
```

    ## [[1]]

![](LTW_code_files/figure-markdown_github/unnamed-chunk-13-1.png)

Compared to the previous case, the fluctuation in the ratio of marginal utilities seem smaller. This reuslts in less fluctuating consumptions as the table shows.

``` r
kable(summary_table, caption = "Overlapping region (summary statistics)")
```

<table class="kable_wrapper">
<caption>
Overlapping region (summary statistics)
</caption>
<tbody>
<tr>
<td>
|                 |   Mean|    Std|
|-----------------|------:|------:|
| Consumption (1) |  1.246|  0.166|
| Consumption (2) |  1.221|  0.158|
| Income (1)      |  1.267|  0.211|
| Income (2)      |  1.200|  0.281|

</td>
</tr>
</tbody>
</table>
### First-best region

Finally, I try to see the transitions in the first-best region.

``` r
n <- 10
lambda_init <- exp(0)
delta <- 0.98
set.seed(1)
result <- lambda_change_plot(lambda_init, delta, n)
plot_figure <- result[1]
result_table <- result[2]
summary_table <- result[3]

kable(result_table, caption = "First-best region")
```

<table class="kable_wrapper">
<caption>
First-best region
</caption>
<tbody>
<tr>
<td>
| Period | ln(lambda) | Income shocks | Net transfer (1 -&gt; 2) | Consumption (1) | Consumption (2) |
|:-------|:-----------|:--------------|:-------------------------|:----------------|:----------------|
| 0      | 0          | NaN           | NaN                      | NaN             | NaN             |
| 1      | 0          | High, High    | 0                        | 1.333           | 1.333           |
| 2      | 0          | High, High    | 0                        | 1.333           | 1.333           |
| 3      | 0          | High, High    | 0                        | 1.333           | 1.333           |
| 4      | 0          | High, Low     | 0.333                    | 1               | 1               |
| 5      | 0          | High, High    | 0                        | 1.333           | 1.333           |
| 6      | 0          | Low, High     | -0.333                   | 1               | 1               |
| 7      | 0          | High, Low     | 0.333                    | 1               | 1               |
| 8      | 0          | High, High    | 0                        | 1.333           | 1.333           |
| 9      | 0          | High, High    | 0                        | 1.333           | 1.333           |
| 10     | 0          | High, High    | 0                        | 1.333           | 1.333           |

</td>
</tr>
</tbody>
</table>
``` r
plot_figure
```

    ## [[1]]

![](LTW_code_files/figure-markdown_github/unnamed-chunk-15-1.png)

As the name of the region suggests, the ratio of marginal utilities stays at the initial position and never changes. This means that full insurance is achieved, although due to the aggregate shocks consumptions fluctuate as the table shows:

``` r
kable(summary_table, caption = "First-best region (summary statistics)")
```

<table class="kable_wrapper">
<caption>
First-best region (summary statistics)
</caption>
<tbody>
<tr>
<td>
|                 |   Mean|    Std|
|-----------------|------:|------:|
| Consumption (1) |  1.233|  0.161|
| Consumption (2) |  1.233|  0.161|
| Income (1)      |  1.267|  0.211|
| Income (2)      |  1.200|  0.281|

</td>
</tr>
</tbody>
</table>
### Comparison with static limited commitment model (Coate & Ravallion (1993))

What is another fun exercise? A comparison with static limited commitment model! In static limited commitment model, transfers do not depend on the history of states but only on the current state. Presumably this restricts the range of transfer contracts and thus weakens the performance of risk sharing. I try to see this numerically.

Again, I need to begin from derivation of end-points. The update rule is derived by Ligon, Thomas, and Worrall (2002) (equation (15)). This equation indicates the importance of initial *λ*. Thus, let me define this first:

``` r
lambda_0 <- 1
```

### Case 1: Autarky case

Under autarky, there is no transfer and thus the intervals become degenerate. Therefore, *λ*<sub>*h**h*, *l**o**w**e**r*</sub> = *λ*<sub>*h**h*, *u**p**p**e**r*</sub> = 1, *λ*<sub>*h**l*, *l**o**w**e**r*</sub> = *λ*<sub>*h**l*, *u**p**p**e**r*</sub> = 2, and *λ*<sub>*l**h*, *l**o**w**e**r*</sub> = *λ*<sub>*l**h*, *u**p**p**e**r*</sub> = 1/2.

### Case 2: Non-overlapping case

To calculate the end-points of the intervals, evaluate profits at the following end-points and the initial point:

-   *λ*<sub>*l**h*, *u**p**p**e**r*</sub>, at which the discounted surplus of 2 is 0 and the next-period ratio of marginal utilities are *λ*<sub>*l**h*, *u**p**p**e**r*</sub>, *λ*<sub>0</sub>, *λ*<sub>*h**l*, *l**o**w**e**r*</sub>, with probabilities *p*(1 − *p*), (1 − *p*)<sup>2</sup> + *p*<sup>2</sup>, and *p*(1 − *p*), respectively.
-   *λ*<sub>*h**h*, *l**o**w**e**r*</sub>, at which the discounted surplus of 1 is 0 and the next-period ratio of marginal utilities are *λ*<sub>*l**h*, *u**p**p**e**r*</sub>, *λ*<sub>0</sub>, *λ*<sub>*h**l*, *l**o**w**e**r*</sub>, with probabilities *p*(1 − *p*), (1 − *p*)<sup>2</sup> + *p*<sup>2</sup>, and *p*(1 − *p*), respectively.
-   *λ*<sub>*h**h*, *u**p**p**e**r*</sub>, at which the discounted surplus of 2 is 0 and the next-period ratio of marginal utilities are *λ*<sub>*l**h*, *u**p**p**e**r*</sub>, *λ*<sub>0</sub>, *λ*<sub>*h**l*, *l**o**w**e**r*</sub>, with probabilities *p*(1 − *p*), (1 − *p*)<sup>2</sup> + *p*<sup>2</sup>, and *p*(1 − *p*), respectively.
-   *λ*<sub>*h**l*, *l**o**w**e**r*</sub>, at which the discounted surplus of 1 is 0 and the next-period ratio of marginal utilities are *λ*<sub>*l**h*, *u**p**p**e**r*</sub>, *λ*<sub>0</sub>, *λ*<sub>*h**l*, *l**o**w**e**r*</sub>, with probabilities *p*(1 − *p*), (1 − *p*)<sup>2</sup> + *p*<sup>2</sup>, and *p*(1 − *p*), respectively.
-   *λ*<sub>0</sub>, at which the discounted surpluses of 1 and 2 are *U*<sub>0</sub> and *V*<sub>0</sub> and the next-period ratio of marginal utilities are *λ*<sub>*l**h*, *u**p**p**e**r*</sub>, *λ*<sub>0</sub>, *λ*<sub>*h**l*, *l**o**w**e**r*</sub>, with probabilities *p*(1 − *p*), (1 − *p*)<sup>2</sup> + *p*<sup>2</sup>, and *p*(1 − *p*), respectively.

The conditions are

$\\begin{aligned}  0 &= v(y\_2(h) + (y\_1(l) - y\_2(h) \\lambda\_{lh, upper}) / (1 + \\lambda\_{lh, upper})) - v(y\_2(h)) + \\\\  &\\delta ((1 - 2 \* p + 2 \* p^2) V\_0 + p (1 - p) V\_{hl, upper}) \\\\  U\_{lh, upper} &= u(y\_1(l) - (y\_1(l) - y\_2(h) \\lambda\_{lh, upper}) / (1 + \\lambda\_{lh, upper})) - u(y\_1(l)) + \\\\  &\\delta ((1 - 2 \* p + 2 \* p^2) U\_0 + p (1 - p) U\_{lh, upper}) \\\\  V\_{hh, upper} &= v(y\_2(h) + (y\_1(h) - y\_2(h) \\lambda\_{hh, lower}) / (1 + \\lambda\_{hh, lower})) - v(y\_2(h)) + \\\\  &\\delta ((1 - 2 \* p + 2 \* p^2) V\_0 + p (1 - p) V\_{hl, upper}) \\\\  0 &= u(y\_1(h) - (y\_1(h) - y\_2(h) \\lambda\_{hh, lower}) / (1 + \\lambda\_{hh, lower})) - u(y\_1(h)) + \\\\  &\\delta ((1 - 2 \* p + 2 \* p^2) U\_0 + p (1 - p) U\_{lh, upper}) \\\\  0 &= v(y\_2(h) + (y\_1(h) - y\_2(h) \\lambda\_{hh, upper}) / (1 + \\lambda\_{hh, upper})) - v(y\_2(h)) + \\\\  &\\delta ((1 - 2 \* p + 2 \* p^2) V\_0 + p (1 - p) V\_{hl, upper}) \\\\  U\_{hh, upper} &= u(y\_1(h) - (y\_1(h) - y\_2(h) \\lambda\_{hh, upper}) / (1 + \\lambda\_{hh, upper})) - u(y\_1(h)) + \\\\  &\\delta ((1 - 2 \* p + 2 \* p^2) U\_0 + p (1 - p) U\_{lh, upper}) \\\\  V\_{hl, upper} &= v(y\_2(l) + (y\_1(h) - y\_2(l) \\lambda\_{hl, lower}) / (1 + \\lambda\_{hl, lower})) - v(y\_2(l)) + \\\\  &\\delta ((1 - 2 \* p + 2 \* p^2) V\_0 + p (1 - p) V\_{hl, upper}) \\\\  0 &= u(y\_1(h) - (y\_1(h) - y\_2(l) \\lambda\_{hl, lower}) / (1 + \\lambda\_{hl, lower})) - u(y\_1(h)) + \\\\  &\\delta ((1 - 2 \* p + 2 \* p^2) U\_0 + p (1 - p) U\_{lh, upper}) \\\\  \\\\  V\_0 &= v(y\_2(h) + (y\_1(h) - y\_2(h) \\lambda\_0) / (1 + \\lambda\_0)) - v(y\_2(h)) + \\\\  &\\delta ((1 - 2 \* p + 2 \* p^2) V\_0 + p (1 - p) V\_{hl, upper}) \\\\  U\_0 &= u(y\_1(h) - (y\_1(h) - y\_2(h) \\lambda\_0) / (1 + \\lambda\_0)) - u(y\_1(h)) + \\\\  &\\delta ((1 - 2 \* p + 2 \* p^2) U\_0 + p (1 - p) U\_{lh, upper}) \\end{aligned}$

From the last two conditions, the discounted surpluses at the initial point are

$\\begin{aligned}  V\_0 &= \\frac{v(y\_2(h) + (y\_1(h) - y\_2(h) \\lambda\_0) / (1 + \\lambda\_0)) - v(y\_2(h)) + \\delta p (1 - p) V\_{hl, upper}}{1 - \\delta (1 - 2 \* p + 2 \* p^2)} \\\\  U\_0 &= \\frac{u(y\_1(h) - (y\_1(h) - y\_2(h) \\lambda\_0) / (1 + \\lambda\_0)) - u(y\_1(h)) + \\delta p (1 - p) U\_{lh, upper}}{1 - \\delta (1 - 2 \* p + 2 \* p^2)} \\end{aligned}$

The following function solves this system of equations:

``` r
# For the case where there is no overlap, solve end-points and surpluses at these points

no_overlap_slc <- function(x, delta, y_l, y_h, p){
  # x <- c(lambda_{lh, upper}, lambda_{hh, lower}, lambda_{hh, upper}, lambda_{hl, lower}, 
  #       U_{lh, upper}, V_{hh, upper}, U_{hh, upper}, V_{hl, upper})
  # First I derive the discounted surpluses at the interior of intervals
  V_0 <- (v(y_h + (y_h - y_h * lambda_0) / (1 + lambda_0)) - v(y_h) + 
          delta * p * (1 - p) * x[8]) / (1 - delta * (1 - 2 * p + 2 * p^2))
  U_0 <- (u(y_h - (y_h - y_h * lambda_0) / (1 + lambda_0)) - u(y_h) + 
          delta * p * (1 - p) * x[5]) / (1 - delta * (1 - 2 * p + 2 * p^2))

  y <- numeric(8)
  y[1] <- v(y_h + (y_l - y_h * x[1]) / (1 + x[1])) - v(y_h) + 
           delta * ((1 - 2 * p + 2 * p^2) * V_0 + p * (1 - p) * x[8])
  y[2] <- u(y_l - (y_l - y_h * x[1]) / (1 + x[1])) - u(y_l) + 
            delta * ((1 - 2 * p + 2 * p^2) * U_0 + p * (1 - p) * x[5]) - x[5]
  y[3] <- v(y_h + (y_h - y_h * x[2]) / (1 + x[2])) - v(y_h) + 
            delta * ((1 - 2 * p + 2 * p^2) * V_0 + p * (1 - p) * x[8]) - x[6]
  y[4] <- u(y_h - (y_h - y_h * x[2]) / (1 + x[2])) - u(y_h) + 
            delta * ((1 - 2 * p + 2 * p^2) * U_0 + p * (1 - p) * x[5])
  y[5] <- v(y_h + (y_h - y_h * x[3]) / (1 + x[3])) - v(y_h) + 
            delta * ((1 - 2 * p + 2 * p^2) * V_0 + p * (1 - p) * x[8])
  y[6] <- u(y_h - (y_h - y_h * x[3]) / (1 + x[3])) - u(y_h) + 
            delta * ((1 - 2 * p + 2 * p^2) * U_0 + p * (1 - p) * x[5]) - x[7]
  y[7] <- v(y_l + (y_h - y_l * x[4]) / (1 + x[4])) - v(y_l) + 
            delta * ((1 - 2 * p + 2 * p^2) * V_0 + p * (1 - p) * x[8]) - x[8]
  y[8] <- u(y_h - (y_h - y_l * x[4]) / (1 + x[4])) - u(y_h) + 
            delta * ((1 - 2 * p + 2 * p^2) * U_0 + p * (1 - p) * x[5])
  return(y)
}

x_start <- rep((y_l + y_h) / 2, 8)
no_overlap_solutions_slc <- nleqslv(x_start, 
                                no_overlap_slc, 
                                delta = 0.9,
                                y_l = y_l, 
                                y_h = y_h, 
                                p = p)$x
```

Given this function to solve each end-point, I obtain these end-points in this region, and the threshold with other regions (*δ*<sub>1, *s**l**c*</sub> (cutoff between autarky region and no-overlapping region) and *δ*<sub>2, *s**l**c*</sub> (cutoff between non-overlapping region and overlapping region)).

``` r
delta_vals <- seq(0.80, 0.999, by = 0.0001)
lambda_lh_upper_vals_no_overlap_slc <- numeric(length(delta_vals))
lambda_hh_lower_vals_no_overlap_slc <- numeric(length(delta_vals))
lambda_hh_upper_vals_no_overlap_slc <- numeric(length(delta_vals))
lambda_hl_lower_vals_no_overlap_slc <- numeric(length(delta_vals))
same_as_autarky_slc <- numeric(length(delta_vals))
x_start <- rep((y_l + y_h) / 2, 8)
for (delta in 1:length(delta_vals)){
  no_overlap_solutions_slc <- nleqslv(x_start, 
                                  no_overlap_slc, 
                                  delta = delta_vals[delta],
                                  y_l = y_l, 
                                  y_h = y_h, 
                                  p = p)$x
  lambda_lh_upper_vals_no_overlap_slc[delta] <- no_overlap_solutions_slc[1]
  lambda_hh_lower_vals_no_overlap_slc[delta] <- no_overlap_solutions_slc[2]
  lambda_hh_upper_vals_no_overlap_slc[delta] <- no_overlap_solutions_slc[3]
  lambda_hl_lower_vals_no_overlap_slc[delta] <- no_overlap_solutions_slc[4]
  same_as_autarky_slc[delta] <- is.logical(all.equal(no_overlap_solutions_slc[1], 
                                                     y_l / y_h, tolerance = 1e-3))
}

delta_1_slc <- delta_vals[max(which(same_as_autarky_slc == 1))]
delta_2_slc <- delta_vals[min(which(lambda_lh_upper_vals_no_overlap_slc
                                         >= lambda_hh_lower_vals_no_overlap_slc))]
```

### Case 3: Overlapping case

To calculate the end-points of the intervals, evaluate profits at the following end-points:

-   *λ*<sub>*l**h*, *u**p**p**e**r*</sub>, at which the discounted surplus of 2 is 0 and the next-period ratio of marginal utilities are *λ*<sub>*l**h*, *u**p**p**e**r*</sub>, *λ*<sub>0</sub>, and *λ*<sub>*h**l*, *l**o**w**e**r*</sub> with probabilities *p*(1 − *p*), (1 − *p*)<sup>2</sup> + *p*<sup>2</sup>, and *p*(1 − *p*), respectively.
-   *λ*<sub>*h**h*, *l**o**w**e**r*</sub>, at which the discounted surplus of 1 is 0 and the next-period ratio of marginal utilities are *λ*<sub>*l**h*, *l**o**w**e**r*</sub>, *λ*<sub>0</sub>, and *λ*<sub>*h**l*, *l**o**w**e**r*</sub> with probabilities *p*(1 − *p*), (1 − *p*)<sup>2</sup> + *p*<sup>2</sup>, and *p*(1 − *p*), respectively.
-   *λ*<sub>*h**h*, *u**p**p**e**r*</sub>, at which the discounted surplus of 2 is 0 and the next-period ratio of marginal utilities are *λ*<sub>*l**h*, *u**p**p**e**r*</sub>, *λ*<sub>0</sub>, and *λ*<sub>*l**l*, *u**p**p**e**r*</sub> with probabilities *p*(1 − *p*), (1 − *p*)<sup>2</sup> + *p*<sup>2</sup>, and *p*(1 − *p*), respectively.
-   *λ*<sub>*h**l*, *l**o**w**e**r*</sub>, at which the discounted surplus of 1 is 0 and the next-period ratio of marginal utilities are *λ*<sub>*l**h*, *u**p**p**e**r*</sub>, *λ*<sub>0</sub>, and *λ*<sub>*h**l*, *l**o**w**e**r*</sub> with probabilities *p*(1 − *p*), (1 − *p*)<sup>2</sup> + *p*<sup>2</sup> and *p*(1 − *p*), respectively.
-   *λ*<sub>0</sub>, at which the discounted surpluses of 1 and 2 are *U*<sub>0</sub> and *V*<sub>0</sub> and the next-period ratio of marginal utilities are *λ*<sub>*l**h*, *u**p**p**e**r*</sub>, *λ*<sub>0</sub>, *λ*<sub>*h**l*, *l**o**w**e**r*</sub>, with probabilities *p*(1 − *p*), (1 − *p*)<sup>2</sup> + *p*<sup>2</sup>, and *p*(1 − *p*), respectively.

Indeed they are exactly identical to the previous case, which is because of the lack of path dependence in the current model.

The conditions are

$\\begin{aligned}  0 &= v(y\_2(h) + (y\_1(l) - y\_2(h) \\lambda\_{lh, upper}) / (1 + \\lambda\_{lh, upper})) - v(y\_2(h)) + \\\\  &\\delta ((1 - 2 \* p + 2 \* p^2) V\_0 + p (1 - p) V\_{hl, upper}) \\\\  U\_{lh, upper} &= u(y\_1(l) - (y\_1(l) - y\_2(h) \\lambda\_{lh, upper}) / (1 + \\lambda\_{lh, upper})) - u(y\_1(l)) + \\\\  &\\delta ((1 - 2 \* p + 2 \* p^2) U\_0 + p (1 - p) U\_{lh, upper}) \\\\  V\_{hh, upper} &= v(y\_2(h) + (y\_1(h) - y\_2(h) \\lambda\_{hh, lower}) / (1 + \\lambda\_{hh, lower})) - v(y\_2(h)) + \\\\  &\\delta ((1 - 2 \* p + 2 \* p^2) V\_0 + p (1 - p) V\_{hl, upper}) \\\\  0 &= u(y\_1(h) - (y\_1(h) - y\_2(h) \\lambda\_{hh, lower}) / (1 + \\lambda\_{hh, lower})) - u(y\_1(h)) + \\\\  &\\delta ((1 - 2 \* p + 2 \* p^2) U\_0 + p (1 - p) U\_{lh, upper}) \\\\  0 &= v(y\_2(h) + (y\_1(h) - y\_2(h) \\lambda\_{hh, upper}) / (1 + \\lambda\_{hh, upper})) - v(y\_2(h)) + \\\\  &\\delta ((1 - 2 \* p + 2 \* p^2) V\_0 + p (1 - p) V\_{hl, upper}) \\\\  U\_{hh, upper} &= u(y\_1(h) - (y\_1(h) - y\_2(h) \\lambda\_{hh, upper}) / (1 + \\lambda\_{hh, upper})) - u(y\_1(h)) + \\\\  &\\delta ((1 - 2 \* p + 2 \* p^2) U\_0 + p (1 - p) U\_{lh, upper}) \\\\  V\_{hl, upper} &= v(y\_2(l) + (y\_1(h) - y\_2(l) \\lambda\_{hl, lower}) / (1 + \\lambda\_{hl, lower})) - v(y\_2(l)) + \\\\  &\\delta ((1 - 2 \* p + 2 \* p^2) V\_0 + p (1 - p) V\_{hl, upper}) \\\\  0 &= u(y\_1(h) - (y\_1(h) - y\_2(l) \\lambda\_{hl, lower}) / (1 + \\lambda\_{hl, lower})) - u(y\_1(h)) + \\\\  &\\delta ((1 - 2 \* p + 2 \* p^2) U\_0 + p (1 - p) U\_{lh, upper}) \\\\  \\\\  V\_0 &= v(y\_2(h) + (y\_1(h) - y\_2(h) \\lambda\_0) / (1 + \\lambda\_0)) - v(y\_2(h)) + \\\\  &\\delta ((1 - 2 \* p + 2 \* p^2) V\_0 + p (1 - p) V\_{hl, upper}) \\\\  U\_0 &= u(y\_1(h) - (y\_1(h) - y\_2(h) \\lambda\_0) / (1 + \\lambda\_0)) - u(y\_1(h)) + \\\\  &\\delta ((1 - 2 \* p + 2 \* p^2) U\_0 + p (1 - p) U\_{lh, upper}) \\end{aligned}$

From the last two conditions, the discounted surpluses at the initial point are

$\\begin{aligned}  V\_0 &= \\frac{v(y\_2(h) + (y\_1(h) - y\_2(h) \\lambda\_0) / (1 + \\lambda\_0)) - v(y\_2(h)) + \\delta p (1 - p) V\_{hl, upper}}{1 - \\delta (1 - 2 \* p + 2 \* p^2)} \\\\  U\_0 &= \\frac{u(y\_1(h) - (y\_1(h) - y\_2(h) \\lambda\_0) / (1 + \\lambda\_0)) - u(y\_1(h)) + \\delta p (1 - p) U\_{lh, upper}}{1 - \\delta (1 - 2 \* p + 2 \* p^2)} \\end{aligned}$

The following function solves this system of equations:

``` r
# For the case where there is no overlap, solve end-points and surpluses at these points

overlap_slc <- function(x, delta, y_l, y_h, p){
  # x <- c(lambda_{lh, upper}, lambda_{hh, lower}, lambda_{hh, upper}, lambda_{hl, lower}, 
  #       U_{lh, upper}, V_{hh, upper}, U_{hh, upper}, V_{hl, upper})
  # First I derive the discounted surpluses at the interior of intervals
  V_0 <- (v(y_h + (y_h - y_h * lambda_0) / (1 + lambda_0)) - v(y_h) + 
          delta * p * (1 - p) * x[8]) / (1 - delta * (1 - 2 * p + 2 * p^2))
  U_0 <- (u(y_h - (y_h - y_h * lambda_0) / (1 + lambda_0)) - u(y_h) + 
          delta * p * (1 - p) * x[5]) / (1 - delta * (1 - 2 * p + 2 * p^2))

  y <- numeric(8)
  y[1] <- v(y_h + (y_l - y_h * x[1]) / (1 + x[1])) - v(y_h) + 
          delta * ((1 - 2 * p + 2 * p^2) * V_0 + p * (1 - p) * x[8])
  y[2] <- u(y_l - (y_l - y_h * x[1]) / (1 + x[1])) - u(y_l) + 
          delta * ((1 - 2 * p + 2 * p^2) * U_0 + p * (1 - p) * x[5]) - x[5]
  y[3] <- v(y_h + (y_h - y_h * x[2]) / (1 + x[2])) - v(y_h) + 
          delta * ((1 - 2 * p + 2 * p^2) * V_0 + p * (1 - p) * x[8]) - x[6]
  y[4] <- u(y_h - (y_h - y_h * x[2]) / (1 + x[2])) - u(y_h) + 
          delta * ((1 - 2 * p + 2 * p^2) * U_0 + p * (1 - p) * x[5])
  y[5] <- v(y_h + (y_h - y_h * x[3]) / (1 + x[3])) - v(y_h) + 
          delta * ((1 - 2 * p + 2 * p^2) * V_0 + p * (1 - p) * x[8])
  y[6] <- u(y_h - (y_h - y_h * x[3]) / (1 + x[3])) - u(y_h) + 
          delta * ((1 - 2 * p + 2 * p^2) * U_0 + p * (1 - p) * x[5]) - x[7]
  y[7] <- v(y_l + (y_h - y_l * x[4]) / (1 + x[4])) - v(y_l) + 
          delta * ((1 - 2 * p + 2 * p^2) * V_0 + p * (1 - p) * x[8]) - x[8]
  y[8] <- u(y_h - (y_h - y_l * x[4]) / (1 + x[4])) - u(y_h) + 
          delta * ((1 - 2 * p + 2 * p^2) * U_0 + p * (1 - p) * x[5])
  return(y)
}

x_start <- rep((y_l + y_h) / 2, 8)
overlap_solutions_slc <- nleqslv(x_start, 
                                overlap_slc, 
                                delta = 0.9,
                                y_l = y_l, 
                                y_h = y_h, 
                                p = p)$x
```

Given this function to solve each end-point, I obtain these end-points in this region, and the cutoff between overlapping region and first-best region (*δ*<sub>3, *s**l**c*</sub>).

``` r
delta_vals <- seq(0.80, 0.999, by = 0.0001)
lambda_lh_upper_vals_overlap_slc <- numeric(length(delta_vals))
lambda_hh_lower_vals_overlap_slc <- numeric(length(delta_vals))
lambda_hh_upper_vals_overlap_slc <- numeric(length(delta_vals))
lambda_hl_lower_vals_overlap_slc <- numeric(length(delta_vals))
x_start <- rep((y_l + y_h) / 2, 8)
for (delta in 1:length(delta_vals)){
  overlap_solutions_slc <- nleqslv(x_start, 
                               overlap_slc, 
                               delta = delta_vals[delta],
                               y_l = y_l, 
                               y_h = y_h, 
                               p = p)$x
  lambda_lh_upper_vals_overlap_slc[delta] <- overlap_solutions_slc[1]
  lambda_hh_lower_vals_overlap_slc[delta] <- overlap_solutions_slc[2]
  lambda_hh_upper_vals_overlap_slc[delta] <- overlap_solutions_slc[3]
  lambda_hl_lower_vals_overlap_slc[delta] <- overlap_solutions_slc[4]
}

delta_3_slc <- delta_vals[min(which(lambda_lh_upper_vals_overlap_slc
                                         >= lambda_hl_lower_vals_overlap_slc))]
```

### Case 4: First-best case

To calculate the end-points of the intervals, evaluate profits at the following end-points and the initial point:

-   *λ*<sub>*l**h*, *u**p**p**e**r*</sub>, at which the discounted surplus of 2 is 0 and the next-period ratio of marginal utility is *λ*<sub>0</sub> with probability 1.
-   *λ*<sub>*h**h*, *l**o**w**e**r*</sub>, at which the discounted surplus of 1 is 0 and the next-period ratio of marginal utility is *λ*<sub>0</sub> and *λ*<sub>*h**l*, *l**o**w**e**r*</sub> with probabilities 1 − *p* + *p*<sup>2</sup> and *p* − *p*<sup>2</sup>, respectively.
-   *λ*<sub>*h**h*, *u**p**p**e**r*</sub>, at which the discounted surplus of 2 is 0 and the next-period ratio of marginal utility is *λ*<sub>0</sub> and *λ*<sub>*l**h*, *u**p**p**e**r*</sub> with probabilities 1 − *p* + *p*<sup>2</sup> and *p* − *p*<sup>2</sup>, respectively.
-   *λ*<sub>*h**l*, *l**o**w**e**r*</sub>, at which the discounted surplus of 1 is 0 and the next-period ratio of marginal utility is *λ*<sub>0</sub> with probability 1.
-   *λ*<sub>0</sub>, at which the discounted surpluses of 1 and 2 are *U*<sub>0</sub> and *V*<sub>0</sub> and the next-period ratio of marginal utility is *λ*<sub>0</sub> with probability 1.

Here, I define the discounted surpluses at the interior of intervals as the followings:

-   *U*<sub>0</sub><sup>\*</sup>, *V*<sub>0</sub><sup>\*</sup>: discounted surpluses when income realization is {Low, High} and the ratio of marginal utilities is *λ*<sub>0</sub>,
-   *U*<sub>0</sub><sup>\*\*</sup>, *V*<sub>0</sub><sup>\*\*</sup>: discounted surpluses when income realization is {High, Low} and the ratio of marginal utilities is *λ*<sub>0</sub>,

Then, the conditions are

$\\begin{aligned}  0 &= v(y\_2(h) + (y\_1(l) - y\_2(h) \\lambda\_{lh, upper}) / (1 + \\lambda\_{lh, upper})) - v(y\_2(h)) + \\\\  &\\delta ((p - p^2) V\_0^\* + (1 - 2 p + 2 p^2) V\_0 + (p - p^2) V\_0^{\*\*}) \\\\  U\_{lh, upper} &= u(y\_1(l) - (y\_1(l) - y\_2(h) \\lambda\_{lh, upper}) / (1 + \\lambda\_{lh, upper})) - u(y\_1(l)) + \\\\  &\\delta ((p - p^2) U\_0^\* + (1 - 2 p + 2 p^2) U\_0 + (p - p^2) U\_0^{\*\*}) \\\\  V\_{hh, upper} &= v(y\_2(h) + (y\_1(h) - y\_2(h) \\lambda\_{hh, lower}) / (1 + \\lambda\_{hh, lower})) - v(y\_2(h)) + \\\\  &\\delta ((p - p^2) V\_0^\* + (1 - 2 p + 2 p^2) V\_0 + (p - p^2) V\_0^{\*\*}) \\\\  0 &= u(y\_1(h) - (y\_1(h) - y\_2(h) \\lambda\_{hh, lower}) / (1 + \\lambda\_{hh, lower})) - u(y\_1(h)) + \\\\  &\\delta ((p - p^2) U\_0^\* + (1 - 2 p + 2 p^2) U\_0 + (p - p^2) U\_0^{\*\*}) \\\\  0 &= v(y\_2(h) + (y\_1(h) - y\_2(h) \\lambda\_{hh, upper}) / (1 + \\lambda\_{hh, upper})) - v(y\_2(h)) + \\\\  &\\delta ((p - p^2) V\_0^\* + (1 - 2 p + 2 p^2) V\_0 + (p - p^2) V\_0^{\*\*}) \\\\  U\_{hh, upper} &= u(y\_1(h) - (y\_1(h) - y\_2(h) \\lambda\_{hh, upper}) / (1 + \\lambda\_{hh, upper})) - u(y\_1(h)) + \\\\  &\\delta ((p - p^2) U\_0^\* + (1 - 2 p + 2 p^2) U\_0 + (p - p^2) U\_0^{\*\*}) \\\\  V\_{hl, upper} &= v(y\_2(l) + (y\_1(h) - y\_2(l) \\lambda\_{hl, lower}) / (1 + \\lambda\_{hl, lower})) - v(y\_2(l)) + \\\\  &\\delta ((p - p^2) V\_0^\* + (1 - 2 p + 2 p^2) V\_0 + (p - p^2) V\_0^{\*\*}) \\\\  0 &= u(y\_1(h) - (y\_1(h) - y\_2(l) \\lambda\_{hl, lower}) / (1 + \\lambda\_{hl, lower})) - u(y\_1(h)) + \\\\  &\\delta ((p - p^2) U\_0^\* + (1 - 2 p + 2 p^2) U\_0 + (p - p^2) U\_0^{\*\*}) \\\\  \\\\  V\_0^\* &= v(y\_2(l) + (y\_1(h) - y\_2(l) \\lambda\_0) / (1 + \\lambda\_0)) - v(y\_2(l)) + \\\\  &\\delta ((p - p^2) V\_0^\* + (1 - 2 p + 2 p^2) V\_0 + (p - p^2) V\_0^{\*\*}) \\\\  U\_0^\* &= u(y\_1(h) - (y\_1(h) - y\_2(l) \\lambda\_0) / (1 + \\lambda\_0)) - u(y\_1(h)) + \\\\  &\\delta ((p - p^2) U\_0^\* + (1 - 2 p + 2 p^2) U\_0 + (p - p^2) U\_0^{\*\*}) \\\\  V\_0^{\*\*} &= v(y\_2(h) + (y\_1(l) - y\_2(h) \\lambda\_0) / (1 + \\lambda\_0)) - v(y\_2(h)) + \\\\  &\\delta ((p - p^2) V\_0^\* + (1 - 2 p + 2 p^2) V\_0 + (p - p^2) V\_0^{\*\*}) \\\\  U\_0^{\*\*} &= u(y\_1(l) - (y\_1(l) - y\_2(h) \\lambda\_0) / (1 + \\lambda\_0)) - u(y\_1(l)) + \\\\  &\\delta ((p - p^2) U\_0^\* + (1 - 2 p + 2 p^2) U\_0 + (p - p^2) U\_0^{\*\*}) \\\\  V\_0 &= v(y\_2(h) + (y\_1(h) - y\_2(h) \\lambda\_0) / (1 + \\lambda\_0)) - v(y\_2(h)) + \\\\  &\\delta ((p - p^2) V\_0^\* + (1 - 2 p + 2 p^2) V\_0 + (p - p^2) V\_0^{\*\*}) \\\\  U\_0 &= u(y\_1(h) - (y\_1(h) - y\_2(h) \\lambda\_0) / (1 + \\lambda\_0)) - u(y\_1(h)) + \\\\  &\\delta ((p - p^2) U\_0^\* + (1 - 2 p + 2 p^2) U\_0 + (p - p^2) U\_0^{\*\*}) \\end{aligned}$

From these we can derive the followings:

$\\begin{aligned}  V\_0^\* &= v(y\_2(l) + (y\_1(h) - y\_2(l) \\lambda\_0) / (1 + \\lambda\_0)) - v(y\_2(l)) - \\\\  &(v(y\_2(h) + (y\_1(l) - y\_2(h) \\lambda\_{lh, upper}) / (1 + \\lambda\_{lh, upper})) - v(y\_2(h))) \\\\  U\_0^\* &= u(y\_1(h) - (y\_1(h) - y\_2(l) \\lambda\_0) / (1 + \\lambda\_0)) - u(y\_1(h)) + U\_{lh, upper} - \\\\  &(u(y\_1(l) - (y\_1(l) - y\_2(h) \\lambda\_{lh, upper}) / (1 + \\lambda\_{lh, upper})) - u(y\_1(l))) \\\\  V\_0^{\*\*} &= v(y\_2(h) + (y\_1(l) - y\_2(h) \\lambda\_0) / (1 + \\lambda\_0)) - v(y\_2(h)) - \\\\  &(v(y\_2(h) + (y\_1(l) - y\_2(h) \\lambda\_{lh, upper}) / (1 + \\lambda\_{lh, upper})) - v(y\_2(h))) \\\\  U\_0^{\*\*} &= u(y\_1(l) - (y\_1(l) - y\_2(h) \\lambda\_0) / (1 + \\lambda\_0)) - u(y\_1(l)) + U\_{lh, upper} - \\\\  &(u(y\_1(l) - (y\_1(l) - y\_2(h) \\lambda\_{lh, upper}) / (1 + \\lambda\_{lh, upper})) - u(y\_1(l))) \\\\  V\_0 &= \\frac{v(y\_2(h) + (y\_1(h) - y\_2(h) \\lambda\_0) / (1 + \\lambda\_0)) - v(y\_2(h)) + \\delta ((p - p^2) V\_0^\* + (p - p^2) V\_0^{\*\*})}{1 - \\delta (1 - 2 p + 2 p^2)} \\\\  U\_0 &= \\frac{u(y\_1(h) - (y\_1(h) - y\_2(h) \\lambda\_0) / (1 + \\lambda\_0)) - u(y\_1(h)) + \\delta ((p - p^2) U\_0^\* + (p - p^2) U\_0^{\*\*})}{1 - \\delta (1 - 2 p + 2 p^2)} \\end{aligned}$

The following function solves the system of equations:

``` r
# For the case where first best can be achieved, solve end-points and surpluses at these points

first_best_slc <- function(x, delta, y_l, y_h, p){
  # x <- c(lambda_{lh, upper}, lambda_{hh, lower}, lambda_{hh, upper}, lambda_{hl, lower}, 
  #       U_{lh, upper}, V_{hh, upper}, U_{hh, upper}, V_{hl, upper})
  # First I derive the discounted surpluses at the interior of intervals and at the initial point

  V_0_star  <-  v(y_h + (y_l - y_h * lambda_0) / (1 + lambda_0)) - v(y_h) - 
                (v(y_h + (y_l - y_h * x[1]) / (1 + x[1])) - v(y_h))
  U_0_star  <-  u(y_l - (y_l - y_h * lambda_0) / (1 + lambda_0)) - u(y_l) + 
                x[5] - (u(y_l - (y_l - y_h * x[1]) / (1 + x[1])) - u(y_l))
  V_0_2star <-  v(y_l + (y_h - y_l * lambda_0) / (1 + lambda_0)) - v(y_l) -
                (v(y_h + (y_l - y_h * x[1]) / (1 + x[1])) - v(y_h))
  U_0_2star <-  u(y_h - (y_h - y_l * lambda_0) / (1 + lambda_0)) - u(y_h) + 
                x[5] - (u(y_l - (y_l - y_h * x[1]) / (1 + x[1])) - u(y_l))
  V_0       <- (v(y_h + (y_h - y_h * lambda_0) / (1 + lambda_0)) - v(y_h) + 
                delta * ((p - p^2) * V_0_star + (p - p^2) * V_0_2star)) / 
                (1 - delta * (1 - 2 * p + 2 * p^2))
  U_0       <- (u(y_h - (y_h - y_h * lambda_0) / (1 + lambda_0)) - u(y_h) + 
                delta * ((p - p^2) * U_0_star + (p - p^2) * U_0_2star)) / 
                (1 - delta * (1 - 2 * p + 2 * p^2))

  y <- numeric(8)
  y[1] <- v(y_h + (y_l - y_h * x[1]) / (1 + x[1])) - v(y_h) + 
          delta * ((p - p^2) * V_0_star + (1 - 2 * p + 2 * p^2) * V_0 + (p - p^2) * V_0_2star)
  y[2] <- u(y_l - (y_l - y_h * x[1]) / (1 + x[1])) - u(y_l) + 
          delta * ((p - p^2) * U_0_star + (1 - 2 * p + 2 * p^2) * U_0 + (p - p^2) * U_0_2star) - x[5]
  y[3] <- v(y_h + (y_h - y_h * x[2]) / (1 + x[2])) - v(y_h) + 
          delta * ((p - p^2) * V_0_star + (1 - 2 * p + 2 * p^2) * V_0 + (p - p^2) * V_0_2star) - x[6]
  y[4] <- u(y_h - (y_h - y_h * x[2]) / (1 + x[2])) - u(y_h) + 
          delta * ((p - p^2) * U_0_star + (1 - 2 * p + 2 * p^2) * U_0 + (p - p^2) * U_0_2star)
  y[5] <- v(y_h + (y_h - y_h * x[3]) / (1 + x[3])) - v(y_h) + 
          delta * ((p - p^2) * V_0_star + (1 - 2 * p + 2 * p^2) * V_0 + (p - p^2) * V_0_2star)
  y[6] <- u(y_h - (y_h - y_h * x[3]) / (1 + x[3])) - u(y_h) + 
          delta * ((p - p^2) * U_0_star + (1 - 2 * p + 2 * p^2) * U_0 + (p - p^2) * U_0_2star) - x[7]
  y[7] <- v(y_l + (y_h - y_l * x[4]) / (1 + x[4])) - v(y_l) + 
          delta * ((p - p^2) * V_0_star + (1 - 2 * p + 2 * p^2) * V_0 + (p - p^2) * V_0_2star) - x[8]
  y[8] <- u(y_h - (y_h - y_l * x[4]) / (1 + x[4])) - u(y_h) + 
          delta * ((p - p^2) * U_0_star + (1 - 2 * p + 2 * p^2) * U_0 + (p - p^2) * U_0_2star)
  return(y)
}
```

Given this function to solve each end-point, I obtain these end-points in this region:

``` r
delta_vals <- seq(0.80, 0.999, by = 0.0001)
lambda_lh_upper_vals_first_best_slc <- numeric(length(delta_vals))
lambda_hh_lower_vals_first_best_slc <- numeric(length(delta_vals))
lambda_hh_upper_vals_first_best_slc <- numeric(length(delta_vals))
lambda_hl_lower_vals_first_best_slc <- numeric(length(delta_vals))
x_start <- rep((y_l + y_h) / 2, 8)
for (delta in 1:length(delta_vals)){
  first_best_solutions_slc <- nleqslv(x_start, 
                               first_best_slc, 
                               delta = delta_vals[delta],
                               y_l = y_l, 
                               y_h = y_h, 
                               p = p)$x
  lambda_lh_upper_vals_first_best_slc[delta] <- first_best_solutions_slc[1]
  lambda_hh_lower_vals_first_best_slc[delta] <- first_best_solutions_slc[2]
  lambda_hh_upper_vals_first_best_slc[delta] <- first_best_solutions_slc[3]
  lambda_hl_lower_vals_first_best_slc[delta] <- first_best_solutions_slc[4]
}
```

### Plot the lambda's in static model

With the information at hand, I plot the end-points in static model:

``` r
delta_vals <- seq(0.80, 0.999, by = 0.0001)

lambda_lh_upper_vals_slc <- numeric(length(delta_vals))
lambda_hh_lower_vals_slc <- numeric(length(delta_vals))
lambda_hh_upper_vals_slc <- numeric(length(delta_vals))
lambda_hl_lower_vals_slc <- numeric(length(delta_vals))
x_start <- rep((y_l + y_h) / 2, 8)
for (delta in 1:length(delta_vals)){
  if (delta_vals[delta] <= delta_1_slc){
    lambda_lh_upper_vals_slc[delta] <- y_l / y_h 
    lambda_hh_lower_vals_slc[delta] <- 1
    lambda_hh_upper_vals_slc[delta] <- 1
    lambda_hl_lower_vals_slc[delta] <- y_h / y_l 
  } else if (delta_vals[delta] > delta_1_slc & delta_vals[delta] <= delta_2_slc){
    lambda_lh_upper_vals_slc[delta] <- lambda_lh_upper_vals_no_overlap_slc[delta]
    lambda_hh_lower_vals_slc[delta] <- lambda_hh_lower_vals_no_overlap_slc[delta]
    lambda_hh_upper_vals_slc[delta] <- lambda_hh_upper_vals_no_overlap_slc[delta]
    lambda_hl_lower_vals_slc[delta] <- lambda_hl_lower_vals_no_overlap_slc[delta]
  } else if (delta_vals[delta] > delta_2_slc & delta_vals[delta] <= delta_3_slc){
    lambda_lh_upper_vals_slc[delta] <- lambda_lh_upper_vals_overlap_slc[delta]
    lambda_hh_lower_vals_slc[delta] <- lambda_hh_lower_vals_overlap_slc[delta]
    lambda_hh_upper_vals_slc[delta] <- lambda_hh_upper_vals_overlap_slc[delta]
    lambda_hl_lower_vals_slc[delta] <- lambda_hl_lower_vals_overlap_slc[delta]
  } else {
    lambda_lh_upper_vals_slc[delta] <- lambda_lh_upper_vals_first_best_slc[delta]
    lambda_hh_lower_vals_slc[delta] <- lambda_hh_lower_vals_first_best_slc[delta]
    lambda_hh_upper_vals_slc[delta] <- lambda_hh_upper_vals_first_best_slc[delta]
    lambda_hl_lower_vals_slc[delta] <- lambda_hl_lower_vals_first_best_slc[delta]
  }
}

lambda_lh_lower_vals_slc <- rep(y_l / y_h, length(delta_vals))
lambda_hl_upper_vals_slc <- rep(y_h / y_l, length(delta_vals))

ggplot() +
  geom_line(aes(delta_vals, log(lambda_lh_lower_vals_slc), color="a")) + 
  geom_line(aes(delta_vals, log(lambda_lh_upper_vals_slc), color="b")) + 
  geom_line(aes(delta_vals, log(lambda_hh_lower_vals_slc), color="c")) + 
  geom_line(aes(delta_vals, log(lambda_hh_upper_vals_slc), color="d")) +
  geom_line(aes(delta_vals, log(lambda_hl_lower_vals_slc), color="e")) +
  geom_line(aes(delta_vals, log(lambda_hl_upper_vals_slc), color="f")) +
  coord_cartesian(xlim = c(0.8, 1.0), ylim = c(log(y_l / y_h), log(y_h / y_l))) +
  geom_ribbon(aes(x = delta_vals,
                  ymin = log(lambda_lh_lower_vals_slc),
                  ymax = log(lambda_lh_upper_vals_slc)),
                  fill = "blue", alpha = "0.2") +
  geom_ribbon(aes(x = delta_vals,
                  ymin = log(lambda_hh_lower_vals_slc),
                  ymax = log(lambda_hh_upper_vals_slc)),
                  fill = "red", alpha = "0.2") +
  geom_ribbon(aes(x = delta_vals,
                  ymin = log(lambda_hl_lower_vals_slc),
                  ymax = log(lambda_hl_upper_vals_slc)),
                  fill = "green", alpha = "0.2") + 
  scale_color_manual(name = "End-points", values = c("blue",
                                                     "purple",
                                                     "brown",
                                                     "red",
                                                     "yellow",
                                                     "green"), 
                     labels = c("lambda_lh_lower",
                                "lambda_lh_upper",
                                "lambda_hh_lower",
                                "lambda_hh_upper",
                                "lambda_hl_lower",
                                "lambda_hl_upper")) +
  xlab("Discount factor (delta)") +
  ylab("log of the ratio of marginal utilities (lambda)")
```

![](LTW_code_files/figure-markdown_github/unnamed-chunk-24-1.png)

The cutoffs of the regions are *δ*<sub>1, *s**l**c*</sub> = 0.9175, *δ*<sub>2, *s**l**c*</sub> = 0.9464, and *δ*<sub>3, *s**l**c*</sub> = 0.9645. Note that the cutoffs of the regions in the dynamic limited commitment model are *δ*<sub>1</sub> = 0.8666, *δ*<sub>2</sub> = 0.9351, and *δ*<sub>3</sub> = 0.9645.

### Plot the figure

I will next see how *λ* changes over time in response to income shocks. This is the function to create the figure:

``` r
lambda_change_plot_slc <- function(lambda_init, delta, n){


  # Sequence of income shocks
  income_realization <- c(1, 2, 3, 4)
  income_realization_label <- c("High, High", "High, Low", "Low, High", "Low, Low")
  income_seq <- sample(income_realization,
                       size = n,
                       replace = TRUE,
                       prob = c((1 - p)^2,
                                p * (1 - p),
                                p * (1 - p),
                                p^2))

  # Sequence of the ratios of marginal utilities under autarky 
  xi_seq <- numeric(n)
  income_realization_label_seq <- rep("", n)
  for (i in 1:n){
    xi_seq[i] <- 1 * (income_seq[i] == 1) + 
                (y_h / y_l) * (income_seq[i] == 2) +
                (y_l / y_h) * (income_seq[i] == 3) +
                1 * (income_seq[i] == 4)
    income_realization_label_seq[i] <- income_realization_label[income_seq[i]]
  }
  
  # Given the initial lambda, derive the sequence of lambda's
  
  if (delta <= delta_1_slc){
    lambda_lh_upper <- y_l / y_h
    lambda_hh_lower <- 1
    lambda_hh_upper <- 1
    lambda_hl_lower <- y_h / y_l
  } else if (delta > delta_1_slc & delta <= delta_2_slc){
    no_overlap_solutions_slc <- nleqslv(x_start, 
                                    no_overlap_slc, 
                                    delta = delta,
                                    y_l = y_l, 
                                    y_h = y_h, 
                                    p = p)$x
    lambda_lh_upper <- no_overlap_solutions_slc[1]
    lambda_hh_lower <- no_overlap_solutions_slc[2]
    lambda_hh_upper <- no_overlap_solutions_slc[3]
    lambda_hl_lower <- no_overlap_solutions_slc[4]
  } else if (delta > delta_2_slc & delta <= delta_3_slc) {
    overlap_solutions_slc <- nleqslv(x_start, 
                                 overlap_slc, 
                                 delta = delta,
                                 y_l = y_l, 
                                 y_h = y_h, 
                                 p = p)$x
    lambda_lh_upper <- overlap_solutions_slc[1]
    lambda_hh_lower <- overlap_solutions_slc[2]
    lambda_hh_upper <- overlap_solutions_slc[3]
    lambda_hl_lower <- overlap_solutions_slc[4]
  } else {
    first_best_solutions_slc <- nleqslv(x_start, 
                                 first_best_slc, 
                                 delta = delta,
                                 y_l = y_l, 
                                 y_h = y_h, 
                                 p = p)$x
    lambda_lh_upper <- first_best_solutions_slc[1]
    lambda_hh_lower <- first_best_solutions_slc[2]
    lambda_hh_upper <- first_best_solutions_slc[3]
    lambda_hl_lower <- first_best_solutions_slc[4]
  }
  
  lambda_seq <- numeric(n)

  if (delta <= delta_1_slc){
    lambda_seq <- xi_seq
  } else if (delta > delta_1_slc & delta <= delta_2_slc){
    for (i in 1:n){
      if (lambda_init >= lambda_hl_lower){
        lambda_seq[i] <- lambda_init * (income_seq[i] == 2) +
          lambda_hh_upper * (income_seq[i] == 1 | income_seq[i] == 4) +
          lambda_lh_upper * (income_seq[i] == 3)
      } else if (lambda_init < lambda_hl_lower & lambda_init >= lambda_hh_upper){
        lambda_seq[i] <- lambda_hl_lower * (income_seq[i] == 2) +
          lambda_hh_upper * (income_seq[i] == 1 | income_seq[i] == 4) +
          lambda_lh_upper * (income_seq[i] == 3)
      } else if (lambda_init < lambda_hh_upper & lambda_init >= lambda_hh_lower){
        lambda_seq[i] <- lambda_hl_lower * (income_seq[i] == 2) +
          lambda_init * (income_seq[i] == 1 | income_seq[i] == 4) +
          lambda_lh_upper * (income_seq[i] == 3)
      } else if (lambda_init < lambda_hh_lower & lambda_init >= lambda_lh_upper){
        lambda_seq[i] <- lambda_hl_lower * (income_seq[i] == 2) +
          lambda_hh_lower * (income_seq[i] == 1 | income_seq[i] == 4) +
          lambda_lh_upper * (income_seq[i] == 3)
      } else {
        lambda_seq[i] <- lambda_hl_lower * (income_seq[i] == 2) +
          lambda_hh_lower * (income_seq[i] == 1 | income_seq[i] == 4) +
          lambda_init* (income_seq[i] == 3)
      }
    }
  } else if (delta > delta_2_slc & delta <= delta_3_slc) {
    for (i in 1:n){
      if (lambda_init >= lambda_hh_upper){
        lambda_seq[i] <- lambda_init * (income_seq[i] == 2) +
          lambda_hh_upper * (income_seq[i] == 1 | income_seq[i] == 4) +
          lambda_lh_upper * (income_seq[i] == 3)
      } else if (lambda_init < lambda_hh_upper & lambda_init >= lambda_hl_lower){
        lambda_seq[i] <- lambda_init * (income_seq[i] == 2) +
          lambda_init * (income_seq[i] == 1 | income_seq[i] == 4) +
          lambda_lh_upper * (income_seq[i] == 3)
      } else if (lambda_init < lambda_hl_lower & lambda_init >= lambda_lh_upper){
        lambda_seq[i] <- lambda_hl_lower * (income_seq[i] == 2) +
          lambda_init * (income_seq[i] == 1 | income_seq[i] == 4) +
          lambda_lh_upper * (income_seq[i] == 3)
      } else if (lambda_init < lambda_lh_upper & lambda_init >= lambda_hh_lower){
        lambda_seq[i] <- lambda_hl_lower * (income_seq[i] == 2) +
          lambda_init * (income_seq[i] == 1 | income_seq[i] == 4) +
          lambda_init * (income_seq[i] == 3)
      } else {
        lambda_seq[i] <- lambda_hl_lower * (income_seq[i] == 2) +
          lambda_hh_lower * (income_seq[i] == 1 | income_seq[i] == 4) +
          lambda_init * (income_seq[i] == 3)
      }
    }
  } else {
    for (i in 1:n){
      if (lambda_init >= lambda_hh_upper){
        lambda_seq[i] <- lambda_init * (income_seq[i] == 2) +
          lambda_hh_upper * (income_seq[i] == 1 | income_seq[i] == 4) +
          lambda_lh_upper * (income_seq[i] == 3)
      } else if (lambda_init < lambda_hh_upper & lambda_init >= lambda_lh_upper){
        lambda_seq[i] <- lambda_init * (income_seq[i] == 2) +
          lambda_init * (income_seq[i] == 1 | income_seq[i] == 4) +
          lambda_lh_upper * (income_seq[i] == 3)
      } else if (lambda_init < lambda_lh_upper & lambda_init >= lambda_hl_lower){
        lambda_seq[i] <- lambda_init * (income_seq[i] == 2) +
          lambda_init * (income_seq[i] == 1 | income_seq[i] == 4) +
          lambda_init * (income_seq[i] == 3)
      } else if (lambda_init < lambda_hl_lower & lambda_init >= lambda_hh_lower){
        lambda_seq[i] <- lambda_hl_lower * (income_seq[i] == 2) +
          lambda_init * (income_seq[i] == 1 | income_seq[i] == 4) +
          lambda_init * (income_seq[i] == 3)
      } else {
        lambda_seq[i] <- lambda_hl_lower * (income_seq[i] == 2) +
          lambda_hh_lower * (income_seq[i] == 1 | income_seq[i] == 4) +
          lambda_init * (income_seq[i] == 3)
      }
    }
  }
  # Sequence of net transfers (1 -> 2), income, and consumptions
  transfer_seq <- numeric(n)
  inc_seq_1 <- numeric(n)
  inc_seq_2 <- numeric(n)
  cons_seq_1 <- numeric(n)
  cons_seq_2 <- numeric(n)
  for (i in 1:n){
    if (income_seq[i] == 1){
      transfer_seq[i] <- (y_h - y_h * lambda_seq[i]) / (1 + lambda_seq[i])
      inc_seq_1[i] <- y_h
      inc_seq_2[i] <- y_h
      cons_seq_1[i] <- y_h - transfer_seq[i]
      cons_seq_2[i] <- y_h + transfer_seq[i]
    } else if (income_seq[i] == 2) {
      transfer_seq[i] <- (y_h - y_l * lambda_seq[i]) / (1 + lambda_seq[i])
      inc_seq_1[i] <- y_h
      inc_seq_2[i] <- y_l
      cons_seq_1[i] <- y_h - transfer_seq[i]
      cons_seq_2[i] <- y_l + transfer_seq[i]
    } else if (income_seq[i] == 3) {
      transfer_seq[i] <- (y_l - y_h * lambda_seq[i]) / (1 + lambda_seq[i])
      inc_seq_1[i] <- y_l
      inc_seq_2[i] <- y_h
      cons_seq_1[i] <- y_l - transfer_seq[i]
      cons_seq_2[i] <- y_h + transfer_seq[i]
    } else {
      transfer_seq[i] <- (y_l - y_l * lambda_seq[i]) / (1 + lambda_seq[i])
      inc_seq_1[i] <- y_l
      inc_seq_2[i] <- y_l
      cons_seq_1[i] <- y_l - transfer_seq[i]
      cons_seq_2[i] <- y_l + transfer_seq[i]
    }
  }
  plot_test <- ggplot() +
                geom_line(aes(delta_vals,
                              log(lambda_lh_lower_vals_slc),
                              color="a")) + 
                geom_line(aes(delta_vals,
                              log(lambda_lh_upper_vals_slc),
                              color="b")) + 
                geom_line(aes(delta_vals,
                              log(lambda_hh_lower_vals_slc),
                              color="c")) + 
                geom_line(aes(delta_vals,
                              log(lambda_hh_upper_vals_slc),
                              color="d")) +
                geom_line(aes(delta_vals,
                              log(lambda_hl_lower_vals_slc),
                              color="e")) +
                geom_line(aes(delta_vals,
                              log(lambda_hl_upper_vals_slc),
                              color="f")) +
                coord_cartesian(xlim = c(0.8, 1.0),
                                ylim = c(log(y_l / y_h), log(y_h / y_l))) +
                geom_ribbon(aes(x = delta_vals, 
                                ymin = log(lambda_lh_lower_vals_slc), 
                                ymax = log(lambda_lh_upper_vals_slc)),
                                fill = "blue", alpha = "0.2") +
                geom_ribbon(aes(x = delta_vals, 
                                ymin = log(lambda_hh_lower_vals_slc), 
                                ymax = log(lambda_hh_upper_vals_slc)),
                                fill = "red", alpha = "0.2") +
                geom_ribbon(aes(x = delta_vals, 
                                ymin = log(lambda_hl_lower_vals_slc), 
                                ymax = log(lambda_hl_upper_vals_slc)),
                                fill = "green", alpha = "0.2") + 
                geom_vline(xintercept = delta, 
                           color = "black", 
                           size = 0.8, 
                           alpha = 0.6) +
                scale_color_manual(name = "End-points", 
                                   values = c("blue", "purple", "brown", "red", "yellow", "green"), 
                                   labels = c("lambda_lh_lower", "lambda_lh_upper", 
                                              "lambda_hh_lower", "lambda_hh_upper", 
                                              "lambda_hl_lower", "lambda_hl_upper")) +
                xlab("Discount factor (delta)") +
                ylab("log of the ratio of marginal utilities (lambda)") +
                geom_point(aes(rep(delta, (n + 1)), c(log(lambda_init), log(lambda_seq)))) +
                geom_label_repel(
                  aes(rep(delta, (n + 1)), c(log(lambda_init), log(lambda_seq)), label = seq(0,n)),
                  box.padding = 0.35, point.padding = 0.5)
  
  result_table <- data.frame(cbind(seq(0,n), 
                                   round(c(log(lambda_init), log(lambda_seq)), 3), 
                                   c(NaN, income_realization_label_seq),
                                   round(c(NaN, transfer_seq), 3),
                                   round(c(NaN, cons_seq_1), 3),
                                   round(c(NaN, cons_seq_2), 3)
                                    ))
  colnames(result_table) <- c("Period", 
                              "ln(lambda)", 
                              "Income shocks", 
                              "Net transfer (1 -> 2)", 
                              "Consumption (1)", 
                              "Consumption (2)")

  summary_table <- matrix(c(mean(cons_seq_1), 
                            std(cons_seq_1), 
                            mean(cons_seq_2), 
                            std(cons_seq_2), 
                            mean(inc_seq_1), 
                            std(inc_seq_1), 
                            mean(inc_seq_2), 
                            std(inc_seq_2)), byrow = TRUE, ncol = 2)
  summary_table <- round(summary_table, 3)
  colnames(summary_table) <- c("Mean", "Std")
  rownames(summary_table) <- c("Consumption (1)", "Consumption (2)", "Income (1)", "Income (2)")
  return(list(plot_test, result_table, summary_table))
}
```

#### Plot the figure and compare!

I plot the figure under the static limited commitment first, and for comparison I plot the figure under the dynamic limited commitment too. I try a discount factor so that in both models we are in the non-overlapping region. This is because in this region the difference is the most salient.

``` r
n <- 10
lambda_init <- lambda_0
delta <- 0.925
set.seed(1)
result_slc <- lambda_change_plot_slc(lambda_init, delta, n)
plot_figure_slc <- result_slc[1]
result_table_slc <- result_slc[2]
summary_table_slc <- result_slc[3]

set.seed(1)
result <- lambda_change_plot(lambda_init, delta, n)
plot_figure <- result[1]
result_table <- result[2]
summary_table <- result[3]
```

First, I show the result in the static model.

``` r
kable(result_table_slc, caption = "Static limited commitment model")
```

<table class="kable_wrapper">
<caption>
Static limited commitment model
</caption>
<tbody>
<tr>
<td>
| Period | ln(lambda) | Income shocks | Net transfer (1 -&gt; 2) | Consumption (1) | Consumption (2) |
|:-------|:-----------|:--------------|:-------------------------|:----------------|:----------------|
| 0      | 0          | NaN           | NaN                      | NaN             | NaN             |
| 1      | 0          | High, High    | 0                        | 1.333           | 1.333           |
| 2      | 0          | High, High    | 0                        | 1.333           | 1.333           |
| 3      | 0          | High, High    | 0                        | 1.333           | 1.333           |
| 4      | 0.592      | High, Low     | 0.046                    | 1.288           | 0.712           |
| 5      | 0          | High, High    | 0                        | 1.333           | 1.333           |
| 6      | -0.592     | Low, High     | -0.046                   | 0.712           | 1.288           |
| 7      | 0.592      | High, Low     | 0.046                    | 1.288           | 0.712           |
| 8      | 0          | High, High    | 0                        | 1.333           | 1.333           |
| 9      | 0          | High, High    | 0                        | 1.333           | 1.333           |
| 10     | 0          | High, High    | 0                        | 1.333           | 1.333           |

</td>
</tr>
</tbody>
</table>
``` r
plot_figure_slc
```

    ## [[1]]

![](LTW_code_files/figure-markdown_github/unnamed-chunk-27-1.png)

Next, I show the result in the dynamic model. Note that the income shock sequence is common in two models.

``` r
kable(result_table, caption = "Dynamic limited commitment model")
```

<table class="kable_wrapper">
<caption>
Dynamic limited commitment model
</caption>
<tbody>
<tr>
<td>
| Period | ln(lambda) | Income shocks | Net transfer (1 -&gt; 2) | Consumption (1) | Consumption (2) |
|:-------|:-----------|:--------------|:-------------------------|:----------------|:----------------|
| 0      | 0          | NaN           | NaN                      | NaN             | NaN             |
| 1      | 0          | High, High    | 0                        | 1.333           | 1.333           |
| 2      | 0          | High, High    | 0                        | 1.333           | 1.333           |
| 3      | 0          | High, High    | 0                        | 1.333           | 1.333           |
| 4      | 0.164      | High, Low     | 0.252                    | 1.082           | 0.918           |
| 5      | 0.057      | High, High    | -0.038                   | 1.372           | 1.295           |
| 6      | -0.164     | Low, High     | -0.252                   | 0.918           | 1.082           |
| 7      | 0.164      | High, Low     | 0.252                    | 1.082           | 0.918           |
| 8      | 0.057      | High, High    | -0.038                   | 1.372           | 1.295           |
| 9      | 0.057      | High, High    | -0.038                   | 1.372           | 1.295           |
| 10     | 0.057      | High, High    | -0.038                   | 1.372           | 1.295           |

</td>
</tr>
</tbody>
</table>
``` r
plot_figure
```

    ## [[1]]

![](LTW_code_files/figure-markdown_github/unnamed-chunk-28-1.png)

There are several things that should be noted. First, in the static limited commitment model, the payment is not path dependence: only current incomes matter for transfers. Second, in the static model, when possible, the marginal utilities are equated. While this appears to be a nice thing, this comes at the cost of larger fluctuation of *λ*. This results in more fluctuating consumptions as shown in the table.

``` r
summary_table_both <- cbind(summary_table[[1]], summary_table_slc[[1]])
colnames(summary_table_both) <- c("Mean (dynamic)", 
                                  "Std (dynamic)", 
                                  "Mean (static)", 
                                  "Std (static)")
kable(summary_table_both, caption = "Summary statistics")
```

|                 |  Mean (dynamic)|  Std (dynamic)|  Mean (static)|  Std (static)|
|-----------------|---------------:|--------------:|--------------:|-------------:|
| Consumption (1) |           1.257|          0.165|          1.262|         0.194|
| Consumption (2) |           1.210|          0.170|          1.205|         0.260|
| Income (1)      |           1.267|          0.211|          1.267|         0.211|
| Income (2)      |           1.200|          0.281|          1.200|         0.281|
