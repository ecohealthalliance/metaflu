# Design Notes

Approach:

-  Recreate model in https://dx.doi.org/10.1371/journal.pone.0080091
    -  Original is in C, available on Aegypti server but difficult to decipher
-  Simulation via 2 structures: location by disease class matrix and connectivity matrix
-  Implement simulation via Rcpp, parallelism within R
-  Allow for simulations to stochastically generate network via parameters or take
   existing network as input
-  Add term for external regular source of disease
-  Use a Monte Carlo sensitivity approach to deal with uncertainty in parameters, esp ones
   we do not have. 
-  Ensure that simulations are reproducible with seed, parallel or not, let
   user set seed, memoize seeded simulation calls

Questions:

-  What is our measure of interest?
   -   Total infected birds/year?
   -   Probability of outbreak, overall or by farm?
-  How should patch size distribution be determined? Hopefully from FAO data.
-  Implementing the control program - what's realistic compared to actual
   country measures?
-  What are appropriate parameters for the poultry we have (ducks, chickens, etc.)?  
-  Upstream and downstream
   -   What is probability of initial infection?  Need info on wild birds.  See
       papers by N Gaidet (also looks for anything by Hosseini, Karesh)
   -   What is probability of human infections (primary cases) and human outbreak?  


Outputs:

-  Interactive/App: (Probably) simulation is too compuationally expensive to run
   for easily distributed app.  Save lots of simulations' summary data that can
   be explored?
   
Optimization:

Notes:
-  Implement Cao et al 2007 without the implicit tau-leaping, just the
  critical variables and explicity poisson leaping.  Interim idea may
  be to approximate implicit via Runge-Kutta type correction on the
  explicit approach.
-  Switch the arrangement so that all matrix-specific calculations happen inside
  the rate and update functions.  The algorithm should deal with a vector only
-  As for parallelism, consider doing nsims inside the C++ code using OpenMP
  (but leave multiple parameterizations outside - what about initiating? May
  have to call the generation function in R? Initiate before simulating?)
-  Look at the fastmvn package approach to parallel, reproducible RNG
-  Consider an arma sparse matrix for defining actions and networks
-  Reduction function can also occur in OpenMP, but separately
-  Reduction functions via Xptr?
