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