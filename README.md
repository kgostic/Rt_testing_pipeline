# Tests for R0 and R(t) estimation

* `simulation.R` contains functions to simulate data using a deterministic or stochatic SIR or SEIR model. The simulations can input a vector of time-varying R0 values.
* `experiments/` contains notebooks and their dependencies, which call simulation.R to test the performance of specific $R_t$ estimators.

Best practices:

1. Create a new directory for your experiment.
2. Run your experiment-specific script in that directory.
3. Commit your experiment-specific script to git; do not commit outputs.
