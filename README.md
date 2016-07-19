# Variational Monte Carlo for 1D Heisenberg model
using Haldane-Shastry variational wavefunction 

Compile as: `make vmc`

To run: `./vmc.x [Char(20) UUID of the run] [Int seed]`

Example:
`./vmc.x 'testRun' 3127387 < model.in > sim.out`

Starts simulation for 1D Heisenberg system of **L** sites as given in `model.in`. The simulation performs **nOPT** iterations of optimization cycle, during which
the new updated variational parameters **p** are computed. The magnitude of update is proportional to value of **opt_eps**. Each single iteration consists of **nEQ** equilibration sweeps and **nPROD** production sweeps for   averaging. The data from optimization are stored in `testRun-opt.dat` file.
    
If the value of **nOPT**= **1**, the samples of energy and spin-spin correlation functions are saved to files `testRun-vals.dat` and `testRun-corr.dat` respectively for subsequent analysis. The average values of energy and spin-spin correlation functions are given in `sim.out` as well.

### Data analysis:

##### Sample means

`bining.py` computes average and estimates the error of sampled observable by decorrelating the data using binning analysis

To run: `python3 binning.py [String filename] [Int column] [Int #samples] [Int MaxBinSize]`

Example: `python3 binning.py 'testRun-vals.dat' 1 1000 50 > averages.dat`

Performs binning analysis of samples of energy stored in column No. 1 (starting from 0th column)

##### References
\[1\] A. W. Sandvik, arXiv:1101.3281 

\[2\] S. Sorella, arXiv:cond-mat/0502553