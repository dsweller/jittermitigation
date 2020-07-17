# jittermitigation
Bayesian post-processing methods for sampling jitter mitigation

This software package contains code described in two published papers:

1. An Expectation-Maximization (EM) algorithm ([emalg.m](emalg.m)) described in: Weller, Daniel S. and Goyal, Vivek K. "On the estimation of nonrandom signal coefficients from jittered samples." *IEEE Trans. Signal Process.*, vol. 59, no. 2, pp. 587-597, Feb. 2011. DOI: [10.1109/TSP.2010.2090347](https://doi.org/10.1109/TSP.2010.2090347)

2. A Markov chain Monte Carlo (MCMC) algorithm ([gibbsslicesampler.m](gibbsslicesampler.m)) described in: Weller, Daniel S. and Goyal, Vivek K. "Bayesian post-processing methods for jitter mitigation in sampling." *IEEE Trans. Signal Process.*, vol. 59, no. 5, pp. 2112-2123, May 2011. DOI: [10.1109/TSP.2011.2108289](https://doi.org/10.1109/TSP.2011.2108289)

In addition, the script [test_performance.m](test_performance.m) is provided as an example of how to use this software, as well as reproduce performance comparisons such as those found in the above papers. This script uses in addition the gammar() function from version 4.2 of the StatBox toolbox available for download at http://www.statsci.org/matlab/statbox.html or on the MATLAB Central Link Exchange. StatBox is copyright (c) 1996-2003, Gordon Smyth. *If the toolbox does not already reside on your MATLAB path, you need to adjust the path in this script to point to this toolbox.*

The file [gen_laguerre_rule.m](gen_laguerre_rule.m) provided is based on software by John Burkardt, distributed under the GNU LGPL license.

This program is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public [License](LICENSE) as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

These are MATLAB (The Mathworks, Natick, MA) scripts. They are likely compatible (with the possible exception of certain plotting commands) with [GNU Octave](https://www.gnu.org/software/octave/).
