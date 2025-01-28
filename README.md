# adaptive_cooperative_control_MAS_hyperbolic

This repository implements methods for the simulation and the design of adaptive control algorithms presented in the papers [1,2] for the output regulation of multi-agent systems, in which every agent has hyperbolic distributed-parameter dynamics.

# Dependencies

The programms were tested and developed using Matlab R2024b.
The following Matlab toolboxes and external libraries are used.

### Matlab toolboxes

Frequently, the following Matlab toobloxes are used: control_toolbox, optimization_toolbox, symbolic_toolbox.

### External toolboxes (not included)

* "coni" available at [https://gitlab.com/control-system-tools/coni] (see also [3])
* "hyperbolic" available at [https://gitlab.com/control-system-tools/hyperbolic] (see also [4])
* "chebfun" available at [https://github.com/chebfun/chebfun] (see also [5, 6])

# Usage

Include all folders, including the external toolboxes, to your Matlab path.
For an introduction to the frequently employed and very useful quantity.Discrete framework provided by the coni toolbox, see the Example +quantity\example.m of in the coni folder.
Note that the calculation of the backstepping kernel stores all computed kernels in a build-folder in the working directory. 
Hence, it is advised to always work in the same root working folder.
The skrpits that were used to generate the examples that appear in my papers and presentations can be found in the directory +publication_examples, other skripts are in +examples.

# References

[1] Enderes, T.; J. Deutscher: Cooperative robust output regulation for networks of hyperbolic systems
with unknown signal models. accepted for Int. J. Robust Nonlin. Control (2024).
[2] Enderes, T.; J. Gabriel; J. Deutscher: Cooperative output regulation for networks of hyperbolic
systems using adaptive cooperative observers. Automatica 162, Art. no. 111506 (2024).
[3] F. Fischer, J. Gabriel, S. Kerschbaum (2022). coni - a Matlab toolbox facilitating the solution of control problems (1.2). Zenodo.
[4] Gabriel, J.: hyperbolic – A Matlab Toolbox for the Robust Cooperative Output Regulation of Hyperbolic PIDE–ODE Systems. Zenodo. Version v1.0. 2023.
[5] L. N. Trefethen. Spectral Methods in MATLAB. Philadelphia: SIAM, 2000
[6] T. A. Driscoll, N. Hale, and L. N. Trefethen, Chebfun Guide, Pafnuty Publications, Oxford, 2014.
