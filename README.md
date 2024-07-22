Code for Experiment in Publication “Determining the Dominant Factors Controlling Mineralization in Three-Dimensional Fracture Networks”
========

All code used to perform the experiment in the cited publication.  For questions, issues, or clarifications please reach out to Murph: <murph@lanl.gov>.

This repository is organized according to the different directories needed to perform the studies in the paper.  These directories are:

* generate_dfns: code to perform the DFN simulations to create the data used in the experiment;

* fixed_dfn: code to perform the experiments on a single DFN (see paper);

* variable_dfn: code to perform the experiments the underlying DFN is allowed to vary (see paper), on a reparameterized time scale;

* variable_dfn_regulartime: code to perform the experiments the underlying DFN is allowed to vary (see paper), on regular time;

* paper_stuff: the file run_all_plots.R write all figures to this directory.

## Note on Reproducibility

We suggest cloning this repository into a folder entitled GitHub in your home directory.  Filepaths within this directory structure should then work out of the box.

The figures in the paper can be recreated by running the file run_all_plots.R.  All the data from dfnWorks is pre-compiled into .pkl files in the directories  
* fixed_dfn/data

* variable_dfn/data

* variable_dfn_regulartime/data

This is due to the very large size of a typical dfnWorks output; we compiled this data using scripts in this repository, removing extraneous output for space considerations.  The raw output can be produced using the dfnWorks software suite, using the files in generate_dfns.  The corresponding seeds a specific inputs can be found in these directories.

## Software Required
To create underground particle transport simulation data, one will need access to the dfnWorks simulation suite, available for download [here](https://dfnworks.lanl.gov/). The dfnWorks data files are included in this repository.

## Citation
Jeffrey D. Hyman, Alexander C. Murph, Lawrence Boampong, Alexis Navarre-Sitchler, James W. Carey, Phil Stauffer, Hari S. Viswanathan. (2024). Determining the Dominant Factors Controlling Mineralization in Three-Dimensional Fracture Networks.  _In Review._ 

## Release

This software has been approved for open source release and has been assigned **O4756** 

## Copyright

© 2024. Triad National Security, LLC. All rights reserved.
This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S. Department of Energy/National Nuclear Security Administration. All rights in the program are reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear Security Administration. The Government is granted for itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare derivative works, distribute copies to the public, perform publicly and display publicly, and to permit others to do so.

## License

This code repository is distributed under a MIT License:

Copyright 2024
Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

