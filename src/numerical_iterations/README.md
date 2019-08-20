# Parental effects and fecundity versus viability selection

Bram Kuijper, a dot l dot w dot kuijper at exeter.ac.uk

## Numerical iterations in C++

Overview of the various files (all a bit clunky but try to upscale running mathematica code over 1000s of nodes):

- `numsolve.cpp`: the skeleton `C++` file in which mathematica equations will be added using the `Python3` script `generate_cpp.py` to produce `numsolve2.cpp` that is the complete and compilable code
- `numsolve2.cpp`: The main `C++` file that contains the numerical simulations: compilable code generated from`numsolve.cpp` by `generate_cpp.py`. This file does all the numerical work (once compiled that is).
- `auxiliary.h`: `C++` header file with helper functions to generate unique output file names (so that simultaneous runs do not overwrite each other)
- `generate_cpp.py`: `Python3` script to combine the exported `mathematica` equations (present in various `*.txt` files) to `numsolve.cpp` in order to create `numsolve2.cpp`
- `generate_iters.py`: `Python3` script to make a batch file, in which  line contains a command to run a single instance of a numerical iteration
- `plot_sims.py`: `Python3` script to plot a single iteration over time
- `summarize_sims.py`: `Python3` script to summarize output and parameters from multiple files
- `testruns.sh`: example batch file to generate one of the figures in the manuscript
- `makefile` : compile the `C++` code using GNU `make`
- `variables.txt` : (exported from `mathematica` notebook, see there). All the variables that reflect the ecological or evolutionary dynamics
- `params.txt` : (exported from `mathematica` notebook, see there). List with parameters of the model
- `relatedness.txt` : (exported from `mathematica` notebook, see there). Equations of the relatedness coefficients
- `repvals.txt` : (exported from `mathematica` notebook, see there). Equations of the reproductive values
- `selgrad.txt` : (exported from `mathematica` notebook, see there). Equations of the selection gradients
- `patchrecurs.txt` : (exported from `mathematica` notebook, see there). Equations of the patch frequencies 

### Requirements

To run the `C++` code:

- a `C++` compiler (e.g., g++) that supports `C++11`
- the GNU scientific library (`GSL`) 
- GNU `make`
- A terminal like `bash`

Simulations ran on a UNIX-based system with gcc 
and the GNU scientific library (GSL) installed. To run everything
smoothly on MS-Windows, take a look at MinGW http://www.mingw.org
and build GNU GSL within that system (see https://stackoverflow.com/questions/30015656/using-gnu-scientific-library-gsl-under-windows-x64-with-mingw). 

## How to compile

Open a terminal and type

    cd <top-level-directory>/src/numerical_iterations
    make 

The directory should now contain the executable file `numsolve`
which can be run as 

    ./numsolve 1 1 1 1.0 1.0 2.0 2.0 1 1 0.0 0.007088812050083358 1.4106735979665885 0.0 0.0 0.05 0.05 0.05 0.0001 0.16666666666666666 0.16666666666666666 0.16666666666666666 0.16666666666666666 0.16666666666666666 0.16666666666666666 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 0.5 0.5 0.5 0.5 0.5 0.5 0.1 0.05 0.05 
i.e., with a whole sleuth of command line parameters (see the
`main()` function in `numsolve2.cpp` file for
more info about those). I wouldn't type all those values 
out yourself, rather use the `Python3` script `generate_iters.py`.

## How to run multiple instances of the simulation
To generate combinations of parameters provided to the command
line I have written a `Python3` script
`generate_sims.py`, for which you need `Python3` installed with
the `numpy` and `pandas` packages. On the command line:

    ./generate_sims.py > runs.sh

The file `runs.sh` is a batch file, each line of which contains 
a command to run the simulation with a combination of parameters
To run all the parameter combinations in the bash shell file 
`runs.sh`, execute

    bash runs.sh

If too much hassle to install python3 with packages to get a
batch file create it yourself using `R` or whatever language
you prefer.

## How to plot output
Use the `Python3` script `plot_sims.py` (requires the
`Python3` packages `pandas` and `matplotlib`). On the command line

    ./plot_sims.py sim_px_20_8_2019_103126_300641676

The plot is then written to the file `graph_sim_px_20_8_2019_103126_300641676.pdf`

If you have multiple files there is also a script `summarize_sims.py` to summarize
output from the last line of multiple files into a single aggregrate output file
