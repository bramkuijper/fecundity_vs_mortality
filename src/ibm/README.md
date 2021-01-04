# Parental effects and fecundity versus mortality selection

Bram Kuijper, a dot l dot w dot kuijper at exeter.ac.uk

## Individual-based simulations (Gillespie model)

Overview of the various files:

- `p1p2_normal_n.cpp`: the simulation itself
- `generate_sims.py`: `Python3` script to create batch files with which to run the simulations
- `plot_sims.py`: `Python3` script to plot a single simulation run
- `auxiliary.h`: helper functions to generate unique file names 
- `makefile`: required to compile the code using GNU `make`. 

### Requirements

To run the `C++` code:

- a `C++` compiler (e.g., g++) that supports `C++11`
- the GNU scientific library (`GSL`) for random number generation
- GNU `make`
- A terminal like `bash`

Simulations ran on a UNIX-based system with gcc 
and the GNU scientific library (GSL) installed. To run everything
smoothly on MS-Windows, take a look at MinGW http://www.mingw.org
and build GNU GSL within that system (see https://stackoverflow.com/questions/30015656/using-gnu-scientific-library-gsl-under-windows-x64-with-mingw). 

## How to compile

Open a terminal and type

    cd <top-level-directory>/src/ibm
    make 

The directory should now contain the executable file `xp1p2n`
which can be run as 

    ./xp1p2n 0.08164965809277261 0.1224744871391589 1.0 1.0 1.0 1.0 1.0 1.0 2.0 2.0 0  0.02 0.01 0.25 0 0.1 0 1.0 

i.e., with quite a bunch of command line parameters (see the
`init_arguments()` function in the `p1p2_normal_n.cpp` file for
more info about those).

## How to run multiple instances of the simulation
To generate combinations of parameters provided to the command
line I have written a `Python3` script
`generate_sims.py`, for which you need `Python3` installed with
the `numpy` package. On the command line:

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
