# CANDY (Chemistry ANd DYnamics)

CANDY is a protoplanetary disk model for determining the time-dependent chemistry and a 1D vertical slice of a dynamic disk. CANDY includes chemistry, vertical diffusion, pebble growth, and ice sequestration. Details of the model are given in Van Clepper et al. 2022.

CANDY contains two submodules, a modified astrochem (Maret & Bergin 2015) directory and the chemdiff directory, where the python wrapper is contained.

## Installation and Requirements

Astrochem dependencies are given in the astrochem wiki, while chemdiff dependencies include openmpi and mpi4py. Recommended installation is using conda to create a virtual environment with all dependencies. A new environment can be created with all needed dependencies with

	conda create --name candy -c conda-forge sundials python cython numpy matplotlib h5py autoconf automake libtool mpi4py openmpi

Active the enviornment with

	conda activate candy

or

	source activate candy

To compile Astrochem run the following commands *from the Astrochem subdirectory*

	./bootstrap
	./configure CPPFLAGS="-I$CONDA_PREFIX/include" LDFLAGS="-Wl,-rpath,$CONDA_PREFIX/lib -L/$CONDA_PREFIX/lib"  --prefix=$CONDA_PREFIX
	make
	make install

## Astrochem

The astrochem folder is copied from the (astrochem)[https://github.com/smaret/astrochem] github page, complete with documentation and installation instructions. Changes have been made within the `src` folder, and Astrochem can be installed as usual. If you already have astrochem installed on your machine, you can simply copy the `src` folder from here into your astrochem directory (replacing the default astrochem/src folder), then reinstall as usual.

Chemical network files (`.chm`) should follow the same format as astrochem but with added chemical reactions for:

|description | reacno|
|------------:|:-------|
| HD formation on grains | 100 |
| Shielded dissociation of H2 | 14 |
| Shielded dissociation of HD | 114 |
| Shielded photo-ionization of C | 15 |
| Shielded dissociation of C16O | 16 |
| Shielded dissociation of C17O | 17 |
| Shielded dissociation of C18O | 18 |
| Hydrogenation | 25 |
| Cosmic-ray desorption of CO | 42 |
| Cosmic-ray reactions involving He | 43 |
| secondary xray ionization of ... H | 60 |
| ... H2 | 61 |
| ... other molecules | 62 |
| creation of excited H2 | 90 |
| de-excitation of H2 | 91 | 
| reactions with excited H2 | 92 |

Some reaction rate calculations are also adjusted to be more similar to the DALI calculations (Burderer et al 2012). A full description of the handling of the chemical calculations is in the appendix of Van Clepper et al. 2022.

## chemdiff

The python wrapper for for running astroCHEM with DIFFusion. This is entirely written in python, and successively calls `astrochem` in parallel at different vertical locations ranging from disk midplane to 5 scale heights. Input parameters can be read in from the `cdinput.in` file. Full documentation is located in the chemdiff directory.

To run, make sure both chemdiff and astrochem are in your path, then place the `run_parallel_growth.py` and `cdinput.in` files into your desired directory and run:

	python run_parallel_growth.py

or

	mpiexec python run_parallel_growth.py

An example sbatch submission script is also included.

The outputs will be placed in a subdirectory `./r00/` from where `run_parallel_growth.py` was run. This will take up a lot of space. the `get_all_abundances.py` module in the `post_processing` subdirectory will consolidate the abundances, species, times, heights, gas densities, and visual extinctions into an npz file. Just call

	python get_all_abundances.py {path to directory}

You can also do

	python get_all_abundances.py -h

to see more options. Once this has been called the `./r00/` directory can safely be removed.