##Olfactory White Algorithm##

Old set of scripts I wrote during my PhD to build an olfactory white recipe based on the vibrational spectra
of a number of volatile molecules. 

Starting from a training dataset of molecular spectra, and assuming their sum is an olfactory white,  
the algorithm iteratively adds or subtracts molecules (with a metropolis acceptance/regection algorithm)
and tries to minimise the distance between the newly built total spectrum and the spectrum 
obtained from the training set according to the chosen metric. 

A set of auxiliary scripts to test the validity of the results are also here.

List of files:

- Owl2.py the main file to built the olfactory white;
- Owl_core.py contains the functions used in the other scripts;
- Owl2_multirun.py to repeat the search a number of times and get some statistics;
- OwlCheck2.py to re-evaluate the error on an ouput recipe;
- OwlErrors2.py to evaluate the error between different given white recipies;
- OwlRandom2.py to evaluate the error on a random set of molecules;
