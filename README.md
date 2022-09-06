The codes contained in this repository can evaluate the momentum-dependent
superconducting (SC) gap for the Kagome model as presented in the arxiv preprint [arXiv 2208.13521](https://arxiv.org/abs/2208.13521) and as described in this notes:
https://drive.google.com/file/d/1Yo49f6vaBZoAiL0a5J1kWhHoxUHX4__i/view?usp=sharing

Individual functions are written in 3 different programming languages:
Python (version 2.7), Fortran90, MATLAB (version 2019a)

The codes should be run in the following chronological order:

1. chio_parent.py: It is a parent glue script that takes in inputs
   (qxsus,kxsus,orbitals,Ef,MatsuT,nproc) and returns the output (baresus.mat)
   for evaluating the bare susceptibility. For each q-point evaluation,
   it is parallelized over 'nproc' number of processors using 'subprocess
   Popen' module.

   For a 80 X 80 q-grid and one orbital model parallelized over 45 processors,
   time required for the completion of this part of the code was ~4 hours.

2. pairing_parent.py: Takes in inputs (ncont,U,J,Vnn) and yields output
   (gapstructure.mat) for the final SC gap structure on the Fermi surface.

   For each set of U-J-V values and for ncont=150 Fermi surface points, time
   required for the one-orbital model evaluation was ~3 hours.

Note that if every other parameter is kept unchanged except for the number of
orbitals 'orbitals', time and resource requirement scales as (orbitals)^4.
Going through each of the above scripts will give further idea of what other
functions and individual scripts are doing.

Once all results are evaluated, once can use 'plot_finalgap.m' MATLAB function
to plot the SC gap structure on the Fermi surface.
The folders 'bareresults', 'functions', 'images', and 'input_files' have
additional scripts/images/input/output files.
