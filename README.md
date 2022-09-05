# Buenzli2022_CurvatureReacDiff

### Author
Pascal R Buenzli, 2022

### Description
Buenzli PR and Simpson MJ (2022) Curvature dependences of wave propagation in reaction-diffusion models;
Preprint available at https://arxiv.org/abs/2112.00928

### Requirements
- scid library 0.3.x, https://code.dlang.org/packages/scid
- D compiler: ldc2 or dmd
- Visualisation done with Gnuplot 5.4.1 (others possible)
- Tested on MacOSX and Linux (Kubuntu 21.10)

### Compilation
Edit src/model_inputs.d to select the model to run, then, from the root directory:

    ldc2 src/*.d src/pb/*.d -od=build
    
or (optimised):

    ldc2 src/*.d src/pb/*.d -od=build -O3 -release -mcpu=native -ffast-math

### Execution

	  ./main
    
=> data created in `datadir/*.dat` (filenames are suffixed with time frame number)

### Visualisation
Edit plot.py (e.g. frame=70) and feed to gnuplot

    ./plot.py | gnuplot
    pdflatex frame-0070 (twice)
    
=> frame-0070.pdf

This produces the following figure:
![Square pore t=7](https://github.com/prbuen/Buenzli2022_CurvatureReacDiff/blob/master/frame-0070.png)
