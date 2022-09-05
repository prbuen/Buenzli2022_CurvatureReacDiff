# Buenzli2022_CurvatureReacDiff

![Square pore t=7](https://github.com/prbuen/Buenzli2022_CurvatureReacDiff/blob/main/frame-0070.png)

### Author
Pascal R Buenzli, 2022

### Description
Buenzli PR and Simpson MJ (2022) Curvature dependences of wave propagation in reaction-diffusion models;
Preprint available at https://arxiv.org/abs/2112.00928

### Requirements
- scid library 0.3.x, https://code.dlang.org/packages/scid
- D compiler: ldc2 v1.24.0 (DMD v2.094.1, LLVM 11.0.1) or equivalent dmd, see https://wiki.dlang.org/LDC
- Visualisation script written for Gnuplot 5.4.1
- Tested on MacOSX and Linux (Kubuntu 21.10)

### Compilation
Edit src/model_inputs.d to select the model to run, then, from the root directory:

    ldc2 src/*.d -od=build
    
or (optimised):

    ldc2 src/*.d -od=build -O3 -release -mcpu=native -ffast-math

### Execution

    ./main
    
=> data created in `datadir/*.dat` (filenames are suffixed with time frame number)

### Visualisation
Edit plot.py (e.g. frame=70) and feed to gnuplot

    ./plot.py | gnuplot
    pdflatex frame-0070
    
=> frame-0070.pdf

This produces the above figure.
