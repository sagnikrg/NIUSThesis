# NIUSThesis
This repository contains all the relevant files, codes and graphs related to the National Initiative on Undergraduate Science (NIUS) Fellowship project pursued at Homi Bhabha Centre for Science Education (TIFR) between December 2016 and December 2019.

We use a mixture of Mathematica and C. Using Mathematica we simplify all the analytical computation to express the Path Integral in terms of known functions of the classical path, then using C we perform the quadratures.

## Files

1. **Path Action.nb:** This Mathematica file contains the general algorithm, where we numerically solve the Classical Path as a BVP on a grid and compute the Propagator using the same.
2. **Module 2:** Sample notebook on the computation of wave function, using pre-store classical path.
3. **Module 3:** Plot File. Sample notebook to generate plots of P(x,t) and Q(x,t) using pregenerated wave functions from Module 2.
4. **ClassicalPert.nb:** Notebook generating plots of the analytical perturbative treatment for the Anharmonic Oscilator. 
