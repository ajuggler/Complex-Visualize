# Complex-Visualize
 
A package for visualizing arbitrary complex functions.  Implements the strategy of domain coloring on the Riemann sphere.

Authors:  Ma. de los Angeles Sandoval-Romero and Antonio Hernandez-Garduno.

## Introduction

*Domain coloring* is a technique for constructing a tractable visual object of the graph of a complex function.  The package *complexVisualize.wl* improves on existing domain coloring techniques by rendering a global picture fo the Riemann sphere (the compactification of the complex plane).  Additionally, the package allows dynamic visualization of families of Moebius transformations.

To learn more about this package and view some nice examples, refer to [this article](https://www.mathematica-journal.com/2015/11/30/domain-coloring-on-the-riemann-sphere/).

## Quick reference

### Requirements
Mathematica version 8.0 or newer

### Installation instructions
- Clone this repository
- Use the menu "File > Install..." to install the package `complexVisualize.wl` found in the cloned directory.

### Usage
- In a *Mathematica* session, load the package with:

      <<complexVisualize.wl

- To plot the function *f(z)*:

      complexVisualize[f[z], z]

  For example, the identity function (which is useful to try as a first example) is visualized with

      complexVisualize[z, z]
    
- To get further usage instructions, type

      ?complexVisualize
