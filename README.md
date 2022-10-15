# Brownian Surface Mediated Diffusion Simulations

## Description

This Rust based python module is built with the help of Rust's **maturin** binder. It provides a few functions callable from any python interpreter to perform numerical simulations of surface mediated diffusion of a Brownian Particle in 3 dimensionnal spherical geometries. 

The module also includes the possibility to add spherical cheese-holes in the interior of the domain, onto which the Particle can absorb and perform surface diffusion until desorption. Observables among which splitting probabilities, first passage times and steady-state state distribution are computed using a mixture of time-driven and event-driven algorithms.  

## Installation and documentation

Instalation is straighforward. Make sure you have Maturin installed, or run 
```bash
pip install maturin
```
to install. Clone the repository, cd into it and run :
```bash
maturin develop
```
If all dependencies are satisfied, it should work! 😇 To find documentation, run:
```python
import smd-3d as smd
help(smd)
```

## disclaimer:

The code is under constant developpement, hence the astonishing volume of commented-out content. Please bear with me on this one 😬😬
