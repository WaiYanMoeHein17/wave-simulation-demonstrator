# wave-simulation-demonstrator
L3 Project

Working on at the moment: 

1. Adding CAD support to the wave simulation demonstrator to act as a layer between the depth camera taking in the topgography/bathymetry data and the wave simulation engine. 

GOAL: Ideally, to enable users to add things like houses and flood barries to the simulation without having to model them with sand in the tank first. 

2. Extending ADER-DG method to VisualPDE which is in GPU either through SYCL or CUDA. 

GOAL: To enable GPU support for the wave simulation demonstrator

3. Optimising the current wave simulation demonstrator to run faster and more efficiently on CPU. 

GOAL: Quicker simulations

Methods being tried: 
- Kernel Fusion 
- SYCL Offloading

4. Plan to add additional features: 
    - Gravity Acoustics Coupling
    - Non-newtonian Fluids
    - Diffusion for sharp edges
