<p align="center">
  <img src="https://bjwmills.com/wp-content/uploads/2023/09/SCION_banner.png" alt="SCION" style="width:700px;"/>
</p>

SCION (Spatial Continuous IntegratiON) is a global climate-biogeochemical model that runs over geological timescales. It is a mixed dimensional system with a steady state 3D interpolated climate, 2D land surface and nondimensional ocean. It runs forwards in time and computes the Earth’s major elemental cycles and surface climate. It also predicts the values of a suite of geochemical tracers to aid in hypothesis testing. 

<p align="center">
  <img src="http://bjwmills.com/wp-content/uploads/2023/09/SCION_outline.png" alt="SCION" style="width:700px;"/>
</p>

v1.3 - October 2025

### Requires MATLAB. 
Current MATLAB package requirements (local machine only – HPC use does not return plots)
•	M_Map plotting package: https://www.eoas.ubc.ca/~rich/map.html
•	TTCmap.m colormaps from topotoolbox: https://github.com/wschwanghart/topotoolbox/blob/master/colormaps/ttcmap.m 

Calling SCION_initialise(0) runs the model and plots full output. 
Calling SCION_initialise(-1) runs model and plots only fluxes for brevity. 
Calling SCION_initialise(-2) runs model for fixed present day forcings, use to check the present day steady state if modifying the model. Note that due to constant supply of carbon from the mantle and conversion to organic C, crustal organic C increases throughout all model runs.


For more information on model derivation and running the model, see the Guidebook in the documentation folder.

For tutorial videos see the [Earth Evolution Modelling Group code page](https://earthevolutionmodelling.com/code)

This code is free to use. The model is under continual revision/extension by my research group and collaborators. For any queries or collaboration ideas please email b.mills@leeds.ac.uk
