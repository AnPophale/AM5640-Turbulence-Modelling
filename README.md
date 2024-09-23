# AM5640-Turbulence-Modelling
FVM implementation of Eddy Viscosity models for fully developed turbulent channel flow in MATLAB

In this project, we use the Finite Volume Method (FVM) to numerically solve the governing equations for a fully developed turbulent channel flow using MATLAB. The turbulence closure is modelled using Eddy Viscosity Models including the $\kappa-\epsilon$, $\kappa-\omega$ models as well as the Reynolds Stress Model. The results are compared with Direct Numerical Simulaton (DNS) data from Kim, Moser and Mansour (1999) as well as 2D ANSYS Fluent simulations using periodic boundary conditions. MATLAB codes developed for the project can be found in the MATLAB Codes folder and the results are summarized in the report. Following is a brief explanation of the governinng equations, details of the turbulence models used, as well as the FVM discretized equations.

**Governing Equations:**  
For a fully developed turbulent channel flow where, x is the stream wise direction, y is the wall normal direction and z is the spanwise direction, we have the following conditions:
* Statistically stationary: $\frac{\partial \bar{\phi}}{\partial t} = 0$ for any mean quantity $\bar{\phi}$
* Statistically homogeneous in the span-wise direction: $\frac{\partial \bar{\phi}}{\partial z} = 0$ for any mean quantity $\bar{\phi}$

Applying these conditions, the 3D Reynolds Averaged Navier Stokes (RANS) equations with the Boussinesq approximation along with the continuity equation reduce to:  
```math
\frac{\partial U}{\partial x} + \frac{\partial V}{\partial y} = 0 \\
```

```math
U \frac{\partial U}{\partial x} + V \frac{\partial U}{\partial y} = -\frac{1}{\rho} \frac{\partial P}{\partial x} 
+ \frac{\partial}{\partial x} \left( (\nu + \nu_t) \frac{\partial U}{\partial x} \right) 
+ \frac{\partial}{\partial y} \left( (\nu + \nu_t) \frac{\partial U}{\partial y} \right) \\
```  

```math
U \frac{\partial V}{\partial x} + V \frac{\partial V}{\partial y} = -\frac{1}{\rho} \frac{\partial P}{\partial y} 
+ \frac{\partial}{\partial x} \left( (\nu + \nu_t) \frac{\partial V}{\partial x} \right) 
+ \frac{\partial}{\partial y} \left( (\nu + \nu_t) \frac{\partial V}{\partial y} \right)
```
 
Further, for a fully developed flow, $\frac{\partial \bar{U}}{\partial x} = 0$  and $V = W = 0$
Hence, the equations reduce to a single equation of the form  
```math
0 = -\frac{1}{\rho} \frac{\partial P}{\partial x} + \frac{\partial}{\partial y} \left( (\nu + \nu_t) \frac{\partial U}{\partial y} \right)
```

**Turbulence modelling:**
Consider the RANS equations in the tensor form using index notation given below

The term which is also called as Reynolds Stresses need to be modelled as they are unknown and this is known as the turbulence closure model. Based on the modelling approach used, we have different types of turbulence closure models.

**Eddy Viscosity models:**  
In this class of models, the Reynolds Stresses are modelled using the Boussinesq approximation which is given as 
Hence, the 6 unknown Reynolds Stresses are reduced to 2 unknowns, k and vt

This leads to a general form of RANS equations for Eddy Viscosity based models 

Based on the approaches used to model the unknown turbulent viscosity, there are various types of Eddy viscosity models
We focus on two equations models such as the $\kappa-\epsilon$ and $\kappa-\omega$


**Reynolds Stress Model:**  

**FVM Discretization:**


**References:**  
[1] R. D. Moser, J. Kim, and N. N. Mansour. Direct numerical simulation of turbulent channel flow up to Reτ =590. Physics of Fluids, 11(4):943–945, 1999.  
[2] H. Versteeg and W. Malalasekera. An Introduction to Computational Fluid Dynamics - The Finite Volume Method. Longman Scientific & Technical, Harlow, England, 1st edition, 1995.

