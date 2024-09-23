# AM5640-Turbulence-Modelling
### FVM implementation of Eddy Viscosity and Reynolds Stress models for fully developed turbulent channel flow

In this project, we use the Finite Volume Method (FVM) to numerically solve the governing equations for a fully developed turbulent channel flow using MATLAB. The turbulence closure is modelled using Eddy Viscosity Models including the $\kappa-\epsilon$, $\kappa-\omega$ models as well as the Reynolds Stress Model. The results are compared with Direct Numerical Simulaton (DNS) data from Kim, Moser and Mansour (1999) as well as 2D ANSYS Fluent simulations using periodic boundary conditions. MATLAB codes developed for the project can be found in the MATLAB Codes folder and the results are summarized in the report. Following is a brief explanation of the governinng equations, details of the turbulence models used, as well as the FVM discretized equations.

**Governing Equations:**  
For a fully developed turbulent channel flow where, x is the stream wise direction, y is the wall normal direction and z is the spanwise direction, we have the following conditions:
* Statistically stationary: $\frac{\partial \bar{\phi}}{\partial t} = 0$ for any mean quantity $\bar{\phi}$
* Statistically homogeneous in the span-wise direction: $\frac{\partial \bar{\phi}}{\partial z} = 0$ for any mean quantity $\bar{\phi}$

Applying these conditions, the Reynolds Averaged Navier Stokes (RANS) equations with the Boussinesq approximation along with the continuity equation reduce to (All variables are ensemble averaged unless otherwise specified):  
```math
\frac{\partial U}{\partial x} + \frac{\partial V}{\partial y} = 0 \quad \quad (1)\\
```

```math
U \frac{\partial U}{\partial x} + V \frac{\partial U}{\partial y} = -\frac{1}{\rho} \frac{\partial P}{\partial x} 
+ \frac{\partial}{\partial x} \left( (\nu + \nu_t) \frac{\partial U}{\partial x} \right) 
+ \frac{\partial}{\partial y} \left( (\nu + \nu_t) \frac{\partial U}{\partial y} \right) \quad \quad (2)\\
```  

```math
U \frac{\partial V}{\partial x} + V \frac{\partial V}{\partial y} = -\frac{1}{\rho} \frac{\partial P}{\partial y} 
+ \frac{\partial}{\partial x} \left( (\nu + \nu_t) \frac{\partial V}{\partial x} \right) 
+ \frac{\partial}{\partial y} \left( (\nu + \nu_t) \frac{\partial V}{\partial y} \right) \quad \quad (3)
```
 
Further, for a fully developed flow, $\frac{\partial \bar{U}}{\partial x} = 0$  and $V = W = 0$.
Hence, the equations reduce to a single 1D equation of the form  
```math
0 = -\frac{1}{\rho} \frac{\partial P}{\partial x} + \frac{\partial}{\partial y} \left( (\nu + \nu_t) \frac{\partial U}{\partial y} \right) \quad \quad (4)
```

**Turbulence modelling:**  
Consider the RANS equations in the tensor form using index notation given below
```math
\frac{\partial \bar{u}_i}{\partial t} + \bar{u}_j \frac{\partial \bar{u}_i}{\partial x_j}  = -\frac{1}{\rho} \frac{\partial \bar{p}}{\partial x_i} + \nu \frac{\partial^2 \bar{u}_i}{\partial x_j \partial x_j} + \frac{\partial \overline{(u_i' u_j')}}{\partial x_j} \quad \quad (5)
```

The term $\overline{(u_i' u_j')}$  which is also called as Reynolds Stresses need to be modelled as they are unknown and this is known as the turbulence closure model. Based on the modelling approach used, we have different types of turbulence closure models.

**Eddy Viscosity models:** 
One of the approaches for this uses the Boussinesq approximation which is given below and this leads to a class of turbulence models known as eddy visocity models
```math
\overline{u_i' u_j'} = -\nu_t \left( \frac{\partial \bar{u}_i}{\partial x_j} + \frac{\partial \bar{u}_j}{\partial x_i} \right) + \frac{2}{3} k \delta_{ij} \quad \quad (6)
```
Hence, the 6 unknown Reynolds Stresses are reduced to 2 unknowns, the turbulence kinetinc energy k and the turbulent or eddy viscosity $\nu_t$. Substituting the Boussinesq approximation into the RANS equations leads to a general form for Eddy Viscosity based models which is similar to equations (2) and (3). The two unknows still need to be modelled appropriatley which gives various eddy viscosity models such as the Prandtl's one equation model or 2 equation models such as the $\kappa-\epsilon$, $\kappa-\omega$ models and in this project, we focus on the later. The drawback of such models is that the inherent anisotropic nature of turbulece is not completely captured and replacing 6 unknowns with 2 promotes isotropy in the model. 

In the $\kappa-\epsilon$ model, a model transport equation in solved for the turbulence kinetic energy $\kappa$ which is derived by using model approximations in the exact transport equation for $\kappa$. Analogus to this, a transport equation is derived for the dissipation rate of turbulence kinetic energy known as $\epsilon$. The turbulence viscosity is calculated based on the values of $\kappa$ and $\epsilon$. The equations for the $\kappa-\epsilon$ model, simplified for the case of fully developed turbulent channel flow are given as follows:
```math
\frac{\partial}{\partial y}\left[(\nu + \nu_t)\frac{\partial u}{\partial y}\right] - \frac{1}{\rho}\frac{\partial P}{\partial y} = 0 \quad \quad (7)
```

```math
\frac{\partial}{\partial y}\left[\left(\nu + \frac{\nu_t}{\sigma_k}\right)\frac{\partial k}{\partial y}\right] + P_k  - \epsilon = 0 \quad \quad (8)
```

```math
\frac{\partial}{\partial y}\left[\left(\nu + \frac{\nu_t}{\sigma_\epsilon}\right)\frac{\partial \epsilon}{\partial y}\right] + C_1\frac{\epsilon}{k} P_k - C_2\frac{\epsilon^2}{k} = 0 \quad \quad (9)
```

```math
P_k = \nu_t\left(\frac{\partial u}{\partial y}\right)^2 \quad \quad (10)
```

```math
\nu_t = C_\mu\frac{k^2}{\epsilon} \quad \quad (11)
```



**Reynolds Stress Model:**  
Another approach to model the Reynolds stresses is formulating a seperate transport equation for each of the Reynolds stresses which solves the problem of isotropy in eddy viscosity models. This class of models are known as Reynolds Stress Models.

 





**FVM Discretization:**


**References:**  
[1] R. D. Moser, J. Kim, and N. N. Mansour. Direct numerical simulation of turbulent channel flow up to Reτ =590. Physics of Fluids, 11(4):943–945, 1999.  
[2] H. Versteeg and W. Malalasekera. An Introduction to Computational Fluid Dynamics - The Finite Volume Method. Longman Scientific & Technical, Harlow, England, 1st edition, 1995.

