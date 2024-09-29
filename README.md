# AM5640-Turbulence-Modelling
### FVM implementation of Eddy Viscosity and Reynolds Stress models for fully developed turbulent channel flow

In this project, we use the Finite Volume Method (FVM) to numerically solve the governing equations for a fully developed turbulent channel flow using MATLAB. The turbulence closure is modelled using Eddy Viscosity Models including the $k-\varepsilon$, $k-\omega$ models as well as the Reynolds Stress Model. The results are compared with Direct Numerical Simulaton (DNS) [data](https://github.com/AnPophale/AM5640-Turbulence-Modelling/tree/main/DNS%20Data) from Kim, Moser and Mansour (1999) as well as 2D ANSYS Fluent simulations using periodic boundary conditions. MATLAB codes developed for the project can be found in the [MATLAB Codes](https://github.com/AnPophale/AM5640-Turbulence-Modelling/tree/main/MATLAB%20codes) folder and the results are summarized in the [reports](https://github.com/AnPophale/AM5640-Turbulence-Modelling/tree/main/Reports). Following is a brief explanation of the governinng equations, details of the turbulence models used, as well as the FVM discretized equations, detailed explanations for each section can be found in [2].  

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

**Eddy Viscosity Models:**   
One of the approaches for this uses the Boussinesq approximation which is given below and this leads to a class of turbulence models known as eddy visocity models
```math
\overline{u_i' u_j'} = -\nu_t \left( \frac{\partial \bar{u}_i}{\partial x_j} + \frac{\partial \bar{u}_j}{\partial x_i} \right) + \frac{2}{3} k \delta_{ij} \quad \quad (6)
```
Hence, the 6 unknown Reynolds Stresses are reduced to 2 unknowns, the turbulence kinetinc energy k and the turbulent or eddy viscosity $\nu_t$. Substituting the Boussinesq approximation into the RANS equations leads to a general form for Eddy Viscosity based models which is similar to equations (2) and (3). The two unknows still need to be modelled appropriatley which gives various eddy viscosity models such as the Prandtl's one equation model or 2 equation models such as the $k-\varepsilon$, $k-\omega$ models and in this project, we focus on the later. The drawback of such models is that the inherent anisotropic nature of turbulece is not completely captured and replacing 6 unknowns with 2 promotes isotropy in the model. 

In the $k-\varepsilon$ model, a model transport equation in solved for the turbulence kinetic energy $k$ which is derived by using model approximations in the exact transport equation for $k$. Analogus to this, a transport equation is derived for the dissipation rate of turbulence kinetic energy known as $\varepsilon$. The turbulence viscosity is calculated based on the values of $k$ and $\varepsilon$. The equations for the $k-\varepsilon$ model, simplified for the case of fully developed turbulent channel flow are given as follows while details of other eddy viscosity models and their model assumptions can be found in [2].
```math
\frac{\partial}{\partial y}\left[(\nu + \nu_t)\frac{\partial u}{\partial y}\right] - \frac{1}{\rho}\frac{\partial P}{\partial y} = 0 \quad \quad (7)
```

```math
\frac{\partial}{\partial y}\left[\left(\nu + \frac{\nu_t}{\sigma_k}\right)\frac{\partial k}{\partial y}\right] + P_k  - \varepsilon = 0 \quad \quad (8)
```

```math
\frac{\partial}{\partial y}\left[\left(\nu + \frac{\nu_t}{\sigma_\varepsilon}\right)\frac{\partial \varepsilon}{\partial y}\right] + C_1\frac{\varepsilon}{k} P_k - C_2\frac{\varepsilon^2}{k} = 0 \quad \quad (9)
```

```math
P_k = \nu_t\left(\frac{\partial u}{\partial y}\right)^2 \quad \quad (10)
```

```math
\nu_t = C_\mu\frac{k^2}{\varepsilon} \quad \quad (11)
```

**Reynolds Stress Models:**    
Another approach to model the Reynolds stresses is formulating a seperate transport equation for each of the Reynolds stresses which solves the problem of isotropy in eddy viscosity models. This class of models are known as Reynolds Stress Models. The transport equation for the Reynolds Stresses is given as

```math
\frac{\partial}{\partial x_k} \left( \rho U_k \overline{u_i' u_j'} \right) = \mu \frac{\partial^2 \overline{u_i' u_j'}}{\partial x_k \partial x_k} + P_{ij} + \Phi_{ij} + D_{ij} - \rho \varepsilon_{ij}
```
The terms on the RHS are the viscous diffusion, production, pressure dedistribution and dissipation rate of the reynolds stresses which need further modelling. The equations for these models are complicated and hence, we do not present them here, further details can be found in [2] For example, the pressure redristributon term is split into two components, a slow term and a fast term which are modelled using the Rotta and the IP model respectively.
The terms such as need further modelling, but as the equations are complicated we do not present them here, they can be found in [2]. Along with the transport equation for the reynolds stresses, a transport equation for the disspiation of the turbulent kinetic energy is also needed to use this model. In this project, we have used wall functions for the near wall treatmnet in the reynolds stress model which is explained in the attached codes.

**FVM Discretization:**  
Here, we only describe the FVM discretization for the k epsilon model. Details of the Finite Volume Method and applications in turbulence modelling can be found in [2]. In FVM, the governing equations are integrated over a finite control volume and the domain is discretized into several such control volumes which converts the governing differential equations to a system of linear equations. This system is further solved using the Gauss Seidel method with unbder relaxation. Here, we use a central difference scheme for all diffusive terms and the non linear source terms are linearzied. Here, the superscript old refers to the values from the previous iteration which are used to decouple the equations at each iteration as well as to linearizr the source terms. All the other terms are expressed in standard FVM notation such as subscripts P, N, S refering to the parent, north and south nodes and $S_u$, $S_p$ denoting the source terms, etc. which has been followed in [2]

* Discretized u momentum equation:
```math
a_{P_u} u_P = a_{N_u} u_N + a_{S_u} u_S + Su_u
```

```math
a_{N_u} = \frac{(\nu + \nu_t^{\text{old}})_n}{\Delta y_n}, \quad 
a_{S_u} = \frac{(\nu + \nu_t^{\text{old}})_s}{\Delta y_s}, \quad 
a_{P_u} = a_{N_u} + a_{S_u}, \quad 
Su_u = \frac{\Delta y}{\rho}
```

* Discretized $k$ model equation:
```math
(a_{P_k} + Sp_k) k_P = a_{N_k} k_N + a_{S_k} k_S + Su_k
```

```math
a_{N_k} = \frac{(\nu + \nu_t^{\text{old}} / \sigma_k)_n}{\Delta y_n}, \quad 
a_{S_k} = \frac{(\nu + \nu_t^{\text{old}} / \sigma_k)_s}{\Delta y_s}, \quad 
a_{P_k} = a_{N_k} + a_{S_k}, \quad 
Sp_k = \frac{\varepsilon_P^{\text{old}}}{k_P^{\text{old}}}
```

```math
Su_k = P_{k_P} \Delta y = \left[ \nu_t^{\text{old}} \left( \frac{\partial u^{\text{old}}}{\partial y} \right)^2 \right]_P \Delta y = \nu_{t_P}^{\text{old}} \left( \frac{u_N^{\text{old}} - u_S^{\text{old}}}{\Delta y_n + \Delta y_s} \right) \Delta y
```

* Discretized $\varepsilon$ model equation:
```math
(a_{P_\varepsilon} + Sp_\varepsilon) \varepsilon_P = a_{N_\varepsilon} \varepsilon_N + a_{S_\varepsilon} \varepsilon_S + Su_\varepsilon
```

```math
a_{N_\varepsilon} = \frac{(\nu + \nu_t^{\text{old}} / \sigma_\varepsilon)_n}{\Delta y_n}, \quad 
a_{S_\varepsilon} = \frac{(\nu + \nu_t^{\text{old}} / \sigma_\varepsilon)_s}{\Delta y_s}, \quad 
a_{P_\varepsilon} = a_{N_\varepsilon} + a_{S_\varepsilon}, \quad 
Sp_\varepsilon = \frac{\varepsilon_P^{\text{old}}}{k_P^{\text{old}}}
```

```math
Su_\varepsilon = P_{\varepsilon_P} \Delta y = \left[ \nu_t^{\text{old}} \left( \frac{\partial u^{\text{old}}}{\partial y} \right)^2 \right]_P \Delta y = \nu_{t_P}^{\text{old}} \left( \frac{u_N^{\text{old}} - u_S^{\text{old}}}{\Delta y_n + \Delta y_s} \right) \Delta y
```

**Results**:  
All the results as well as inferences from the simulation can be found in the project report, here we present the comparion of velocity profiles and the turbulence kinetic energy for different models compared with the DNS data. \\
Figures 1 and 2 show the comparison of the velocity profiles using the $k-\omega$ and Reynolds Stress models.

<p align="center">
  <img src="https://github.com/user-attachments/assets/9c7977fc-bb26-4c6f-8c80-047cd652a7b5" alt="Comparison of velocity profile using κ-ω model and DNS data" style="width: 50%;">
</p>
<p align="center">
  <em>Figure 1: Comparison of velocity profile using κ-ω model and DNS data</em>
</p>

<p align="center">
  <img src="https://github.com/user-attachments/assets/8f32242e-ed64-48c1-99ae-15645ea407fc" alt="Comparison of velocity profile using Reynolds Stress model with wall function and DNS data" style="width: 50%;">
</p>
<p align="center">
  <em>Figure 2: Comparison of velocity profile using Reynolds Stress model with wall function and DNS data</em>
</p>

Figures 3 and 4 show the comparison of turbulence kinetic energy for both the models with the DNS data.

<p align="center">
  <img src="https://github.com/user-attachments/assets/8aaecaa3-a44e-425e-b097-cad17d223a6e" alt="Comparison of turbulence kinetic energy using κ-ω model and DNS data" style="width: 50%;">
</p>
<p align="center">
  <em>Figure 3: Comparison of turbulence kinetic energy using κ-ω model and DNS data</em>
</p>

<p align="center">
  <img src="https://github.com/user-attachments/assets/959fec8f-649d-48de-ac8b-84a3b29b223c" alt="Comparison of turbulence kinetic energy using Reynolds Stress model with wall function and DNS data" style="width: 50%;">
</p>
<p align="center">
  <em>Figure 4: Comparison of turbulence kinetic energy using Reynolds Stress model with wall function and DNS data</em>
</p>



**References:**  
[1] R. D. Moser, J. Kim, and N. N. Mansour. Direct numerical simulation of turbulent channel flow up to Reτ =590. Physics of Fluids, 11(4):943–945, 1999.  
[2] H. Versteeg and W. Malalasekera. An Introduction to Computational Fluid Dynamics - The Finite Volume Method. Longman Scientific & Technical, Harlow, England, 1st edition, 1995.

