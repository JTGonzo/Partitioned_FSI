# DISCLAIMER!!
These codes were very recently uploaded (March 20th, 2024) from my system and have yet to be vetted for general distribution/use. I am currently working through commenting, cleaning and re-testing these solver which were put together approximately 1 year ago. Please allow me the next few weeks to clean up the codes so as to make them more instructive, easy to use, and bug-free. 

If the reader is still compelled to experiment I encourgae them to do so, however, `while this disclaimer is being displayed` I suggest that one exercises ***reasonable*** caution when using the partitoned solvers provided.

# Strongly Coupled Partitioned FSI Solver
<ins>***2D & 3D redbKIT-based solvers***</ins><br>
This repository contains 2D and 3D high-fidelity computational fluid-structure interaction solvers built (mainly) in MATLAB. They are primarily focused on solving two-field problems in a partitioned fashion, where the incompressible Navier-Stokes equations - written in the Arbitrary Lagrangian Eulerian (ALE) form - are coupled with the nonlinear elastodynamics structural equations. For 2D FSI problems the user has the option to choose either a linear elastic constitutive material model for the structural mechanics solver (**2D_Part_FSI_LE**) or a nonlinear St. Venant Kirchhoff material model (**2D_Part_FSI_SVK**). For 3D FSI problems, the user can additionally choose a nearly-incompressible NeoHookean constitutive material model for the structural mechanics solver (**3D_Part_FSI_NH**). 

These solvers use a Galerkin least-squares finite element strategy to spatial discretize both the fluid and structural equations. In these solvers we exclusively use higher-order P<sub>2</sub> finite elements for the structural equations and stable Taylor-Hood P<sub>2</sub>/P<sub>1</sub> mixed finite elements for the fluid equations. A second-order backward difference scheme is adopted for the temporal discretization of the fluid domian whilst the generalized-&alpha; method is used for the structural domian. 

In all cases except the 2D rigid-body solver each code wwas developed through  
 
 was largely built starting from the redbkit monolithic FSI code which was
% deconstructed to its baseline elements, reassembled and added to (as needed), and used 
% to solve specific FSI benchmarks in a Partitioned fashion for subsequent coupling analysis

The redbKIT-based solvers presented above are bare-bones deconstructions (/simplifications) of the original [redbKIT](https://github.com/redbKIT/redbKIT).


As described in detail below; while the [**redbKIT**](https://github.com/redbKIT/redbKIT) toolkit should be the primary reference point for the interested reader, it is my belief that the condensed monolithic solvers (as presented here in their simplified intermediary forms) hold their own unique value to the wider student/researcher community. In particular, for the purposes of streamlined learning, efficient tool adoption, and additional/alternative technique development. 



<ins>***2D Rigid-Body FSI solver***</ins><br>
The partitioned FSI scheme was developed at the [University of British Columbia](https://cml.mech.ubc.ca/) by the graduate students of my superviosor [Dr. Rajeev Jaiman](https://scholar.google.com/citations?user=iofAU68AAAAJ&hl=en&oi=ao). This technique uses a finite-element spatial discretization for both the fluid and structural equations. It uses stabilized streamline-upwind Petrov-Galerkin (SUPG) P<sub>m</sub>/P<sub>m-1</sub>/P<sub>m</sub> iso-parametric finite elements to satisfy the inf–sup condition. Temporal discretization is achieved via the second-order bacward difference formula in each domain. In this solver class (2D & 3D) we exclusive use a **St. Venant Kirchhoff** constitutive material model in the nonlinear elastodynamics equations. 

Further reading on this approach can be found in 
>[**[RKJ16] J. Liu, R. Jaiman, P. Gurugubelli. A stable second-order scheme for fluid–structure interaction with strong added-mass effects **, Journal of Computational Physics, 2014.](https://doi.org/10.1016/j.jcp.2014.04.020)



## Download, Installation and Running
-------

For detailed configuartion, installion and operation instructions, in particular for the redbKIT based solvers, please follow the steps contained in the brief [INSTALL.md](INSTALL.md) file provided.

## Further Clarifying Details
For absolute clarity, the [**redbKIT**](https://github.com/redbKIT/redbKIT) toolkit by itself is an ***AMAZING*** resource for (young) academics seeking an education on the finite element method and reduced-order modelling. As such I **strongly** recommend the reader fork that repository and work through redbKIT's numerous well crafted example problems. *I personally give this resource the lion's share of credit in helping me complete my graduate studies.*

Having said that, however, as is true for (almost) all **well developed / robust** software packages used to tackle complex problems in computational physics; its sophistication tends to be its own undoing for those novice students beginning their journey in the simulation sciences or those young numerical methods researchers who are focusing on a niche element of the problem/framework/technique and just want the machinery around the issue to be (first) working and (secondly) easily digestible for peace of mind, debugging, development etc. purposes. As such, the codes provided above seek to shallow the learning curve and expedite the toolkit's adoption/integration. 

In the case of my graduate studies, for example, having simplified the redbKIT solvers into just those most essential numerical elements needed for running a reasonably complex fluid-structure interaction simulation, I was **then** able to more easily/quickly deconstruct the monolthic system of equations into a partitioned framework for solving FSI problem. From this point I built out all the supporting tools needed to perform [partitioned simulations](https://github.com/JTGonzo/Partitioned_FSI) (interface mapping, iterative acceleration etc.) and then was able to **BEGIN** my research into ["*stabilization strategies for low mass-ratio partitioned FSI simulations*"](https://jtgonzo.github.io/).  


## redbKIT : a MATLAB(R) library for reduced-order modeling of parametrized PDEs

redbKIT is a MATLAB library for finite element simulation and reduced-order modeling of Partial Differential Equations developed at [EPFL](https://www.epfl.ch/) - [Chair of Modeling and Scientific Computing](http://cmcs.epfl.ch/)). 

In particular, it includes straightforward implementations of many of the algorithms presented in the companion book:

>[**[QMN16] A. Quarteroni, A. Manzoni, F. Negri. Reduced Basis Methods for Partial Differential Equations. An Introduction**, Springer, 2016.](http://www.springer.com/us/book/9783319154305#aboutBook)

For further details please visit [redbKIT website](http://redbkit.github.io/redbKIT/).

### License
-------

**redbKIT is distributed under BSD 2-clause license**

Copyright (c) 2015-2017, Ecole Polytechnique Fédérale de Lausanne (EPFL)
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

* Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


### **redbkit** Development
-------

redbKIT was developed and maintained by [`Federico Negri`](https://www.linkedin.com/in/negrifederico/). Paola Gervasio (Università degli Studi di Brescia) is gratefully acknowledged for granting the use of parts of the finite element code MLife.


## Contact
-------
If you have any questions regarding the contents of this repository or the use of these solvers, please feel to email me at <josetg@vt.edu>.
