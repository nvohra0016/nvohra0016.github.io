---
title: "Thermo-hydro-mechanical Analysis of Permafrost."
excerpt: "Mathematical analysis of the coupling governing subsidence of thawing permafrost by using Biot's poroelasticity equations and the Stefan problem.<br/><img src='/images/thaw_cartoon3.png'  width='500' height='500'>"
collection: portfolio
---

To study the land subsidence caused by thawing permafrost, one needs to analyze the coupling of the different mechanisms involved. Namely, these are deformation, flow, and energy. Analyzing these together, however, leads to numerically expensive complications and hence they must be studied individually and then coupled in a simpler fashion.

<div align="center">
<img src='/images/thaw_cartoon5.png' width='500' height='500'>
</div>
## Thermo-hydro-mechanical coupling

To couple flow [H], deformation[M], and energy [T], we need to consider thermo-hydro-mechanical models that take into account phase change [T] as well.

<img src='/images/ResearchSketchDiagram1.png' width='500' height='500'>

### [HM]: Biot's poroelasticity equations:
Flow [H] and deformation [M] can be modeled using Biot's poroelasticity equations given below

<img src='/images/Biot_system.png' width='500' height='300'>

### Example: 

**Terzaghi's Problem; 1D soil compaction**: Consider the compaction ​of a column of saturated soil with impermeable walls and bottom but with a free draining top (see figure on right for boundary conditions). When an external unit stress is applied on the top, the column of soil consolidates. The resulting pressure distribution can be solved for analytically using Biot's equations. The simulation below shows the compaction of a column of soil using a sinusoidal force applied at the top of the of the column. Starting from the top left, the 4 sub-simulations show, in a clockwise order,  the domain deformation, pore pressure, y-displacement and the x-displacement. 

<a href="http://www.youtube.com/watch?feature=player_embedded&v=watch?v=yGoINILFoo0
" target="_blank"><img src="http://img.youtube.com/vi/watch?v=yGoINILFoo0/0.jpg" 
alt="IMAGE ALT TEXT HERE" width="240" height="180" border="10" /></a>

### [Tp]: Energy with phase change: Stefan problem
The classical formulation of the Stefan problem allows one to solve for the domain and temperature of the different phases (eg. ice and water).
<img src='/images/Stefan_problem.png' width='500' height='500'>

By defining the enthalpy, the weak form can be derived.
<img src='/images/Stefan_weak_form.png' width='500' height='500'>
