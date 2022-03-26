---
title: "Thermo-hydro-mechanical Analysis of Permafrost."
excerpt: "Mathematical analysis of the coupling governing subsidence of thawing permafrost by using Biot's poroelasticity equations and the Stefan problem.<br/><img src='/images/thaw_cartoon3.png'  width='500' height='500'>"
collection: portfolio
---

To study the land subsidence caused by thawing permafrost, one needs to analyze the coupling of the different mechanisms involved. Namely, these are deformation, flow, and energy. Analyzing these together, however, leads to numerically expensive complications and hence they must be studied individually and then coupled in a simpler fashion.

<img src='/images/thaw_cartoon5.png' width='500' height='500'>

## Thermo-hydro-mechanical coupling

To couple flow [H], deformation[M], and energy [T], we need to consider thermo-hydro-mechanical models that take into account phase change [T] as well.

<img src='/images/ResearchSketchDiagram1.png' width='500' height='500'>

Flow [H] and deformation [H] can be modeled using Biot's poroelasticity equations.

<img src='/images/Biot_system.png' width='500' height='300'>

### Example: 
<section>
**Terzaghi's Problem; 1D soil compaction**: Consider the compaction ​of a column of saturated soil with impermeable walls and bottom but with a free draining top (see figure on right for boundary conditions). When an external unit stress is applied on the top, the column of soil consolidates. The resulting pressure distribution can be solved for analytically using Biot's equations. The simulation below shows the compaction of a column of soil using a sinusoidal force applied at the top of the of the column. Starting from the top left, the 4 sub-simulations show, in a clockwise order,  the domain deformation, pore pressure, y-displacement and the x-displacement. 
</section> 
<section>
<img src='/images/ResearchSketchDiagram1.png' width='500' height='500'>
</section>

<iframe width="420" height="315"
src="https://www.youtube.com/watch?v=yGoINILFoo0&t=2s">
</iframe>



