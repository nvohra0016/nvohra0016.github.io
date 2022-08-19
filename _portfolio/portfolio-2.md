---
title: "Surface-Subsurface Flow Modeling"
excerpt: "Well balanced discretization of shallow water equations and their coupling with subsurface flow.<br/><img src='/images/SW_img.png' width='500' height='500'>"
collection: portfolio
---

## Shallow water equations

Surface flow can be modeled using the shallow water equations. They are used to model rivers, floods, flow in channels; in dam break problems; in oceanic modeling for tsunamis, tides etc.

<div class='wrapper' align='center'>
<section>
    <img id='gif-click' src='/images/SW2.gif'  width='600' height='600'/>
</section>
</div>
  
I have worked on the shallow water equation solver in [Amanzi](https://github.com/amanzi/amanzi). ​Amanzi is a framework which allows users to model phenomena involving coupling of flow and reactive transport. ​ 

<div align='center'>
<img src='/images/AmanziLogo.png' width='400' height='400'>
</div>
  
Here are some aspects of the numerical implementation of shallow water equations which I have worked on:
- **Well-balancedness**: a well-balanced numerical scheme is that which preserves the stationary steady-state i.e. lake at rest solutions
- **Extension to polygonal meshes**: achieved using piecewise linear interpolation of bathymetry after it is defined on cell nodes.
- **Central upwind flux**: more accurate than the *Rusanov* flux.
- **Higher order time stepping**.

<div class='wrapper' align='center'>
<section>
    <img id='gif-click' src='/images/SW.gif'  width='600' height='600'/>
</section>
</div>

- Kurganov, [*Finite-volume schemes for shallow-water equations*](https://www.semanticscholar.org/paper/Finite-volume-schemes-for-shallow-water-equations-Kurganov/0919165796fb0929c8fd316d159d1ad2b39664a7), 2018
- Beljadid, Mohammadian, Kurganov, [*Well- balanced positivity preserving cell-vertex central-upwind scheme for shallow water flows*](https://www.infona.pl/resource/bwmeta1.element.elsevier-32db18a1-2c8c-3ba6-8723-6cfa3e616926), 2016

## Surface-subsurface coupling

The coupling between the shallow water equations and the Richards equation can be achieved through subcycling.

<div align='center'>
<img src='/images/SurfaceSubsurfaceCoupling1.png' width='600' height='600'>
</div>

Through subcycling the different time scales at which the two flows occur in most practical situations is exploited. Essentially, the more physically transient nature of surface flow compared to subsurface flow is taken into account.

### Example:

An example showing 1D dam break over a fully saturated subsurface is shown below. Notice how the pressure distribution changes as the water flows over the dry subsurface.

<div align='center'>
<img src='/images/SurfaceSubsurfaceDamBreak.png' width='600' height='600'>
</div>
