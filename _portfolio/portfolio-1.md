---
title: "Battery dispatch modelling"
excerpt: "Analysis of the equations for battery energy storage systems (BESS) and their role in optimal dispatch.<br/><img src='/images/3_bus_system.png'  width='600' height='600'>"
collection: portfolio
---

## Battery energy storage systems (BESS)

BESS holds an important role in various avenues of power systems, from complementing renewable energy sources (such as Solar) to participating in electricity markets for arbitrage and ancillary services. The basic set of linearized equations for BESS are as follows [1]

$$SOC_n = SOC_{n-1} + \tau\left(\eta P^c_n - \frac{1}{\eta}P^d_n\right)$$

where $SOC_n$ [MWh] is the state of charge at time step $n$, $P^c_n, P^d_n$ [MW] are the charge and discharge rates of the battery, respectively, $\eta$ [-] is the round trip efficiency, and $\tau$ [h] is the time step. The variables are further bounded $\forall n$

$$ P_{min} \leq P^c_n, \; P^d_n \leq P_{max}$$, 
$$ SOC_{min} \leq SOC_n \leq SOC_{max}$$.


## References
<a id="1">[1]</a> 
Pozo, David (2022). 
Linear battery models for power system analysis,
Electric Power Systems Research, 212.
