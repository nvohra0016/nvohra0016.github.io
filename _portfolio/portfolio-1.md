---
title: "Battery dispatch modelling"
excerpt: "Analysis of the equations for battery energy storage systems (BESS) and their role in optimal dispatch.<br/><img src='/images/3_bus_system.png'  width='600' height='600'>"
collection: portfolio
---

## Battery energy storage systems (BESS) model

BESS holds an important role in various avenues of power systems, from complementing renewable energy sources (such as Solar) to participating in electricity markets for arbitrage and ancillary services. The basic set of linearized equations for BESS are as follows [1]

\begin{equation}
\label{eq:SOC}
SOC_n = SOC_{n-1} + \tau\left(\eta P^c_n - \frac{1}{\eta}P^d_n\right)
\end{equation}

where $SOC_n$ [MWh] is the state of charge at time step $n$, $P^c_n, P^d_n$ [MW] are the charge and discharge rates of the battery, respectively, $\eta$ [-] is the round trip efficiency, and $\tau$ [h] is the time step. The variables are further bounded $\forall n$

$$ P_{min} \leq P^c_n, \; P^d_n \leq P_{max}, \; SOC_{min} \leq SOC_n \leq SOC_{max}$$.

### Case 1: Arbitrage in wholesale electricity market

We consider a basic scenario where a given battery performs arbitrage in the wholesale electricity markets: that is, it charges up when the electricity prices are low and discharges when the prices are high, thereby marking a profit. If $\theta_n$ [Rs/MWh] is the wholesale market price at time step $n$, then we wish to find the optimial solution to 

$$ \text{minimize} \sum_{n} \left(P^c_n \theta_n - P^d_n \theta_n \right) \tau $$

constrained by 

## References
<a id="1">[1]</a> 
Pozo, David (2022). 
Linear battery models for power system analysis,
Electric Power Systems Research, 212.
