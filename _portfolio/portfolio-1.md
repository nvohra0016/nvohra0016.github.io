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

\begin{equation}
\label{eq:BESS_constraints}
P_{min} \leq P^c_n, \; P^d_n \leq P_{max}, \; SOC_{min} \leq SOC_n \leq SOC_{max}.
\end{equation}

### Arbitrage in wholesale electricity market

We consider a basic scenario where a given battery performs arbitrage in the wholesale electricity markets: that is, it charges up when the electricity prices are low and discharges when the prices are high, thereby marking a profit. If $\theta_n$ [Rs/MWh] is the wholesale market price at time step $n$, then we wish to find the optimial solution to 

\begin{equation}
\text{minimize} \sum_{n} \left(P^c_n \theta_n - P^d_n \theta_n \right) \tau
\end{equation}

given by eq. \ref{eq:SOC}-eq. \ref{eq:BESS_constraints}.

We assume perfect foresight of the market prices, i.e., we assume prior knowledge of the wholesale market prices over the entire day at any given time. 

We consider a $1$[MWh] battery with maximum charging/discharging rates $1$ [MW], and a round trip efficiency of $\eta = 0.9$. We take 15 [min] real time market prices (RTM) from IEX [2] for the randomly chosen day of 15th February, 2025; see Fig. 1 below.

<div align="center">
<img src='/images/BESS_project_images/wholesaleprices_RTM_15022025.png' width='500' height='500'>
</div>
<div align="center">
Figure 1. Real time market prices used in the simulation.
</div>

The results are shown in Fig. 2 below.

<div align="center">
<img src='/images/BESS_project_images/example1_SOC.png' width='370' height='370'>
<img src='/images/BESS_project_images/example1_Pd_Pc.png' width='370' height='370'>
</div>
<div align="center">
Figure 2. State of charge and charge/discharge rates for the numerical example. 
</div>

The battery follows two complete cycles following the peak electricity prices during the peak demand periods $t \approx 6$[h] and $t \approx 18$[h]. The battery first charges around $t \approx 4$[h] when the prices are at the lowest and discharges when the prices are high around $t\ approx 8$[h]. Later in the day, the battery charges around $t \approx 13$[h] and discharges when the prices are highest around $t\ approx 20$[h]. In this example, there are also two small charging/discharging periods following the local maximums of the prices around $t \approx 9$[h] and $t \approx 17$[h]. 




## References
<a id="1">[1]</a> 
Pozo, David (2022),
Linear battery models for power system analysis,
Electric Power Systems Research, 212.

<a id="2">[2]</a> 
Indian Energy Exchange (retrieved in 2025),
[https://www.iexindia.com/market-data/real-time-market/market-snapshot](https://www.iexindia.com/market-data/real-time-market/market-snapshot).
