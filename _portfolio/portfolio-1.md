---
title: "Battery dispatch modelling"
excerpt: "Analysis of the equations for battery energy storage systems (BESS) and their role in optimal dispatch.<br/><img src='/images/3_bus_system.png'  width='600' height='600'>"
collection: portfolio
---

# Battery energy storage systems (BESS) model

BESS holds an important role in various avenues of power systems, from complementing renewable energy sources (such as Solar) to participating in electricity markets for arbitrage and ancillary services. The basic set of linearized equations for BESS are as follows[^1]

\begin{equation}
\label{eq:SOC}
SOC_n = SOC_{n-1} + \tau\left(\eta P^c_n - \frac{1}{\eta}P^d_n\right)
\end{equation}

where $SOC_n$ [MWh] is the state of charge at time step $n$, $P^c_n, P^d_n$ [MW] are the charge and discharge rates of the battery, respectively, $\eta$ [-] is the round trip efficiency, and $\tau$ [h] is the time step. The variables are further bounded $\forall n$

\begin{equation}
\label{eq:BESS_constraints}
P_{min} \leq P^c_n, \; P^d_n \leq P_{max}, \; 0 < SOC_{min} \leq SOC_n \leq SOC_{max} < E_B,
\end{equation}

where $E_B$ [MWh] is the capacity of the battery.

### Example: Arbitrage in wholesale electricity market

We consider a basic scenario where a given battery performs arbitrage in the wholesale electricity markets: that is, it charges up when the electricity prices are low and discharges when the prices are high, thereby marking a profit. If $\theta_n$ [Rs/MWh] is the wholesale market price at time step $n$, then we wish to find the optimial solution to 

\begin{equation}
\label{eq:objective}
\text{minimize} \sum_{n} \left(P^c_n \theta_n - P^d_n \theta_n \right) \tau
\end{equation}

given by eq. \ref{eq:SOC}-\ref{eq:BESS_constraints}.

We assume perfect foresight of the market prices, i.e., we assume prior knowledge of the wholesale market prices over the entire day at any given time. 

We consider a $1$[MWh] battery with maximum charging/discharging rates $1$ [MW], and a round trip efficiency of $\eta = 0.9$. We take 15 [min] real time market prices (RTM) from IEX[^2] for the randomly chosen day of 15th February, 2025; see Fig. 1 below.

<div align="center">
<img src='/images/BESS_project_images/wholesaleprices_RTM_15022025.png' width='420' height='420'>
</div>
<div align="center">
Figure 1. Real time market prices used in the simulation.
</div>

<br>

**Results and discusion.** The results are shown in Fig. 2 below.

<div align="center">
<img src='/images/BESS_project_images/example1_SOC.png' width='380' height='380'>
<img src='/images/BESS_project_images/example1_Pd_Pc.png' width='380' height='380'>
</div>
<div align="center">
Figure 2. State of charge (left) and charge/discharge rates (right) for the numerical example. 
</div>

<br>

The battery follows two complete cycles following the peak electricity prices during the peak demand periods $t \approx 6$[h] and $t \approx 18$[h]. The battery first charges around $t \approx 4$[h] when the prices are at the lowest and discharges when the prices are high around $t \approx 8$[h]. Later in the day, the battery charges around $t \approx 13$[h] and discharges when the prices are highest around $t \approx 20$[h]. In this example, there are also two small charging/discharging periods following the local maximums of the prices around $t \approx 9$[h] and $t \approx 17$[h]. 

In this case the objective value comes out to be -28654.91 [Rs], indicating a net profit. 

## Effect of battery cycles

We now consider the degradation associated with battery cycling. Typically, the life of a battery (i.e., the number of cycles it is able to perform) depends on the depth of discharge (DOD) at which it is cycled, where the DOD is defined as the absolute difference between the state of charge in consecutive time steps divided by the capacity. We denote the cycles by $\alpha(DOD)$; see Fig. 3 for an illustration (adapted from[^3]).

<div align = "center">
<img src='/images/BESS_project_images/cycle_life.png' width='380' height='380'>
</div>
<div align = "center">
Figure 3. Plot showing the number of cycles $\alpha$ as a function of DOD.
</div>

<br>

If the cost of replacing the battery is $R$ [Rs], then the cost of one cycle is given by[^4]

\begin{equation}
\label{eq:cycle_cost}
C^{cyc} = \frac{R}{\alpha\left(\delta \right)}.
\end{equation}

Or, noting that the discharge rate associated with a DOD $\delta$ can be obtained from \ref{eq:SOC} as

\begin{equation}
 P^{d} = \frac{\eta E_B \delta}{\tau},
\end{equation}

we can rewrite \ref{eq:cycle_cost} as

\begin{equation}
\label{eq:cycle_cost_pd}
C^{cyc}\left(P^{d} \right) = \frac{R}{\alpha\left(\eta^{-1} {E_B}^{-1} \tau P^{d} \right)}.
\end{equation}

Taking the cost of cycling into account, our objective function \ref{eq:objective} now becomes

\begin{equation}
\label{eq:objective_cycles}
\text{minimize} \sum_{n} \left(P^c_n \theta_n - P^d_n \theta_n \right) \tau + C^{cyc}\left(P^d_n \right).
\end{equation}

The issue (slight!) with \ref{eq:objective_cycles} is that it introduces a nonlinearity in the objective function. To this end, we linearize the function $1/\alpha$ to provide an estimate of the solution to \ref{eq:objective_cycles}. Let $\beta_1, \beta_2 \approx 1/\alpha$ be two linear approximations such that $\beta_1 <= 1/\alpha  <= \beta_2; see Fig. 4 for an illustration.

<div align = "center">
<img src='/images/BESS_project_images/inverse_cycle_life.png' width='380' height='380'>
</div>

<div align = "center">
Figure 4. Plot showing the inverse of the cycles, $1/\alpha$, and its two linear approximations $\beta_1$ and $\beta_2$.
</div>

<br>

Note that the choice of $\beta_1$ and $\beta_2$ provide a bound for our objective function, and help us to under- and overestimate the net profits. We make this clear with the next example.

### Example: Arbitrage with battery degradation

We now return to our [example](example). We consider a replacement cost of $R = 2 \times 10^5$ [Rs] (see note below) and $E_B = 1$ [MWh]. We now consider the quarterly RTM prices taken from 1/8/2025 - 3/8/2025 from IEX[^2]. 

<div align = "center">
<img src='/images/BESS_project_images/soc_cycle_degradation1.png' width='500' height='500'>
</div>

<div align = "center">
Figure 5. Plot showing the wholesale price (blue) and SOC (black) profile when using the linear function $\beta_1$. Also shown is the SOC profile when no cycling costs are taken (faded black).
</div>

<div align = "center">
<img src='/images/BESS_project_images/soc_cycle_degradation2.png' width='500' height='500'>
</div>
 

<div align = "center">
Figure 5. Plot showing the wholesale price (blue) and SOC (black) profile when using the linear function $\beta_2$. Also shown is the SOC profile when no cycling costs are taken (faded black).
</div>

**Results and discussion.** The results are shown in Fig. 5 and Fig. 6. 


## Code

The code for the dispatch model has been written by the author using GAMS and Python.


## References
[^1]: D. Pozo (2022), Linear battery models for power system analysis, Electric Power Systems Research, 212.
[^2]: Indian Energy Exchange (retrieved in 2025), [https://www.iexindia.com/market-data/real-time-market/market-snapshot](https://www.iexindia.com/market-data/real-time-market/market-snapshot).
[^3]: N. Padmanabhan, M. Ahmed, K. Bhattacharya (2020), Battery Energy Storage Systems in Energy and Reserve Markets, IEEE Transactions on Power Systems, Vol. 35, No. 1. 
[^4]: Y. Shi et al. (2019), Optimal Battery Control Under Cycle Aging Mechanisms in Pay for Performance Settings, IEEE Transactions on Automatic Control, Vol. 64, No. 6.
