---
title: "Neural Networks I: First Network from Scratch"
collection: talks
excerpt: "We describe the framework of neural networks and build an example from scratch using NumPy."
date: 2025-11-15
---

<script type="text/javascript" async
  src="https://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-MML-AM_CHTML">
</script>

# 1. Introduction

Neural Networks may have been around since the 1940s, but their importance and popularity has permeated the current era. Although a lot of articles and web pages talk about their applications, here we discuss the building blocks of neural networks. In particular, we are interested in their mathematical framework, and how they can be used in interpolation and forecasting. 

## 1.1. Components of a Neural Network
A basic neural network may be viewed as a black-box solver that can be "trained" to take an input $x$ and returns an output $y$. For example, suppose we have measured the outside temperature at every hour during the day, but are interested in finding the temperature at $1:45$PM, or around noon the next day. A solution to this is to train a neural network on the hourly temperature values that have been measured to predict the value at any time (down to seconds, minutes etc.). 

A basic neural network has an architecture similar to that shown in Fig. 1. Here we consider a network that takes $x \in \mathbb{R}$ as an input and gives $y \in \mathbb{R}$ as an output. The hidden layer consists of weights and biases $w_{1,i}, b_{1,i}, w_{2,i} \in \mathbb{R}, \; 1 \leq i \leq 4$, and $b_2 \in \mathbb{R}$, along with a nonlinear smooth activation function $\sigma : \mathbb{R} \rightarrow \mathbb{R}$.

<div align="center">
<img src='/images/neural_network1_blog/neural_network_diagram.png' width='500' height='500'>
</div>

<div align = "center">
 Figure 1. Illustration of a single layer neural network with 4 neurons in the hidden layer. The weights are denoted by $w$ and the biases by $b$.
</div>

<br>

In this case, the output $y$ is given by 

\begin{equation}
\label{eq:nn_ex1}
    y = \sum_{i=1}^{4} \sigma\left( w_{2,i} \sigma\left(x w_{1,j} + b_{1,i} \right) + b_{2} \right).
\end{equation}

Or, more succinctly, using matrix vector notation, we rewrite \ref{eq:nn_ex1} as

\begin{equation}
    y = \sigma \left( {w_2}^T \sigma \left(w_1 x + b_1 \right) + b_2 \right)
\end{equation}

where $w_1 = [w_{1,1} \; w_{1,2} \; w_{1,3} \; w_{1,4} ]^T$ (similar definition for $b_1$ and $w_2$), and $\sigma$ is applied component wise to the vector $w_1 x + b_1$. 

Now let us assume that we are given $N \in \mathbb{Z}$ measurements of $y$, denoted by $\hat{y}_n, \; 1 \leq n \leq N$ for corresponding inputs (known) $x_n$. Consider the following problem statement. 

*Problem statement*. How do we find the parameters $w_1, b_1, w_2$ and $b_2$ such that each $y_n$ obtained using $x_n$ from \ref{eq:nn_ex1} is `sufficiently close' to $\hat{y}_n$? 

