---
title: "[WIP]: Neural Networks I: First Network from Scratch"
collection: talks
excerpt:"We describe the framework of neural networks and build an example from scratch using NumPy."
date: 2025-11-15
---

<script type="text/javascript" async
  src="https://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-MML-AM_CHTML">
</script>

# 1. Introduction

Neural Networks may have been around since the 1940s, but their importance and popularity has permeated the current era. Although a lot of articles and web pages talk about their applications, here we discuss the building blocks of neural networks. In particular, we are interested in their mathematical framework, and how they can be used in interpolation and forecasting. 

## 1.1. Components of a Neural Network
A basic neural network may be viewed as a black-box solver that can be "trained" to take an input $x$ and returns an output $y$. For example, suppose we have measured the outside temperature at every hour during the day, but are interested in finding the temperature at $1:45$PM, or around noon the next day. A solution to this is to train a neural network on the hourly temperature values that have been measured to predict the value at any time (down to seconds, minutes etc.). 

A basic neural network has an architecture as shown in Fig. 1. Here we consider a network that takes $x \in \mathbb{R}$ as an input and gives $y \in \mathbb{R}$ as an output. The hidden layer consists of weights and biases $w_{1,i}, b_{1,i}, w_{2,i} \in \mathbb{R}, \; 1 \leq i \leq 4$, and $b_2 \in \mathbb{R}$, along with a nonlinear smooth activation function $\sigma : \mathbb{R} \rightarrow \mathbb{R}$.

<div align="center">
<img src='/images/Neural_network1_blog/neural_network_diagram.png' width='450' height='450'>
</div>

<div align = "center">
 Figure 1. Illustration of a single layer neural network with 4 neurons in the hidden layer. The weights are denoted by $w$ and the biases by $b$.
</div>

<br>

In this case, the output $y$ is given by 

\begin{equation}
    y = \sum_{i=1}^{4} \sigma\left( w_{2,i} \sigma\left(x w_{1,j} + b_{1,i} \right) + b_{2} \right
\end{equation}