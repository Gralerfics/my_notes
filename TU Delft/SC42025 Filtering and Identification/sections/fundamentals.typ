#import "../generic.typ": *

#import "@preview/cetz:0.4.2"
#import "@preview/fletcher:0.5.8" as fletcher: diagram, node, edge

= Fundamentals

== Linear Algebra

#Cre("TODO")

// 矩阵、向量运算

// 向量空间、子空间

// 秩

// 基本子空间

// 逆、行列式、特征值、特征向量

// （半）正定

// 特征值分解，奇异值分解（代入最小二乘，注意到 UUT 为投影阵），QR分解

// linear equations; solve.

// deterministic least-squares problem; solution; numerical way.

== Probabilistic Theory and Random Process

=== Random Variables

#Cre("TODO")

// 随机变量

// 均值方差相关，及其样本估计，以及性质

对于随机变量 $theta$，其均值（mean）记为 $mu_theta$，定义及样本估计为：

$
mu_theta = E[theta] approx 1/N sum_(i=1)^N theta_i
$

其协方差（covariance）矩阵记为 $P_theta$，定义及样本估计为（#Cre("TODO") 如果是用的样本均值大概是要除以 N - 1 吧）：

$
P_theta = E[(theta - mu_theta)(theta - mu_theta)^T] approx 1/N sum_(i=1)^N (theta_i - mu_theta)(theta_i - mu_theta)^T
$

其相关（correlation）矩阵记为 $R_theta$，定义及样本估计为：

$
R_theta = E[theta theta^T] approx 1/N sum_(i=1)^N theta_i theta_i^T
$

对于 $theta in RR^n$，$P_theta$ 和 $R_theta$ 都是半正定的 $RR^(n times n)$ 矩阵。

=== Gaussian Distribution

// （多维）高斯分布，记号，σ，PDF

对于一个高斯随机变量 $theta in RR^n$，记为：

$
theta ~ cal(N)(theta; mu_theta, P_theta)
$

其概率密度函数为：

$
p(theta) = 1/sqrt((2 pi)^n abs(P_theta)) exp(-1/2 (theta - mu_theta)^T P_theta^(-1) (theta - mu_theta))
$

可验证其期望 $E[theta] = mu_theta$，其协方差 $E[(theta - mu_theta)(theta - mu_theta)^T] = P_theta$。

=== Random Processes

TODO

// 随机过程，（多过程）WSS

// 白噪声性质，记号

// 数值计算均值方差，ergodic，noiseProperties.m（？）
