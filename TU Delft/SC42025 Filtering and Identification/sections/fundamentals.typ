#import "../generic.typ": *

#import "@preview/cetz:0.4.2"
#import "@preview/fletcher:0.5.8" as fletcher: diagram, node, edge

= Fundamentals

== Linear Algebra

#Cre("TODO") 关于舒尔补。

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

对于随机变量 $theta$，其均值（mean）记为 $mu_theta$，定义及样本估计为：

$
mu_theta = E[theta] approx 1/N sum_(i=1)^N theta_i
$

其协方差（covariance）矩阵记为 $P_theta$，定义及样本估计为（此处 $mu_theta$ 为均值真值而非样本均值，故前系数为 $1/N$ 而不需要 $1/(N-1)$）：

$
P_theta = E[(theta - mu_theta)(theta - mu_theta)^T] approx 1/N sum_(i=1)^N (theta_i - mu_theta)(theta_i - mu_theta)^T
$

其相关（correlation）矩阵记为 $R_theta$，定义及样本估计为：

$
R_theta = E[theta theta^T] approx 1/N sum_(i=1)^N theta_i theta_i^T
$

对于 $theta in RR^n$，$P_theta$ 和 $R_theta$ 都是半正定的 $RR^(n times n)$ 矩阵。

=== Gaussian Distribution

对于一个高斯随机变量 $theta in RR^n$，记为：

$
theta ~ cal(N)(theta; mu_theta, P_theta)
$

其概率密度函数为：

$
p(theta) = 1/sqrt((2 pi)^n abs(P_theta)) exp(-1/2 (theta - mu_theta)^T P_theta^(-1) (theta - mu_theta))
$

可验证其期望 $E[theta] = mu_theta$，其协方差 $E[(theta - mu_theta)(theta - mu_theta)^T] = P_theta$。

#Cre("TODO") 联合高斯、边缘分布、一阶/二阶约束、相关性和独立性、反例、……
// 两个反例：各自高斯 + 二阶不相关，不能说明独立
// #image("/assets/image-19.png")
// #image("/assets/image-20.png")
// 继续问：二阶以上的特征如何具象化，具体是什么样子的存在？
// 以及 kalman 中题设是联合高斯，这才能推导出不相关 -> 独立，否则不行

#Cre("TODO") 相乘、除、线性变换、……

=== Random Processes

略。
