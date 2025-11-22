#import "../generic.typ": *

#import "@preview/cetz:0.4.2"
#import "@preview/fletcher:0.5.8" as fletcher: diagram, node, edge

= Fundamentals

== Linear Algebra

#Cre("TODO")

// 列向量；复转置（Hermitian）

// 离散时间信号 x(n) 的表示，一种正着数 N 个，另一种从 n 开始倒数 N 个

// 范数（2、1、无穷）

// 内积，正交，Cauchy-Schwarz 不等式

// 线性无关，向量空间，基矢量（正交基，单位正交基）

// 矩阵，复转置（Hermitian），赫米特矩阵

// 秩，逆，性质，矩阵逆引理

// 行列式，和逆的关系，形式，迹

// 线性方程组，解存在性，伪逆（欠定和超定）

// 对角矩阵和分块对角矩阵，Toeplitz 和 Hankel 矩阵，酉矩阵

// Hermitian form（复二次型），（半）正/负定，相似矩阵

// 特征值和特征向量，性质

// 特征值分解

== Optimization Theory

#Cre("TODO")

// 局部/全局极小值点

// 复变情况

// 多元情况

== Notation

#Cre("TODO")



$
& "一些特殊集合：" && \
& RR & quad quad quad quad quad & "Real numbers" \
& RR^n && "Real" n"-vectors ("n times 1 "matrices)" \
& RR^(m times n) && "Real" m times n "matrices" \
& RR_+ && "Nonnegative real numbers" \
& RR_(++) && "Positive real numbers" \
& CC && "Complex numbers" \
& ZZ && "Integers" \
& SS^n && "Symmetric" n times n "matrices" \
& SS_+^n && "Symmetric positive semidefinite" n times n "matrices" \
& SS_(++)^n && "Symmetric positive definite" n times n "matrices" \
\ \

& "向量和矩阵：" && \
& lambda_i (X) && i"th largest eigenvalue of symmetric matrix" X \
\ \

& "范数和距离：" && \
& B(c, r) && "Ball with center" c "and radius" r \
\ \

& "一般不等式：" && \
& x prec.eq y && "Componentwise inequality between vectors" x "and" y \
& x prec y && "Strict componentwise inequality between vectors" x "and" y \
\ \

& "拓扑和凸分析：" && \
& bold("card") C && "Cardinality of set" C \
& bold("int") C && "interior of set" C \
& bold("relint") C && "Relative interior of set" C \
& bold("cl") C && "Closure of set" C \
& bold("bd") C && "Boundary of set" C": " bold("bd") C = bold("cl") C \\ bold("int") C \
& bold("conv") C && "Convex hull of set" C \
& bold("aff") C && "Affine hull of set" C \
& K^* && "Dual cone associated with" K \
\ \

& "函数和导数：" && \
& bold("dom") f && "Domain of function" f \
& bold("epi") f && "Epigraph of function" f \

& "" && \
&  && "" \
$

// $
// &  & quad quad quad quad quad & "" \
// &  && "" \
// $

= Introduction

== Mathematical Optimization

数学优化问题（mathematical optimization problem），或就称为优化问题（optimization problem），可以表述为如下形式：

$
&"minimize"& quad &f_0 (x) \
&"subject to"& &f_i (x) < b_i, quad i = 1, dots, m
$

其中：
- 优化变量（optimization variable）：$x = (x_1, dots, x_n)$；
- 目标函数（objective function）：$f_0: RR^n -> RR$；
- 约束函数（constraint functions）：$f_i: RR^n -> RR, i = 1, dots, m$；
- 约束边界（limits or bounds）：$b_1, dots, b_m$。

而以上优化问题的解（solution）或最优（optimal）记为 $x^*$，满足：

$
forall z "with" f_1 (z) <= b_1, dots, f_m (z) <= b_m, "we have" f_0 (z) >= f_0 (x^*)
$

即 $x^*$ 是所有满足约束条件的优化变量取值中使得目标函数取值最小的那个。

#Cre("TODO")

// 线性和非线性规划

// 凸优化；凸性和线性

=== Applications

#Cre("TODO")

// array processing：天线阵列的每个天线接收平面波，接收到的数据根据天线方位不同存在相位偏移，分别加复权重线性组合后得到整个阵列的增益（每个方向 theta 上的增益）；问题是通过选择合适的一系列权重 w 来使最后的增益图达到预期；首先约束 theta_target 方向增益为 1，并离散化角度值方便处理，然后可以考虑：1、最小化所有其他方向增益模长平方之和（最小二乘，2-norm），或 2、最小化所有其他方向增益中的最大值。

// machine learning：用超平面分隔两类数据点；超平面向两方向分别平移得到两个超平面，其距离和参数 a 的模长成反比，在满足约束（两类点分别在两个超平面两侧）情况下，两个平面的距离越远越好，即最小化 a 的模长。（就是 SVM 吧）

// portfolio optimization

// device sizing in electronic circuits

// data fitting

=== Solving Optimization Problems

一般优化问题（general optimization problem）不容易求解，求解方法往往有妥协的成分，例如计算时间长、可能无解等。例外的是，最小二乘问题、线性规划问题等往往存在非常高效的求解方法，凸优化问题（包含了前面二者）亦然。

== Least-Squares Problems

最小二乘问题（least-squares problem）不考虑约束（即 $m = 0$），

$
"minimize" quad f_0 (x) = norm(A x - b)_2^2 = inline(sum_(i=1)^k (a_i^T x - b_i)^2)
$

其中，$A in RR^(k times n), k >= n$，$a_i^T$ 是 $A$ 的行向量，$x in RR^n$ 为优化变量。

求解最小二乘问题比较简单：
+ 存在解析解：$x^* = (A^T A)^(-1) A^T b$；
+ 有大量可靠、高效的算法和软件，是成熟的技术；
+ 计算复杂度在 $O(n^2 k)$ 水平（和样本数量线性相关），结构恰当可以更低。

该类问题：
+ 很容易识别；
+ 可以适当拓展，例如加权、添加正则项（regularization terms）等：

$
sum_(i=1)^k w_i (a_i^T x - b_i)^2 quad "and" quad sum_(i=1)^k (a_i^T x - b_i)^2 + rho sum_(i=1)^n x_i^2
$

== Linear Programming Problems

线性规划问题（linear programming problem）指定目标函数和约束函数都是线性函数：

$
&"minimize"& quad &c^T x \
&"subject to"& &a_i^T x < b_i, quad i = 1, dots, m
$

解决线性规划问题：
+ 没有解析解；
+ 但是也有大量可靠、高效的算法和软件，是成熟的技术；
+ 计算复杂度在 $O(n^2 m), m >= n$ 水平（和约束数量线性相关），结构恰当可以更低。

该类问题：
+ 不太容易规范化为最小二乘问题；
+ 有一些标准化的技巧可以将包含 $l_1$- 和 $infinity$-norms 或分段线性函数等的问题，转化为线性规划问题，例如，切比雪夫估计问题（Chebyshev approximation problem）：

$
"minimize" inline(max_(i = 1, dots, k)) abs(a_i^T x - b_i)
$

可以等价为求解：

$
&"minimize"& quad &t \
&"subject to"& &a_i^T x - t <= b_i, quad i = 1, dots, k \
&& &-a_i^T x - t <= -b_i, quad i = 1, dots, k
$

实际上就是令 $t = max_(i = 1, dots, k) abs(a_i^T x - b_i)$，因为是 $max$，其中隐含的约束条件即对于所有 $i$ 都有 $t >= $，这样就把目标函数中的信息转化到了约束函数中；这里最后得到的都是线性方程，优化变量不只是 $x$ 了，还有 $t$，所以写成下面的标准线性规划问题形式会更明显一点：

$
&"minimize"& quad &vec(bold(0), 1)^T vec(x, t) \
&"subject to"& &vec(a_i, -1)^T vec(x, t) <= b_i, quad i = 1, dots, k \
&& &vec(-a_i, -1)^T vec(x, t) <= -b_i, quad i = 1, dots, k
$

== Convex Optimization Problems

凸优化问题则有如下形式：

$
&"minimize"& quad &f_0 (x) \
&"subject to"& &f_i (x) < b_i, quad i = 1, dots, m \
&"where"& &f_0, dots, f_m: RR^n -> RR "are convex"
$

目标函数和约束函数都是凸函数（convex function），指它们满足：

$
f_i (alpha x + beta y) <= alpha f_i (x) + beta f_i (y), quad forall alpha + beta = 1, alpha >= 0, beta >= 0
$

所以，最小二乘问题和线性规划问题都属于凸优化问题。

求解凸优化问题：
+ 没有解析解；
+ 也有可靠和高效的算法；
+ 计算复杂度大致在 $O(max{n^3, n^2 m, F})$ 水平，其中 $F$ 表示计算各 $f_i$ 和其一阶、二阶导数的开销。

该类问题：
+ 通常较难识别出来；
+ 存在大量技巧将各种问题转化为凸优化问题，所以也有大量问题可以神奇地用凸优化技术解决。

== Example: Illumination Problem

用 $m$ 盏灯等照亮 $n$ 个小而平坦的表面，如 @fig:opt_intro_lm_prob 所示。

#figure(
    caption: [A schematic diagram for the illumination problem],
)[
    #image(
        "../figures/illumination_problem.png",
        width: 60%
    )
] <fig:opt_intro_lm_prob>

照到小表面 $k$ 的光照强度视为均匀的，记为 $I_k$，其与各灯的功率 $p_j$ 线性相关（系数来自距离平方反比和方向余弦）：

$
I_k = sum_(j=1)^m a_(k j) p_j, quad a_(k j) = r_(k j)^(-2) max{cos theta_(k j), 0}
$

现在希望通过调整各灯的功率（有界）使得各处光强都尽量接近预期的 $I_"des"$：

$
&"minimize"& quad &inline(max_(i = 1, dots, k)) abs(log I_k - log I_"des") \
&"subject to"& &0 <= p_j <= p_max, quad i = 1, dots, m
$

=== Solutions

#Cre("TODO")

// 1、所有灯都用同一个功率 p，然后去解；

// 2、最小二乘最小化各 I_k 和 I_"des" 的均方误差，然后直接把超范围的功率 p_j 截断到可行范围内；

// 3、用加权最小二乘，在最小二乘的目标函数基础上加一项每盏灯功率和 p_max / 2 的加权均方误差，旨在让功率尽量接近可行范围的中心，通过反复调整权重来使所有功率满足约束；

// 4、用线性规划，解 L1-norm 问题；

// 显然以上都是非最优解（线性规划为何）。

// 5、凸优化求解，目标函数为 f_0 (p) = inline(max_{k = 1, dots, n}) h(I_k / I_"des")，其中 h(u) = max{u, 1/u}，约束就是 0 <= p_j <= p_max, quad j = 1, dots, m。这里选的 h(u) 是凸函数，目标函数是凸函数的最大值也是凸函数，约束也是线性的，所以这是一个凸优化问题。

// (?) exact solution obtained with effort ≈ modest factor × least-squares effort

=== Discussions

#Cre("TODO")

// 1、若添加约束：任何 10 盏灯的功率不超过总功率的一半，问题将依旧好解。

// 2、若添加约束：不超过半数的灯开启，问题将变得很难解。

// 要注意有时简单的问题和困难的问题看起来非常类似，直觉并不一定有效。

== Course Goals and Topics

课程掌握的主要目标是：
+ 识别并规范问题为凸优化问题；
+ 在实验作业中使用优化求解工具（CVX、YALMIP 等）；
+ 描述最优解，分析性能极限等。

包含如下主题：
+ 背景介绍和优化理论基础；
+ 凸集和凸函数；
+ 规范凸优化问题（LP、QP、SDP）；
+ 二阶方法（无约束和带约束优化）；
+ 一阶方法（梯度、子梯度）。
