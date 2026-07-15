#import "../generic.typ": *

#import "@preview/cetz:0.4.2"
#import "@preview/fletcher:0.5.8" as fletcher: diagram, node, edge

= $cal(l)_1$ Norm Minimization (Uniqueness & Algorithms)

继续讨论 $cal(l)_1$ 范数最小化的细节。首先前面提到过，$cal(l)_1$ 范数最小化的解不一定唯一，这种情况下可能解非稀疏解。

另外，$cal(l)_1$ 范数最小化也不一定能正确还原出原信号，而是依赖于 $bold(A)$ 的选取。举例说明，若 $bold(A) = [1 quad 2]$，而原稀疏信号为 $vec(2, 0)$，那么会得到测量 $bold(y) = [1 quad 2] vec(2, 0) = 2$，对应直线 $z_1 + 2 z_2 = 2$，如此得到的最优解会是 $bold(x)^\# = vec(0, 1)$。该解和原信号都满足约束，但前者对应的目标函数值更小。

== Geometric Interpretation of Null Space

分析之前先讨论一下零空间的几何解释。对于我们这里行数 $M$ 小于列数 $N$ 的矩阵 $bold(A)$ 存在至少 $N-M$ 维的零空间 $cal(N)(bold(A)) = {bold(z): bold(A) bold(z) = bold(0)}$。再考虑解空间 $cal(M)(bold(A), bold(y)) = {bold(z): bold(A) bold(z) = bold(y)}$，或者说由于在我们的问题设置中 $bold(y)$ 是测量，称为测量空间/超平面。

从几何上来说，零空间平行于测量超平面，并且经过原点。// 平行是因为对于任意 $bold(z) in cal(M)(bold(A), bold(y))$ 和 $bold(z)_0 in cal(N)(bold(A))$，有 $bold(A) (bold(z) + bold(z)_0) = bold(y)$

== Uniqueness of Solution

#Cre("TODO")

有如下*定理*，在 $RR^N$ 中，$cal(l)_1$ 范数最小化可以唯一辨识出所有 $cal(J)$ 支持的稀疏向量的充要条件是：

$
norm(bold(v)_(cal(J)))_1 < norm(bold(v)_(cal(J)^c))_1, quad forall bold(v) in cal(N) (bold(A)) \\ {0}
$

其中 $cal(J)$ 是集合，其中为稀疏分量的索引，例如 $cal(J) = {4, 7, 8}$ 表示稀疏向量的非零项只出现在第 ${4, 7, 8}$ 分量上；$bold(v)_(cal(J))$ 是 $cal(J)$ 中索引的分量拼起来的部分，$bold(v)_(cal(J)^c)$ 是剩下的。

上述条件成立也称 $bold(a)$ 对于集合 $cal(J)$ 满足*零空间性质*（null-space property，NSP）。举例说明，例如 $RR_3$ 下，$cal(J) = {1, 2}$，则条件为 $abs(v_1) + abs(v_2) < abs(v_3), forall bold(v) in cal(N) (bold(A)) \\ {0}$。

== Practical Settings

#Cre("TODO") 信号不是完全 sparse 的（not exactly sparse），但是可压缩（compressible）；系数可以不一定大量为零，快速衰减即可。

#Cre("TODO") 测量存在噪声：

条件不写严格的 $bold(y) = bold(A) bold(z)$ 了，改成 $norm(bold(y) - bold(A) bold(z))_2 <= gamma$，称为 quadratically constrained basis pursuit。

或者反过来约束系数程度 $norm(z)_1 <= tau$、最小化误差 $norm(bold(y) - bold(A) bold(z))_2$，称为 least absolute shrinkage & selection operator（LASSO）。

或者用正则项，优化 $norm(bold(y) - bold(A) bold(z))_2^2 + lambda norm(bold(z))_1$，称为 basis pursuit denoising/$cal(l)_1$ norm regularized LS。$lambda in (0, infinity)$ 时该形式可以解决 CS 问题，$lambda in [0, infinity]$ 时该问题是凸优化问题。

== Iterative Shrinkage Thresholding Algorithm (ISTA)

#Cre("TODO") 来解前面的 $cal(l)_1$ norm regularized LS。

#Cre("TODO") 子梯度法：

$
bold(p)_(t+1) &= bold(z)_t - alpha_t ((partial norm(bold(y) - bold(A) bold(z)_t)_2^2)/(partial bold(z)_t) + lambda bold(g)), quad forall bold(g) in underbrace(bold(partial) norm(bold(z)_t)_1, "subgradients") \
&= bold(z)_t - alpha_t (2 bold(A)^top (bold(A) bold(z)_t - bold(y)) + bold(g))
$

#Cre("TODO") 近端梯度下降（proximal gradient descent）到 ISTA，先在 smooth 的部分梯度下降：

$
bold(p)_(t+1) = bold(z)_t - eta dot 2 bold(A)^top (bold(A) bold(z)_t - bold(y))
$

然后在关于 $cal(l)_1$ 范数的 non-smooth 的部分做 prox：

$
bold(z)_(t+1) = S_(eta lambda) [bold(p)_(t+1)], quad
S_(eta lambda) [w_i] = cases(
    w_i - eta lambda"," quad &w_i > eta lambda,
    w_i + eta lambda"," & w_i < -eta lambda,
    0"," & abs(w_i) <= eta lambda
)
$

其步长 $eta = 1/(lambda_"max" (2 bold(A)^top bold(A)))$，收敛速率 $~O(1/t)$。
