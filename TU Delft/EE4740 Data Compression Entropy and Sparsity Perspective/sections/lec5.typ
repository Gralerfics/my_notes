#import "../generic.typ": *

#import "@preview/cetz:0.4.2"
#import "@preview/fletcher:0.5.8" as fletcher: diagram, node, edge

= Convex Relaxation for Sparse Recovery in CS

== $cal(l)_0$ Norm Minimization for Sparse Recovery in CS

$cal(l)_p$ 范数定义为：

$
norm(bold(x))_p = (sum_(i=1)^N abs(x_i)^p)^(1/p), quad p > 0
$

特别地，$cal(l)_0$ 范数定义为 $bold(x)$ 中非零项的个数。

对于给定的 $bold(y) = bold(A) bold(x)$ 以及 "$bold(x)$ 稀疏" 的先验，我们的目标是尽可能还原出一个最稀疏即非零项最少的 $hat(bold(x))$。前述 $cal(l)_0$ 范数定义上刚好就是统计向量中的非零项个数，所以可以直接写为优化问题：

$
hat(bold(x)) = arg min_(bold(z)) norm(bold(z))_0 quad "s.t." bold(y) = bold(A) bold(z)
$

但问题是 $cal(l)_0$ 范数是非凸的，导致该优化问题不是一个凸优化问题，求解困难：先要从小到大遍历假设 $s$ 稀疏成立，对于每个 $s$ 最坏还要扫描所有 $vec(N, s)$ 种非零项分布的可能性，每种情况下还要求解相应的方程。

== $cal(l)_1$ Norm Minimization for Sparse Recovery in CS

考虑凸松弛，用 $cal(l)_1$ 范数替代，因为它是最接近 $cal(l)_0$ 范数但又有凸性的范数。具体地，其他 $0<p<1$ 的范数也都是非凸函数，不方便用在这。现在问题写为：

$
bold(x)^\# = arg min_(bold(z)) norm(bold(z))_1 quad "s.t." bold(y) = bold(A) bold(z)
$

该凸优化问题可以用*投影子梯度法*（projected subgradient method）求解。迭代 $t: 1 -> T$ 步，对于第 $t$ 步，先子梯度下降：

$
bold(p)_(t+1) = bold(z)_t - eta_t bold(g), quad forall bold(g) in underbrace(bold(partial) norm(bold(z)_t)_1, "subgradients")
$

然后投影回可行域：

$
bold(z)_(t+1)
&= bold(pi)_({bold(y) = bold(A) bold(z)}) (bold(p)_(t+1)) \
&= arg min_(bold(q)) norm(bold(p)_(t+1) - bold(q))_2 quad "s.t." bold(y) = bold(A) bold(q) \
&= bold(p)_(t+1) - bold(A)^top (bold(A) bold(A)^top)^(-1) (bold(A) bold(p)_(t+1) - bold(y))
$

若考虑逐渐衰减的步长 $eta_t = c/sqrt(t)$，则收敛速率 $~O(1/sqrt(t))$；如果要达成误差 $epsilon$ 的次优解最小所需迭代次数 $T~O(1/epsilon^2)$。总体来说收敛比较慢。

此外，对于 $p>1$ 的范数，比如可以和 $cal(l)_2$ 范数对比一下。$cal(l)_1$ 范数更鼓励稀疏解，比如说假设是二维的 $bold(z)$ 即 $N = 2$，然后 $M = 1$，即 $bold(y) = bold(A) bold(z)$ 是一条直线，那按投影子梯度法模拟一下就会发现用 $cal(l)_1$ 范数会收敛于某个尖角（在坐标轴上，是稀疏的），而用 $cal(l)_2$ 范数则会收敛于某个 sublevel 圆和直线的切点，能量更分散。

我们可以具体证明 $cal(l)_1$ 范数最小化解是稀疏的（仅当解唯一时，否则有反例，比如 $bold(A) = [1 quad 1]$）。

#Cre("TODO")
