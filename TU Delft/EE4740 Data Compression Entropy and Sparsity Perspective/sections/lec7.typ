#import "../generic.typ": *

#import "@preview/cetz:0.4.2"
#import "@preview/fletcher:0.5.8" as fletcher: diagram, node, edge

= Greedy Algorithms for Sparse Recovery

换个视角看稀疏恢复问题，$bold(y) = bold(A) bold(z)$ 等价于 $bold(y) = sum_i z_i bold(a)_i$，即我们要用尽量少的 $bold(A)$ 的列（或称为 atoms）来线性组合出 $bold(y)$。

== Matching Pursuit (MP)

初始化 $bold(z) <- bold(0)$，$bold(r) <- bold(y)$

== Orthogonal Matching Pursuit (OMP)

== Compressive Sampling Matching Pursuit (CoSAMP)

== Iterative Hard Thresholding (IHT)

== Insights Into Matching

*互一致性*（mutual coherence）：

$
mu = max_(k != l) abs(bold(a)_k^top bold(a)_l)
$
