#import "../generic.typ": *

#import "@preview/cetz:0.4.2"
#import "@preview/fletcher:0.5.8" as fletcher: diagram, node, edge

== Stochastic Signal Detection

#Cre("TODO") 假设，P13。

=== White Gaussian Signal Detection

// 信号和噪声都是零均值高斯随机过程，且协方差都已知，且二者独立。
// 照常 NP 搞出 x[n] 平方和作为 T，本质上是在检测能量大小来做决策。
// 最后惯例分析指标。

#Cre("TODO")

=== Estimator-Correlator // Gaussian Signals With Arbitrary Covariance Buried in White Noise

// 信号是已知协方差矩阵的高斯随机过程，噪声还是高斯白噪声。
// 继续 NP 搞出一个 T 的形式，套一个矩阵逆引理收拾好看一点，本质上还是一个二次型。
// 考虑之前 deterministic 信号中是 x^T s，这里就可以把 x^T 后面的部分视为对 s 的估计 hat(s)，整理一下得到类似 LMMSE 的形式。所以像是顺手把信号估计出来然后再套 deterministic 的那套，等价是最优的，实现中直接过一个 Wiener filter。
// 一个 example P25 起，还没看。

#Cre("TODO")

=== Estimator-Dewhitener // Gaussian Signals With Arbitrary Covariance Buried in Colored Noise

// 然后是信号还是任意已知协方差矩阵，但噪声也是已知协方差的有色噪声了。
// 和前面形式是类似的，流程也类似。
// 然后类似 deterministic 的情况，最后可以把中间分解，表示 whitening 过程。

#Cre("TODO")

=== General Gaussian Detection

// 完全统一起来，信号不再零均值，还要有起伏了。

#Cre("TODO")

$
&cal(H)_0: quad bold(x) = bold(w) ~ cal(N)(bold(0), C_w) \
&cal(H)_1: quad bold(x) = bold(s) + bold(w) ~ cal(N)(bold(mu)_s, C_s + C_w)
$

// 还是 NP，算出最终完全版形式，可以分别退化到之前的问题。

#Cre("TODO")
