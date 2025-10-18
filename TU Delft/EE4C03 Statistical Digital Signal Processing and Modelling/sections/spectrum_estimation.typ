#import "../generic.typ": *

#import "@preview/cetz:0.4.2"
#import "@preview/fletcher:0.5.8" as fletcher: diagram, node, edge

= Spectrum Estimation

#text(fill: red, "（TODO）") @sec:fun_dsp_sa。

功率谱的估计显然是个很有用的东西。例如，加性噪声、信号与噪声不相关前提下的非因果 Wiener 平滑滤波器有如下频率响应：

$
H(e^(j omega)) = (P_(d x) (e^(j omega))) / (P_x (e^(j omega))) = (P_d (e^(j omega))) / (P_d (e^(j omega)) + P_v (e^(j omega)))
$

如果目标信号与噪声的功率谱密度已知则可以直接求出其频率响应，若不知道则可以通过谱估计得到结果；再例如对窄频带信号的检测与追踪等。

总之本节考虑宽平稳随机过程的功率谱密度的估计。具体地，给定一个随机过程产生的随机信号，如何根据它有效地估计这个随机过程的谱？

最直接的想法是，由 Wiener–Khinchin Theorem，我们知道自相关函数就可以通过傅里叶变换得到功率谱：

$
P_x (e^(j omega)) = sum_(k=-infinity)^infinity r_x (k) e^(j k omega)
$

而对于一个遍历的随机过程产生的随机信号 $x[n]$，我们可以这样得到该过程的自相关函数：

$
r_x (k) = lim_(N->infinity) {1/(2N+1) sum_(n=-N)^N x[n+k] x^*[n]}
$

这是由于满足遍历性假设，我们用一条无限长的实现就可以无偏地估计出该过程的统计特征。

但显然这么做存在一些问题：第一，我们拥有的样本长度往往是有限的，例如地震波等本身就很短的信号，以及语音信号等只有很短时间内才近似满足平稳性假设的信号；第二，信号样本往往自己还包含噪声。

谱估计的方法可以分为两类：一类是无参数（Nonparametric）的方法，如从估计序列的自相关入手，变换得到功率谱；另一类是在对随机过程的模型有先验了解的情况下可以使用的有参（Parametric）估计的方法，从估计模型参数入手，再由模型计算功率谱。

== Nonparametric Spectrum Estimation

=== Periodogram

本章开头提到从样本中估计自相关函数，然后傅里叶变换得到功率谱的方法。即使现在样本长度有限（例如 $N$ 个），我们也就用这点样本来直接估计自相关函数：

$
hat(r)_x (k) = 1/N sum_(n=0)^(N-1-k) x[n+k] x^*[n], quad k=0, 1, dots, N-1
$

#blockquote[
    这里求和有 $N-k$ 项，但前面除以的是 $N$，最终得到的是一个有偏估计（Biased Estimation）：

    $
    E{hat(r)_x (k)} &= 1/N sum_(n=0)^(N-1-k) E{x[n+k] x^*[n]} \
    &= 1/N sum_(n=0)^(N-1-k) r_x (k) = (N-k)/N r_x (k)
    $

    而 $|k| >= N$ 时该估计直接就是 $0$。不过 $N$ 趋于无穷时就相等了，所以是渐进无偏（Asympotically Unbiased）的。

    #text(fill: red, "（TODO）")不过不太懂为何不直接换成除以 $N-k$，或许有它的道理。
]

此处实际上是假定了 $x[n]$ 在 $[0, N-1]$ 以外的部分值都为 $0$，等效于将求和范围设置为只从 $0$ 到 $N-1-k$。接下来再对它使用傅里叶变换估计功率谱，称为周期图（Periodogram）：

$
hat(P)_"per" (e^(j omega)) = sum_(k=-N+1)^(N-1) hat(r)_x (k) e^(-j k omega)
$

以上这个两步走的过程主要是顺着定义来，了解原理后我们可以化简一下直接通过随机信号 $x[n]$ 求得周期图。首先我们对 $x[n]$ 在 $[0, N-1]$ 以外的部分值都为 $0$ 的假设可以用乘以一个 Dirichlet 核（@equ:fun_dsp_sa_dirichlet_kernel）的过程替代，即 $x_N [n] = x[n] w_R [n]$，其傅里叶变换为：

#emphasis_equbox([
$
X_N (e^(j omega)) = sum_(n=-infinity)^infinity x_N [n] e^(-j omega n) = sum_(n=0)^(N-1) x[n] e^(-j omega n)
$
])

由：

$
hat(r)_x (k) = 1/N sum_(n=-infinity)^infinity x_N [n+k] x_N^* [n] = 1/N x_N [k] * x_N^* [-k]
$

再由 DTFT 的性质，可得功率谱估计为：

#emphasis_equbox([
$
hat(P)_"per" (e^(j omega)) = 1/N X_N (e^(j omega)) X_N^* (e^(j omega)) = 1/N abs(X_N (e^(j omega)))^2
$
])

==== An Equivalent Perspective from Filter Banks

#text(fill: red, "（TODO）")书 396 页。

==== Performance of the Periodogram

可以验证在样本长度 $N -> infinity$ 时周期图是功率谱密度的渐进（只是渐进）无偏估计，即：

$
lim_(N->infinity) E{hat(P)_"per" (e^(j omega))} = P_x (e^(j omega)) \
lim_(N->infinity) "Var"{hat(P)_"per" (e^(j omega))} = 0
$

#text(fill: red, "（TODO）")书 8.2.2 节。

=== Modified Periodogram

=== Averaged Periodogram

== Minimum Variance Spectrum Estimation

== Parametric Spectrum Estimation

=== For Autoregressive Moving Average Models
