#import "../generic.typ": *

#import "@preview/cetz:0.4.2"
#import "@preview/fletcher:0.5.8" as fletcher: diagram, node, edge

== Deterministic Signal Detection

=== In White Noise

问题是和前面例子类似的：

$
&cal(H)_0: quad x[n] = w[n], &n = 0, 1, dots, N - 1 \
&cal(H)_1: quad x[n] = s[n] + w[n], quad &n = 0, 1, dots, N - 1
$

此时我们考虑#underline[已知且 deterministic 的信号 $s[n]$]，以及高斯白噪声 $w[n] ~ cal(N)(0, sigma^2)$。设计 NP 检测器：

$
L(bold(x)) = p(bold(x); cal(H)_1) / p(bold(x); cal(H)_0) = (1 / (2 pi sigma^2)^(N/2) exp[-1 / (2 sigma^2) sum_(n=0)^(N-1) (x[n] - s[n])^2]) / (1 / (2 pi sigma^2)^(N/2) exp[-1 / (2 sigma^2) sum_(n=0)^(N-1) x^2[n]]) > lambda \
=> L(bold(x)) = exp[-1 / (2 sigma^2) (sum_(n=0)^(N-1) (x[n] - s[n])^2 - sum_(n=0)^(N-1) x^2[n])] > lambda
$

注意到检测器中 $(x[n] - s[n])^2$ 这样的结构其实在检测信号 $s[n]$ 引起的观测信号 $x[n]$ 在均值上的变化。取对数并继续整理得：

$
1/sigma^2 sum_(n=0)^(N-1) x[n]s[n] - 1/(2 sigma^2) sum_(n=0)^(N-1) s^2[n] > ln lambda
$

由于 $s[n]$ 已知，能量项 $1/(2 sigma^2) sum_(n=0)^(N-1) s^2[n]$ 可作为常数并入阈值；常系数也一并乘过去，得到：

$
T(bold(x)) := sum_(n=0)^(N-1) x[n]s[n] > lambda'
$

其中 $lambda' = sigma^2 ln lambda + 1/2 sum_(n=0)^(N-1) s^2[n]$。由此，左侧仅留下与 $bold(x)$ 有关的部分，将其作为 $T(bold(x))$。

#blockquote([
    *对此处 $T(bold(x))$ 形式的几种解读*：

    第一种，可以理解为它是在#underline[计算观测信号 $x[n]$ 与信号 $s[n]$ 的相关性]，即 correlator interpretation。观测信号与待检测信号的相关性越强，说明信号越可能存在。

    第二种，将它#underline[视为让 $x[n]$ 通过一个匹配滤波器]，即 matched filter interpretation。滤波器的单位冲激响应为 $h[n] = s[N - 1 - n]$，输入 $x[n]$ 得到的输出为：

    $
    y[n] = sum_(k=0)^n h[n-k] x[k] = sum_(k=0)^n s[N-1-(n-k)] x[k]
    $

    第 $N-1$ 个输出即 $y[N-1] = sum_(k=0)^(N-1) x[k]s[k] = T(bold(x))$。这样实现的 NP 检测器称为匹配滤波器（matched filter），其单位冲激响应来自对已知信号 $s[n]$ 的反转和平移。

    讨论以上这些是考虑在实现层面上，如果我们希望实时、物理地把这个统计量 $T(bold(x))$ 计算出来，那么实现一个匹配滤波器是更加现实有效的。

    此外，匹配滤波器在第 $N-1$ 个输出#underline[最大化信噪比（SNR）]，也即最大化挠度系数和 $P_"D"$，是 deterministic 信号 NP 检测器的#underline[最优实现]。具体一些，我们可以分析匹配滤波器的输出信噪比：

    $
    eta = (EE^2(y[N-1]; cal(H)_1))/("var"(y[N-1]; cal(H)_1))
    $

    记 $bold(s) = [s[0], dots, s[N-1]]^T$，$bold(w) = [w[0], dots, w[N-1]]^T$，$bold(h) = [h[N-1], dots, h[0]]^T$，于是有：

    $
    eta = (bold(h)^T bold(s))^2/(EE[(bold(h)^T bold(w))^2]) = (bold(h)^T bold(s))^2/(bold(h)^T EE(bold(w) bold(w)^T) bold(h)) = (bold(h)^T bold(s))^2/(bold(h)^T sigma^2 I bold(h)) = 1/sigma^2 (bold(h)^T bold(s))^2/(bold(h)^T bold(h))
    $

    由 Cauchy-Schwarz 不等式有 $(bold(h)^T bold(s))^2 <= (bold(h)^T bold(h))(bold(s)^T bold(s))$，当且仅当 $bold(h) = c bold(s)$ 时取等号，于是有：

    $
    eta <= 1/sigma^2 bold(s)^T bold(s) = cal(E)/sigma^2
    $

    取 $c = 1$，则我们证明了当 $h[N-1-n] = s[n], n = 0, 1, dots, N-1$ 时匹配滤波器将在 $N-1$ 出给出最大的输出信噪比，恰好也对应 NP 定理给出的最优 $T$ 函数形式。
])

接下来我们评估匹配滤波器的表现。

#Cre("TODO") P25 起。

$
E&(T; cal(H)_0) = ... \
E&(T; cal(H)_1) = ... \
"var"&(T; cal(H)_0) = ... \
"var"&(T; cal(H)_1) = ... \
$

对于高斯白噪声的情况，$P_"D"$ 和信号 $s[n]$ 的形状无关，只和其能量 $cal(E) = sum_(n=0)^(N-1) s^2[n]$ 有关。

=== In Colored Noise

对于有色噪声，不同时刻的噪声之间可能存在相关性，故似然函数将会变成这样的形式：

$
&p(bold(x); cal(H)_1) = 1/((2 pi)^(N/2) det^(1/2)(C)) exp[-1/2 (bold(x) - bold(s))^T C^(-1) (bold(x) - bold(s))] \
&p(bold(x); cal(H)_0) = 1/((2 pi)^(N/2) det^(1/2)(C)) exp[-1/2 bold(x)^T C^(-1) bold(x)]
$

依旧设计 NP 检测器，令似然比大于一个阈值，取对数并整理：

$
ln L(bold(x)) = bold(x)^T C^(-1) bold(s) - 1/2 bold(s)^T C^(-1) bold(s) > ln lambda \
=> T(bold(x)) := bold(x)^T C^(-1) bold(s) > ln lambda + 1/2 bold(s)^T C^(-1) bold(s) =: lambda'
$

顺便，代入 $C = sigma^2 I$ 即可退化得到同上节中一样的结果：

$
1/sigma^2 bold(x)^T bold(s) > lambda' quad => quad bold(x)^T bold(s) = sum_(n=0)^(N-1) x[n]s[n] > sigma^2 lambda'
$

接下来分析有色噪声情况下匹配滤波器的表现。

#Cre("TODO") P40 起。

// 这是求 T 在不同假设下的分布
$
E&(T; cal(H)_0) = ... \
E&(T; cal(H)_1) = ... \
"var"&(T; cal(H)_0) = ... \
"var"&(T; cal(H)_1) = ...
$

// 然后就是求 P_"FA" 和 P_"D"
$
P_"FA" &= ... \
P_"D" &= ...
$

// 上节中我们总结高斯白噪声情况下 $P_"D"$ 和信号形状无关，而现在由于 $C$ 的存在，在有色噪声的情况下信号形状将会影响指标表现。
// 顺便，同等能量下，信号长什么形状的时候 $P_"D"$ 最优？P45。

// Prewhitening transformation。正定矩阵分解的方法？

=== Summary

#Cre("TODO") P43 起，后面有总 summary 这里还要吗。
