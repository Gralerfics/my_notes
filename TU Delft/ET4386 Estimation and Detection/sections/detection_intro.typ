#import "../generic.typ": *

#import "@preview/cetz:0.4.2"
#import "@preview/fletcher:0.5.8" as fletcher: diagram, node, edge

= Detection Theory

== Introduction

对于一组数据 $bold(x) = [x[0], x[1], dots, x[N-1]]^T$，设计一个测试函数 $T$，然后依 $T(bold(x))$ 的值作出决策。

#Cre("TODO")

// 假设检验；binary detection；一般零假设代表不存在，备择假设代表存在；H0: x[n] = w[n]，H1: x[n] = s[n] + w[n]

// 什么决定我们可以多好地区分两种假设？什么是好的阈值？改变阈值的影响？错误决策的影响（miss？false alarm？）？

// 例子：H0: x[n] = w[n], n = 0, 1, ..., N - 1 和 H1: x[n] = 1 + w[n], n = ...；选择 T(x[n]) = 1 / N sum_(n=0)^(N-1) x[n] >= γ

// detection performance：永远正确的决策？无法永远正确，怎么尽量最优？要找最优，用什么指标？

定义挠度系数（deflection coefficient）为：

$
d^2 = [EE (T; cal(H)_1) - EE (T; cal(H)_0)]^2 / "var"(T; cal(H)_0)
$ <eq:d_intro_defletion_coeff>

从形式上看，挠度系数包含了两方面观察，一是两类数据应用 $T$ 映射后落点中心的距离（越远越好），二是零假设数据映射后方差的大小（越小越好，越不易交叠）。

对于不相关（uncorrelated）高斯数据有 $d^2 = (N A^2) / sigma^2$，也对应了提高 SNR（$A^2 / sigma^2$）或 $N$ 都可以提升检测表现。

=== Important Probability Density Functions

==== Gaussian

概率密度函数（probability density function，PDF）为：

$
p(x) = 1 / sqrt(2 pi sigma^2) exp(-1 / (2 sigma^2) (x - mu)^2), quad -infinity < x < +infinity
$

其中 $mu$ 为均值，$sigma^2$ 为方差。若 $mu = 0$，$sigma^2 = 1$，则称标准正态分布。

一个标准正态分布的累积分布函数（cumulative distribution function，CDF）为：

$
Phi(x) = integral_(-infinity)^x 1 / sqrt(2 pi) exp(-1/2 t^2) dif t
$

我们会经常使用：

$
Q(x) = 1 - Phi(x) = integral_x^infinity 1 / sqrt(2 pi) exp(-1/2 t^2) dif t
$

称为 Q-函数。举个使用的例子，对于：

$
T(x[n]) = 1 / N sum_(n=0)^(N-1) x[n], quad x ~ cal(N)(mu, sigma^2)
$

计算 $Pr(T(x[n]) >= gamma)$，结果可用 Q-函数表示如下：

$
Pr(T(x[n]) >= gamma) = Pr(1 / N sum_(n=0)^(N-1) x[n] >= gamma) = Q((gamma - mu)/sqrt(sigma^2 \/ N))
$

解释一下，$1 / N sum_(n=0)^(N-1) x[n]$ 是多个高斯分布的平均，其服从高斯分布 $T(x[n]) ~ cal(N)(mu, sigma^2 \/ N)$，我们可以两侧标准化，因为是右尾概率（right-tail probability），即求大于或等于某个常数的概率，直接套入针对标准正态分布的 Q-函数：

$
Pr(T(x[n]) >= gamma) = Pr((T(x[n]) - mu)/sqrt(sigma^2 \/ N) >= (gamma - mu)/sqrt(sigma^2 \/ N)) = Q((gamma - mu)/sqrt(sigma^2 \/ N))
$

==== Central Chi-Squared

卡方分布来自多个*独立*的*标准正态分布*随机变量的平方和，即：

$
x = sum_(i=1)^v x_i^2, quad x_i ~ cal(N)(0, 1)
$

自由度为 $v$ 的卡方分布的概率密度函数为：

$
p(x) = cases(
    1 / (2^(v \/ 2) Gamma(v \/ 2)) x^(v \/ 2 - 1) exp(-1 / 2 x)\, quad &x > 0,
    0\, &x < 0
)
$

记为 $x ~ chi_v^2$，其中 $v$ 为整数且 $v >= 1$，$Gamma(u)$ 为伽马函数（Gamma function）：

$
Gamma(u) = integral_0^infinity t^(u-1) exp(-t) dif t
$

类似 Q-函数，一个 $chi_v^2$ 随机变量的右尾概率记为：

$
Q_(chi_v^2) (x) = integral_x^infinity p(t) dif t
$

==== F

#Cre("TODO") 在 Reading Task 里，需要吗？

==== Rayleigh

#Cre("TODO") 在 Reading Task 里，需要吗？

=== Basics of Optimal Binary Detection

// 多数时候我们设定阈值判断是否拒绝或接受假设。

// 如何作出最优的决策？有多种策略。

==== Neyman-Pearson (NP) Theorem

===== Motivation

// 给定一个 false alarm 概率，最大化检测可能性

// 例子：简单的 binary detection，A = 1，sigma^2 = 1，N = 1。

// 零假设 mu = 0，备择假设 mu = 1；检测规则小于 gamma 接受零假设，反之拒绝。

// 两类错误，第一类 false alarm：错误警报（假阳性），即实际是 0，但根据规则检测出了 1；第二类 miss：遗漏（假阴性），即实际是 1，但根据规则检测成了 0。

检测表现通常用两个因素衡量，虚警率（在负样本总体中被错检测成正样本的概率）以及检出率（在正样本总体中被正确检测为正样本的概率）。在这个例子中，虚警率是：

$
P_"FA" = "Pr"(cal(H)_1; cal(H)_0) = "Pr"(x[n] >= gamma; cal(H)_0) = integral_gamma^infinity 1 / sqrt(2 pi) exp(-t^2 / 2) dif t = Q(gamma)
$

检出率是：

$
P_"D" &= "Pr"(cal(H)_1; cal(H)_1) = integral_gamma^infinity 1 / sqrt(2 pi) exp(-(t-1)^2 / 2) dif t = Q(gamma - 1) \
&= 1 - P_"M" = 1 - "Pr"(cal(H)_0; cal(H)_1)
$

其中 $P_"M" = "Pr"(cal(H)_0; cal(H)_1)$ 是漏检（miss）的概率，即在正样本总体中错检成负样本的概率，与 $P_"D"$ 是对立的。

// 图

Neyman-Pearson 法的思路是令 $P_"FA"$ #underline[固定低于]一个较小值，同时#underline[最大化]检出率 $P_"D"$。具体地，给定数据集合 $bold(x) = [x[0], x[1], dots, x[N-1]]^T$，设置检测问题：

$
cal(H)_0: quad T(bold(x)) < gamma \
cal(H)_1: quad T(bold(x)) > gamma
$

要求设计一个 $T(dot)$ 以最大化 $P_"D"$，同时保持 $P_"FA" < alpha$，并决定阈值 $gamma$。

===== Likelihood Ratio Test (LRT)

定义似然比（likelihood ratio）为：

$
L(bold(x)) = p(bold(x); cal(H)_1) / p(bold(x); cal(H)_0)
$

其中，$p(bold(x); cal(H)_1)$ 也就是联合条件分布 $p(x[0], x[1], dots, x[N - 1]; cal(H)_1)$，代表在 $cal(H)_1$ 成立的条件下各个数据的值的分布情况，$p(bold(x); cal(H)_0)$ 类似。

若这两个分布已知，则似然比只是一个已知的用于观察的函数。直观上看，$L(bold(x))$ 代表了 #underline[“观测数据 $bold(x)$ 更像是从哪个假设成立的情况下来的”]。

===== Optimal Detector Design Using NP Theorem

根据对似然比的理解，我们自然地采用如下决策规则：

$
L(bold(x)) > lambda quad => quad cal(H)_1
$ <eq:d_intro_lrt_dr>

意思就是 “当观测数据 $bold(x)$ 足够像是从 $cal(H)_1$ 来的时候，就判为 $cal(H)_1$”。其中阈值 $lambda$ 不是随便选的，而是要服从约束 $P_"FA" = alpha$，即：

$
P_"FA" = integral_(R_1) p(bold(x); cal(H)_0) dif bold(x) = alpha
$

其中，#underline[$R_1$ 代表数据组成的空间（例如 $RR^N$）中被判据 @eq:d_intro_lrt_dr 判为 $cal(H)_1$ 的区域]，也就是 $R_1 = {bold(x): L(bold(x)) > lambda}$。

首先我们要认识到，我们不是在讨论具体的检测问题，而是通用的方法，所以一切都是在计算概率。于是这里的积分含义就是把 $cal(H)_0$ 条件分布下被我们判成 $cal(H)_1$ 的那部分所占概率加起来，代表的就是（对于我们可能遇到的任意数据，应用这一套决策方式总体意义上的）虚警率。

#blockquote([
    *推导*：

    我们欲解决的是一个带约束的优化问题：

    $
    max(P_"D") quad s.t." "P_"FA" = alpha
    $

    应用拉格朗日乘子法，定义目标函数：

    $
    F &= P_"D" + lambda_L (P_"FA" - alpha) \
    &= integral_(R_1) p(bold(x); cal(H)_1) dif bold(x) + lambda_L (integral_(R_1) p(bold(x); cal(H)_0) dif bold(x) - alpha) \
    &= integral_(R_1) #Cpu($[p(bold(x); cal(H)_1) + lambda_L p(bold(x); cal(H)_0)]$) dif bold(x) - lambda_L alpha
    $

    注意这里的 $lambda_L$ 是拉格朗日乘子，不是前面的 $lambda$，这里添加下标区分一下。

    接下来我们要看的就是，把某个具体的数据 $x$ 放进 $R_1$（也就是判为 $cal(H)_1$），会不会使目标函数 $F$ 变大？显然这只和积分里面的那项#Cpu("东西")有关。我们希望最大化 $F$，所以我们希望把会让它变大的 $bold(x)$ 取值放进去，即满足：
    
    $
    p(bold(x); cal(H)_1) + lambda_L p(bold(x); cal(H)_0) > 0
    $
    
    的 $bold(x)$，也就是：

    $
    p(bold(x); cal(H)_1) / p(bold(x); cal(H)_0) =: L(bold(x)) > - lambda_L
    $

    // TODO 为什么最大化 F？

    #underline[不满足的就不积分进去，也就是不放进 $R_1$，也就是说这个条件其实就是判决条件（decision rule）]。由于似然比 $L(bold(x))$ 一定大于 $0$，所以记 $lambda = - lambda_L > 0$，于是得到：

    $
    p(bold(x); cal(H)_1) / p(bold(x); cal(H)_0) =: L(bold(x)) > lambda
    $

    其中 $lambda$ 是如前所述从 $P_"FA" = alpha$ 约束里算出来的。这样就说明了满足这样形式的判决条件给出的就是 Neyman-Pearson 设置下最好的条件。

    // TODO 呃，好像还是有点疑问。
])

举一个*从高斯白噪声中检测直流信号的例子*，即：

$
&cal(H)_0: quad x[n] = w[n], &n = 0, 1, dots, N - 1 \
&cal(H)_1: quad x[n] = s[n] + w[n], quad &n = 0, 1, dots, N - 1
$ <eq:d_intro_np_exp_dc_in_wgn>

其中信号 $s[n] = A > 0$，噪声 $w[n]$ 是方差为 $sigma^2$ 的高斯白噪声（white Gaussian noise，WGN）。现在设计 NP 检测器 $cal(H)_1$ 判据：

$
L(bold(x)) = (1 / (2 pi sigma^2)^(N/2) exp[-1 / (2 sigma^2) sum_(n=0)^(N-1) (x[n] - A)^2]) / (1 / (2 pi sigma^2)^(N/2) exp[-1 / (2 sigma^2) sum_(n=0)^(N-1) x^2[n]]) > lambda
$

取对数是严格单调变换，所以对该判定条件的的两侧取自然对数是等价的，可以简化形式：

$
-1 / (2 sigma^2) sum_(n=0)^(N-1) (x[n] - A)^2 + 1 / (2 sigma^2) sum_(n=0)^(N-1) x^2[n] &> ln lambda
$

整理得：

$
T(bold(x)) := 1 / N sum_(n=0)^(N-1) x[n] > sigma^2/(N A) ln lambda + A / 2 =: lambda'
$

我们将判据化为了样本均值大于一个阈值 $lambda' := sigma^2/(N A) ln lambda + A / 2$ 的形式，按先前检测问题的定义，我们可以就把 $1 / N sum_(n=0)^(N-1) x[n] $ 当作我们设计的 $T$ 函数。由于是一系列高斯随机变量的均值，它服从：

$
T(bold(x)) ~ cases(
    cal(N)(0, sigma^2/N) quad &"under" cal(H)_0,
    cal(N)(A, sigma^2/N) &"under" cal(H)_1
)
$

接下来我们计算 $P_"FA"$：

$
P_"FA" = "Pr"(T(bold(x)) > lambda'; cal(H)_0) = "Pr"((T(bold(x)) - 0)/sqrt(sigma^2\/N) > lambda'/sqrt(sigma^2\/N); cal(H)_0) = Q(lambda'/sqrt(sigma^2\/N))
$

于是倒过来也有 $lambda' = sqrt(sigma^2\/N) dot Q^(-1)(P_"FA")$，我们可以令 $P_"FA"$ 为虚警率上限 $alpha$，计算出最优阈值。然后是 $P_"D"$：

$
P_"D" = "Pr"(T(bold(x)) > lambda'; cal(H)_1) = "Pr"((T(bold(x)) - A)/sqrt(sigma^2\/N) > (lambda' - A)/sqrt(sigma^2\/N); cal(H)_1) = Q((lambda' - A)/sqrt(sigma^2\/N))
$

按照前述 NP 定理，这应当是指定 $P_"FA"$ 条件下最大的 $P_"D"$。顺便，二者由下式相互关联：

$
P_"D" = Q(Q^(-1)(P_"FA") - underbrace(sqrt((N A^2)\/sigma^2), "Signal energy-\nto-noise ratio"))
$

由前 @eq:d_intro_defletion_coeff 我们可以计算挠度系数：

$
d^2 = [EE (T; cal(H)_1) - EE (T; cal(H)_0)]^2 / "var"(T; cal(H)_0) = (A - 0)^2/(sigma^2\/N) = (N A^2) / sigma^2
$

还可以注意到 $P_"D"$ 随挠度系数单调增长，都表征检测策略表现的好坏。

再举一个*两假设为方差不同的高斯随机过程的例子*，即有：

#Cre("TODO")

计算 $P_"FA"$：

$
P_"FA" = "Pr"(sum_(n=0)^(N-1) x^2[n] > gamma'; cal(H)_0) = "Pr"((sum_(n=0)^(N-1) x^2[n])/sigma_0^2 > gamma'/sigma_0^2; cal(H)_0) = Q_(chi_N^2)(gamma'/sigma_0^2)
$

于是有最优阈值 $gamma' = Q_(chi_N^2)^(-1)(P_"FA") sigma_0^2$。

#Cre("TODO")

*总结来说*，由于阈值都是最后倒着算的，所以 NP 定理这一套*其实主要是为了求解 $T$ 函数形式的*。

===== Generalized Likelihood Ratio Test (GLRT)

注意，似然比的计算中我们需要 $p(bold(x); cal(H)_1)$ 和 $p(bold(x); cal(H)_0)$ 是已知的，才能去计算 $L(bold(x))$。若这些分布的形态已知但参数未知，则可以考虑 GLRT。

具体地，先利用估计理论把参数估出来（例如使用最大似然估计量 MLE），然后再用代入最优参数后的分布计算似然比。

#Cre("TODO") 举例说明，见 *Lec 13 即 Detection Exercises Problem 3（P10）*。

再进一步，如果分布完全未知，可能就需要借助（机器）学习方法从数据中将其估计出来了。

==== Minimum Probability of Error

首先，我们这里*假设 $cal(H)_0$ 和 $cal(H)_1$ 的先验概率 $Pr(cal(H)_0)$ 和 $Pr(cal(H)_1)$ 是已知的*，由此定义误差概率（probability of error）：

$
P_e &= Pr(cal(H)_1)"Pr"(cal(H)_0; cal(H)_1) + Pr(cal(H)_0)"Pr"(cal(H)_1; cal(H)_0) \
&= Pr(cal(H)_1) P_"M" + Pr(cal(H)_0) P_"FA"
$

接下来的方法与上节 NP 定理不同之处在于希望最小化的东西不一样，我们限制希望最小化这个 $P_e$。$P_e$ 根据正负样本的大致比例（先验概率）相对均衡地融合了两类错误的概率，#underline[本质上是视两种错误同等重要]而所得到的对总体错误率的描述。

要最小化 $P_e$ 可以这样设计检测器（注意上下顺序）：

$
p(bold(x); cal(H)_1)/p(bold(x); cal(H)_0) > Pr(cal(H)_0)/Pr(cal(H)_1) = lambda
$ <eq:d_intro_mpe_detector>

坐标这个东西其实还是*似然比*，所以接下来证明的思路和 NP 定理也是类似的。

#blockquote([
    *推导*：

    #Cre("TODO") 思路和 NP 定理类似。
])

还是举 @eq:d_intro_np_exp_dc_in_wgn 中从高斯白噪声中检测直流信号的例子。最小错误率法需要知道两类样本出现的先验概率，于是这里我们假设 $Pr(cal(H)_0) = Pr(cal(H)_1) = 0.5$，这种两类样本均衡的情况对应的最小错误率检测器也称为最大似然检测器（maximum likelihood detector）。

根据 @eq:d_intro_mpe_detector 有判决规则：

$
(1 / (2 pi sigma^2)^(N/2) exp[-1 / (2 sigma^2) sum_(n=0)^(N-1) (x[n] - A)^2]) / (1 / (2 pi sigma^2)^(N/2) exp[-1 / (2 sigma^2) sum_(n=0)^(N-1) x^2[n]]) > 1
$

取对数并整理得到：

$
T(bold(x)) := 1 / N sum_(n=0)^(N-1) x[n] > A / 2
$

依此我们可以计算 $P_e$ 评估其效果：

$
P_e &= Pr(cal(H)_1)"Pr"(cal(H)_0; cal(H)_1) + Pr(cal(H)_0)"Pr"(cal(H)_1; cal(H)_0) \
&= 1/2 ["Pr"(cal(H)_0; cal(H)_1) + "Pr"(cal(H)_1; cal(H)_0)] \
&= 1/2 ["Pr"(T(bold(x)) < A / 2; cal(H)_1) + "Pr"(T(bold(x)) > A / 2; cal(H)_0)] \
&= 1/2 ["Pr"((T(bold(x)) - A)/sqrt(sigma^2\/N) < (A \/ 2 - A)/sqrt(sigma^2\/N); cal(H)_1) + "Pr"(T(bold(x))/sqrt(sigma^2\/N) > (A \/ 2)/sqrt(sigma^2\/N); cal(H)_0)] \
&= 1/2 [(1 - Q((A \/ 2 - A)/sqrt(sigma^2\/N))) + Q((A \/ 2)/sqrt(sigma^2\/N))] \
&= Q(sqrt((N A^2)/(4 sigma^2)))
$

计算时注意 $Q$ 函数对应的是右尾概率，分清大于号和小于号。以及记得一些运算性质，例如观察图像即可发现 $Q(-x) = 1 - Q(x)$。

==== Bayesian Detector

上节中提到错误率的定义本质上认为两类错误是同等重要的，但有时我们需要考虑不一样重要的情况。例如，疾病检测中有时假阴性比假阳性造成的后果更严重。

令 $C_(i j)$ 为将 $cal(H)_j$ 的情况判成 $cal(H)_j$ 的成本，然后最小化贝叶斯风险（Bayes risk）：

$
R = EE[C] = sum_(i=0)^i sum_(j=0)^1 C_(i j)"Pr"(cal(H)_i; cal(H)_j) Pr(cal(H)_j)
$

如果 $C_10 > C_00$ 且 $C_01 > C_11$（这两个条件可以理解为错误的成本应当比正确的要高），则该条件下最优检测器的判决条件为：

$
p(bold(x); cal(H)_1)/p(bold(x); cal(H)_0) > (C_10 - C_00)/(C_01 - C_11) Pr(cal(H)_0)/Pr(cal(H)_1) = lambda
$ <eq:d_intro_br_detector>

可以观察到，它相较上节的最小错误率检测器只是多了一项关于系数的校正。

==== Summary

#Cre("TODO")

总结来说，三种检测器的决策规则使用的都是相同的统计量（即似然比），但由于优化指标的不同，需要选择不同的阈值。

// 此外，检测器的设计和评估也可以大致划分为几个步骤：
// + 设计 $T$ 函数；
// + 确定 $T(bold(x))$ 在不同假设下的分布；
// + 计算最优阈值及其和 $P_"FA"$ 的关系；
// + 计算 $P_"D"$ 和 $P_"FA"$ 的关系；
// + 研究 $P_"D"$。
