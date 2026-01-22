#import "../generic.typ": *

#import "@preview/cetz:0.4.2"
#import "@preview/fletcher:0.5.8" as fletcher: diagram, node, edge

== Practical Estimators

由前面的讨论我们有了一个基本的认识：#underline[想要得到 MVU 估计量需要有关于真实模型（概率分布）的知识]，这在实际问题中#underline[通常无法实现]；而且#underline[即便知道了概率分布，也不一定能找到]方法导出 MVUE。

所以在实际情况下我们会采用一些 sub-optimal 的估计量，例如后面要介绍的 MLE、BLUE 和 LS 等。在一些特定情况下，这些 sub-optimal 估计量也可能就是 MVU 估计量，或者方差收敛到 CRLB。

=== Maximum Likelihood Estimator (MLE)

根据某些定理（Kay-I，Theorem 7.1），#underline[若分布 $p(bold(x); theta)$ 满足一系列 “正则条件”]，则未知参数 $theta$ 的最大似然估计量（MLE）渐近（数据量足够大的情况下）分布为 $hat(theta) ~^a cal(N)(theta, I^(-1)(theta))$，也就是渐近无偏且方差逼近 CRLB。

MLE 高效且渐近无偏，若对于一个估计问题存在一个高效的估计量，通过最大似然法基本就能得到它。需要注意的是，#underline[求解 MLE 也需要已知概率分布的形式]。

相较通过 CRLB 定理去寻找 MVUE 的形式，MLE 直接尝试#underline[最大化对数似然函数]，也就是令其导数即得分函数 $s(bold(x); theta) = 0$。

*举例*说明：

$
x[n] = A + w[n], quad n = 0, 1, dots, N - 1, quad w[n] ~ cal(N)(0, A)
$

该例中设定噪声的方差和直流分量大小都为 $A$ 是#underline[为了更好地演示] MLE 和理论最优之间的差异，因为如果和前面一样单独用 $sigma^2$，答案的 MLE 将恰好就是 MVUE，即对所有观测计算样本平均，就无法直观体现 MLE 的渐近最优性了。

言归正传，$bold(x)$ 的 PDF 为：

$
p(bold(x); A) = 1/(2 pi A)^(N\/2) exp[-1/(2A) sum_(n=0)^(N-1) (x[n] - A)^2]
$

取对数求导计算得分函数：

$
(partial ln p(bold(x); A))/(partial A) = -N/(2A) + 1/A sum_(n=0)^(N-1) (x[n] - A) + 1/(2A^2) sum_(n=0)^(N-1) (x[n] - A)^2
$

欲求最大化似然估计量 $hat(A)$，整理并令其为零：

$
hat(A)^2 + hat(A) - 1/N sum_(n=0)^(N-1) x^2[n] = 0
$

求解并选取有效的正数解即为 MLE：

$
hat(A) = -1/2 + sqrt(1/N sum_(n=0)^(N-1) x^2[n] + 1/4)
$

可以注意到该估计量满足如下渐近性质：

$
EE(hat(A)) ->^a A quad "and" quad "var"(hat(A)) ->^a A^2/(N(A + 1\/2))
$

我们可以针对这个问题算一下 CRLB，具体过程略去，总之结果也是 $A^2/(N(A + 1\/2))$，这印证了 MLE 的渐近最优性（再次提醒，该性质是对 PDF 有一定要求的）。
// #image("/assets/image-25.png")
// #image("/assets/image-26.png")
// #image("/assets/image-27.png")
// #image("/assets/image-28.png")
// #image("/assets/image-29.png")

那么很常规地，我们再来考虑*线性高斯模型*的例子，即有：

$
p(bold(x); theta) = 1/((2 pi)^(N/2) det(C)^(1/2)) exp[-1/2 (bold(x) - bold(h) theta)^T C^(-1) (bold(x) - bold(h) theta)]
$

计算 MLE 即：

$
hat(theta) = arg min_theta [(bold(x) - bold(h) theta)^T C^(-1) (bold(x) - bold(h) theta)]
$

常规对代价函数 $J = (bold(x) - bold(h) theta)^T C^(-1) (bold(x) - bold(h) theta)$ 展开并求导再令其为零：

$
(partial J)/(partial theta) = -2 bold(h)^T C^(-1) bold(x) + 2 bold(h)^T C^(-1) bold(h) theta = 0 \
=> hat(theta) = (bold(h)^T C^(-1) bold(h))^(-1) bold(h)^T C^(-1) bold(x)
$

熟悉的最小二乘解，可以认识到对于线性高斯模型 MLE 恰好就是 MVUE。

#blockquote([
    *关于 MLE 的参数变换性质*：

    设 $hat(theta)$ 是最大化 $p(bold(x); theta)$ 的对数似然得到的 MLE，若存在一个单射（one-to-one）函数让 $alpha = g(theta)$，那么 $alpha$ 的 MLE 就不用算了，可以直接 $hat(alpha) = g(hat(theta))$。

    如果 $g(dot)$ 不是单射的，那么求解 $hat(alpha)$ 需要最大化的是某个修正似然函数：

    $
    p_T (bold(x); alpha) = max_{theta: alpha = g(theta)} p(bold(x); theta)
    $

    即如果一个 $alpha$ 值对应多个 $theta$，那么分布函数值在该点取其中最大的那个。
])

=== Best Linear Unbiased Estimator (BLUE)

接下来是最佳线性无偏估计量（BLUE），它的计算不需要完整的 PDF 信息，但需要前两阶矩。

如名称所示，BLUE 首先是一个线性估计量，所以其实我们是将估计量限定成了如下形式：

$
hat(theta) = bold(a)^T bold(x)
$

我们的要求依旧是无偏性和最小方差性。首先满足无偏性代表：

$
EE(hat(theta)) = bold(a)^T EE(bold(x)) = theta
$

#Cre("TODO") 求解；线性模型 bold(h) theta + bold(w)；高斯马尔可夫 H bold(theta) + bold(w)
// #image("/assets/image-30.png")
// #image("/assets/image-31.png")
// #image("/assets/image-32.png")
// #image("/assets/image-33.png")
// #image("/assets/image-34.png")

// 来点具体例子？
