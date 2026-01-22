#import "../generic.typ": *

#import "@preview/cetz:0.4.2"
#import "@preview/fletcher:0.5.8" as fletcher: diagram, node, edge

== Cramér-Rao Lower Bound (CRLB)

设估计量 $hat(theta) = g(bold(x))$，其中 $bold(x) = [x[0], x[1], dots, x[N - 1]]^T$。若其为无偏估计量，则其期望应等于真值：

$
forall theta, quad EE(hat(theta)) := integral g(bold(x)) p(bold(x); theta) dif bold(x) = theta
$

对于满足一些条件的一类参数（待估计量）$theta$，接下来介绍的定理可以用于寻找 $theta$ 所有无偏估计量 $hat(theta)$ 的方差所能达到的最小值，即*下界*。

这也意味着，如果我们可以找到一个无偏估计量，且其方差等于这个下界，那么它就将是我们喜欢的最小方差无偏估计量（MVUE）。

要引出定理，首先来准备一些东西。

=== Observations

首先的首先，回到解决估计问题过程本身，我们实际上是在概率分布的似然函数上求最大值。这里的*似然函数*是指概率关于参数变化的函数，它在表达式上和概率密度函数 $p(bold(x); theta)$ 可以说是相同的，只不过我们此时认为它的*自变量是参数 $theta$* 而非概率密度函数中的 $bold(x)$。

画到图上的话就是横轴是参数值，纵轴是概率，而曲线上每点所描述的是 *“当参数取该值时，有多大概率符合当前的观测数据”*。根据这个定义，似然函数的最大值点所代表的参数取值，就是能让模型最大概率符合观测的参数取值，也即最优估计值。

顺便，我们常对似然函数取对数，得到*对数似然函数*（log-likelihood function）来进行操作。这样主要是可以放大一些极小的概率值，并利用一些较好的性质，如乘法变加法等。

不过由于概率分布或者说模型未知（不然我们也不用想办法从观测数据去估计了），似然函数也当然是*未知*的。但这不妨碍我们通过观察一些似然函数本身、梯度和二阶导数的形态来获得一些*启发*。

// #Cre("TODO") 例如

=== Score Function and Regularity Conditions

定义*得分函数*（score function）为*对数似然函数的梯度*（即对参数的一阶导数）：

$
s(bold(x); theta) = (partial ln p(bold(x); theta))/(partial theta)
$

若 $s(bold(x); theta)$ 存在且有限，并且有（允许将对参数的求导运算与对随机变量的积分运算进行交换）：

$
integral (partial p(bold(x); theta))/(partial theta) dif bold(x) = partial/(partial theta) integral p(bold(x); theta) dif bold(x)
$

则概率分布 $p(bold(x); theta)$ 符合如下正则条件（regularity condition）：

$
forall theta, quad EE[s(bold(x); theta)] = EE[(partial ln p(bold(x); theta))/(partial theta)] = 0
$ <equ:crb_regularity_cond>

// TODO 看起来本来就能交换，实则不然，具体地：
// #image("/assets/image-21.png")
// #image("/assets/image-22.png")
// #image("/assets/image-23.png")
// #image("/assets/image-24.png")

除非该 PDF 的非零定义域取决于 $theta$（#Cre("TODO")）。注意这里的期望是基于 $bold(x)$ 的分布来说的，即：

$
EE[(partial ln p(bold(x); theta))/(partial theta)] = integral (partial ln p(bold(x); theta))/(partial theta) p(bold(x); theta) dif bold(x)
$

#blockquote([
    *关于得分函数和正则条件的直观含义*：

    从定义上讲，得分函数就是对数似然函数的梯度，注意这个梯度是对参数求导。那么从物理意义上来说，它代表了我们在参数 $theta$ 上施加一个微扰时，对数似然函数值的变化趋势。

    具体一些，就是在给定观测 $bold(x)$ 的前提下，如果我们把 $theta$ 调大一点点，那么调整后的模型的合理性（对数似然函数值，或者说模型符合观测的 “概率”）是怎么变化。

    // 需要再次强调的是，现在我们是给定参数 $theta$ 和观测 $bold(x)$，……
    
    // 似然函数的自变量是参数 $theta$，所以这里的 $bold(x)$ 是固定在函数里的具体值，就是说这个函数值是在针对某个具体样本 $bold(x)$ 的，而不是针对作为随机变量的样本总体。

    于是，#underline[得分函数 $s(bold(x); theta)$ 大于 $0$ 时说明，如果观测数据是用参数为 $theta$ 的模型产生的，那么表明具有更大一些的 $theta$ 的模型会让这个观测样本 $bold(x)$ 的出现更合理，反之亦然]。

    举个简单的*例子*，考虑如下模型：
    
    $
    x = theta + w, quad w ~ cal(N)(0, 1)
    $
    
    若参数值 $theta = 0$，由于观测模型中随机噪声 $w$ 的存在，我们可能会得到的观测数据是偏离了 $theta$ 的 $x = 1$。单看这一个样本，在不知情（不知道模型）的人的角度来看，这个 $x = 1$ 从概率上来说更可能是由 $theta = 1$ 的模型生成出来的。反过来说，这个 $x = 1$ 观测的出现，本身在暗示这个 $theta = 0$ 可能调大一点会更自洽。
    
    那么下一个问题，#underline[既然我们都给定了模型参数是 $theta$，为什么由此生成的观测数据还会暗示这个 $theta$ 需要调大或调小呢]？

    因为上面的例子中的只是针对单个样本的情况，如果进行更多次的采样，我们也可能得到 $x = -1$（暗示 $theta$ 应该调小）、$x = 0$（暗示 $theta$ 不用变）等观测数据，它们的得分函数各不相同。#underline[对于一个合理自洽的问题，在观测数据足够多的情况下，观测数据整体的得分函数值应该为 $0$]，这其实就是前面的*正则条件*。
    
    正则条件的字面意思看起来就是：对任意参数值 $theta$ 来说，得分函数#underline[在随机观测 $bold(x)$ 下的期望]都是 $0$，也就是在说观测整体的得分函数平均下来应该是 $0$，不会和给定的 $theta$ 自相矛盾。

    如果还是有一点绕，就反过来看，如果正则条件不满足会怎么样？观测数据是按照给定 $theta$ 生成的，结果观测数据总体却在暗示 $theta$ 不合理，需要调大或调小，这是自相矛盾的。

    #underline[正则条件是对*模型自洽性*的要求，也是后面 Fisher 信息有意义的前提条件。]不过不满足正则条件不意味着就不合逻辑，只是接下来 Fisher 信息和 CRLB 这一套不再适用。// 例如 $x~"Uniform"(0, theta)$，参数在这里决定了分布的边界，“数据能否出现” 本身就是一种信息，且不被包含于后面的 Fisher 信息中，不能使用这一套方法讨论，也确实不满足正规条件。
])

=== Fisher Information

*在满足前述正则条件的前提下*，我们可以证明得分函数的方差等于 Fisher information $I(theta)$：

$
I(theta) = -EE[(partial^2 ln p(bold(x); theta))/(partial theta^2)] = EE[((partial ln p(bold(x); theta))/(partial theta))^2] =: EE[s(bold(x); theta)^2]
$

证明如下，对正则条件两侧求导有：

$
partial/(partial theta) EE[(partial ln p(bold(x); theta))/(partial theta)] = 0 quad => quad partial/(partial theta) integral (partial ln p(bold(x); theta))/(partial theta) p(bold(x); theta) dif bold(x) = 0
$

把导数交换进去（#Cre("TODO") 正则条件允许的？）并展开：

$
integral [(partial^2 ln p(bold(x); theta))/(partial theta^2) p(bold(x); theta) + (partial ln p(bold(x); theta))/(partial theta) (partial p(bold(x); theta))/(partial theta)] dif bold(x) = 0
$

移项并代入对数函数链式求导的结果 $(partial ln p(bold(x); theta))/(partial theta) = 1/p(bold(x); theta) (partial p(bold(x); theta))/(partial theta)$，整理得：

$
-integral (partial^2 ln p(bold(x); theta))/(partial theta^2) p(bold(x); theta) dif bold(x) = integral ((partial ln p(bold(x); theta))/(partial theta))^2 p(bold(x); theta) dif bold(x)
$

#blockquote([
    *关于 Fisher 信息的直观含义*：

    再次注意，这里 Fisher 信息的这些等价定义*依赖于正则条件的成立*。

    若 Fisher 信息较大，说明似然函数在参数真值附近曲率较大，即参数微小偏移会导致模型合理性快速下降，这使得不同参数值之间更易区分；反之若 Fisher 信息较小，则说明似然函数平坦，参数附近难以分辨，易受噪声影响。
    
    总之，#underline[这里 Fisher 信息可以大致理解为 “观测数据对参数的可辨识程度”]。
])

由定义，可以注意到 Fisher 信息有一些性质，例如*非负性*，以及*对独立观测的可加性*，即对于独立的 $x[n]$，若 $ln p(bold(x); theta) = sum_(n=0)^(N-1) ln p(x[n]; theta)$，则有：

$
-EE[(partial^2 ln p(bold(x); theta))/(partial theta^2)] = sum_(n=0)^(N-1) -EE[(partial^2 ln p(x[n]; theta))/(partial theta^2)]
$

于是进一步对于同分布的观测 $x[n]$，即有 $I(theta) = N i(theta)$，其中 $i(theta) = -EE[(partial^2 ln p(x[n]; theta))/(partial theta^2)]$。

=== CRLB Theorem

如果参数（待估计量）$theta$ 的分布 $p(bold(x); theta)$ 满足 @equ:crb_regularity_cond 的正则条件，则该参数的任何无偏估计量 $hat(theta)$ 的方差存在下界：

$
"var"(hat(theta)) >= 1 / I(theta) = 1 / (-EE[(partial^2 ln p(bold(x); theta))/(partial theta^2)]) = 1 / (EE[((partial ln p(bold(x); theta))/(partial theta))^2])
$

也就是说，#underline[在满足正则条件的前提下，任意无偏估计量的方差都不可能小于 Fisher 信息的倒数]。这符合前文对 Fisher 信息的理解：$I(theta)$ 越大，参数越易辨识，无偏估计量可能达成的方差下界越低。

如果某个无偏估计量的方差恰好等于该下界，则说明其已经充分有效地利用了观测数据中关于参数的所有信息，即可被称为*最小方差无偏估计量*。

#underline[总结一下思路]，由 CRLB 定理我们需要 Fisher 信息取倒数得到方差下界；Fisher 信息可以按定义（对数概率密度二阶导期望的相反数）来求，前面还证明了 Fisher 信息也可以用得分函数平方的期望来求（通常简单一些）；而得分函数就是对数概率密度的梯度。

#underline[更简单地]，对于一些可以化为特定形式的得分函数 $s(bold(x); theta)$ 可以省去显式求 Fisher 信息的步骤，直接得到它的表达式，详见 @sec:ed_crb_find_the_mvu_estimator。

=== Find the MVU Estimator <sec:ed_crb_find_the_mvu_estimator>

一个对#underline[任意] $theta$，方差都达到 #underline[CRLB] 的#underline[无偏]估计量可被找到，#underline[当且仅当]其得分函数可写成如下形式：

$
s(bold(x); theta) = (partial ln p(bold(x); theta))/(partial theta) = I(theta) (g(bold(x)) - theta)
$

对于一些 $g$ 和 $I$，$hat(theta) = g(bold(x))$ 就是满足 $EE(hat(theta)) = theta$ 和 $"var"(hat(theta)) = 1/I(theta)$ 的估计量，即 MVUE，这就是 CRLB 定理的取等条件。

自洽性可以通过求一下 $g(bold(x))$ 的方差验证，具体地：

$
s(bold(x); theta) &= I(theta) (g(bold(x)) - theta) \
EE[s^2(bold(x); theta)] &= EE[I(theta)^2(g(bold(x)) - theta)^2] \
I(theta) &= I(theta)^2 "var"(g(bold(x))) \
therefore "var"(g(bold(x))) &= 1 / I(theta)
$

举几个应用的例子，*例如*对于：

$
x[n] = A + w[n], quad n = 0, 1, dots, N - 1
$

其中 $w[n] ~ cal(N)(0, sigma^2)$ 是零均值高斯白噪声，即有：

$
p(bold(x); A) &= product_(n=0)^(N-1) 1/sqrt(2 pi sigma^2) exp[-(x[n] - A)^2/(2 sigma^2)] \
&= 1/(2 pi sigma^2)^(N\/2) exp[-(sum_(n=0)^(N-1) (x[n] - A)^2)/(2 sigma^2)]
$

对其取对数求导（log-likelihood）得到得分函数：

$
s(bold(x); A) = (partial ln p(bold(x); A))/(partial A) &= partial/(partial A) [-1/(2 sigma^2) sum_(n=0)^(N-1) (x[n] - A)^2] \
&= 1/sigma^2 sum_(n=0)^(N-1) (x[n] - A)
$

整理一下就能化为前述形式：

$
s(bold(x); A) = underbrace(N/sigma^2, I(theta)) (underbrace(1/N sum_(n=0)^(N-1) x[n], g(bold(x))) - underbrace(A, theta))
$

从形式中提取出各部分，即有 $hat(A) = g(bold(x))$，且 $"var"(hat(A)) >= 1/I(theta) = sigma^2/N$。

// #Cre("TODO") *再举一个例子*：伯努利分布。

=== CRLB for the General Gaussian Model

高斯模型是非常常见的信号模型，由确定模型以及和一个高斯噪声叠加而成。接下来我们讨论在高斯假设下的 CRLB，设噪声为高斯分布：

$
bold(w) ~ cal(N)(bold(0), C_w) quad <=> quad p(bold(w)) = 1/((2 pi)^(N/2) det(C_w)^(1/2)) exp[-1/2 bold(w)^T C_w^(-1) bold(w)]
$

由此定义高斯模型：

$
bold(x) = bold(h)(theta) + bold(w), quad bold(x) ~ cal(N)(bold(h(theta)), C_w)
$

或者写成概率密度函数为：

$
p(bold(x)) = 1/((2 pi)^(N/2) det(C_w)^(1/2)) exp[-1/2 (bold(x) - bold(h)(theta))^T C_w^(-1) (bold(x) - bold(h)(theta))]
$

有了分布函数，取对数并对参数求导得到得分函数：

$
s(bold(x);theta) = (partial ln p(bold(x); theta))/(partial theta) = (partial bold(h)^T (theta))/(partial theta) C_w^(-1) (bold(x) - bold(h)(theta))
$

以及二阶导：

$
(partial^2 ln p(bold(x); theta))/(partial theta^2) = (partial^2 bold(h)^T (theta))/(partial theta^2) C_w^(-1) (bold(x) - bold(h)(theta)) - (partial bold(h)^T (theta))/(partial theta) C_w^(-1) (partial bold(h) (theta))/(partial theta)
$

从而求得 Fisher 信息：

$
I(theta) = -EE[(partial^2 ln p(bold(x); theta))/(partial theta^2)] = (partial bold(h)^T (theta))/(partial theta) C_w^(-1) (partial bold(h) (theta))/(partial theta)
$

于是 CRLB 为：

$
"var"(hat(theta)) >= 1 / ((partial bold(h)^T (theta))/(partial theta) C_w^(-1) (partial bold(h) (theta))/(partial theta))
$

=== CRLB for the Linear Gaussian Model

进一步，考虑线性高斯模型：

$
bold(x) = bold(h) theta + bold(w), quad bold(w) ~ cal(N)(bold(0), C_w)
$

代入得到 CRLB：

$
"var"(hat(theta)) >= 1 / (bold(h)^T C_w^(-1) bold(h))
$

同时参考 @sec:ed_crb_find_the_mvu_estimator，将线性高斯模型的得分函数整理为如下形式：

$
s(bold(x);theta) = (partial ln p(bold(x); theta))/(partial theta) &= bold(h)^T C_w^(-1) (bold(x) - bold(h) theta) \
&= underbrace(bold(h)^T C_w^(-1) bold(h), I(theta)) [underbrace((bold(h)^T C_w^(-1) bold(h))^(-1) bold(h)^T C_w^(-1) bold(x), g(bold(x))) - theta]
$

即线性高斯模型的 MVU 存在并且为 $hat(theta) = (bold(h)^T C_w^(-1) bold(h))^(-1) bold(h)^T C_w^(-1) bold(x)$，其方差达到 CRLB。

=== Blackboard bkup

$
x = A + omega, quad omega ~ cal(N)(0, y_2 I)
$

$
ln 1(?): ln 1/(pi A)^(1/2) - (x-A)^2/A
$

$
ln l(?): ln 1/(pi A)^(N/2) - (sum (x[n] - A)^2)/A
$

$
s(x, A)(?): -N/(2A) + [(sum x^2[n])/A^2-N] attach(=, t: ?) T(A) [g(x) - A], quad g(x) ~ hat(A)_"MVUE"(?)
$
