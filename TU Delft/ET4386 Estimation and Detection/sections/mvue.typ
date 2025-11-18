#import "../generic.typ": *

#import "@preview/cetz:0.4.2"
#import "@preview/fletcher:0.5.8" as fletcher: diagram, node, edge

= Minimum Variance Unbiased Estimator (MVUE)

== Philosophy

Supposing to fetch a set of *random samples* $X = {X_1, X_2, dots, X_N}, forall 1<=n<=N$ from the *probabiliy distribution* $p_(X_n) (x_n; theta)$.

The subscript $X_n$ in $p$ indicates whose distribution it is; $x_n$ is the virtual variable; the variables after the semicolon such as $theta$ represent the *parameters* of the distribution.

// 假设从概率分布 $p_(X_n) (x_n; theta)$ 中采样得一系列样本 $X = {X_1, X_2, dots, X_N}, forall 1<=n<=N$。其中，$p$ 下标的 $X_n$ 代表是谁的概率分布，$x_n$ 是作为自变量的形式变量，分号后的参量如 $theta$ 表示分布的参数。

Then what we want is to:
+ Recover the unknown parameters $theta$ from the measurements $X$;
+ Provide a performance measure of the estimated parameters $hat(theta)$;
+ Discuss the statistical *optimality*.

这分别是*解决问题*、*解决得如何*以及*能不能解决得更好*三件事。

== Estimates and Estimators

对估计问题的解答方式是设计一个估计量（estimator）。这是一种统计量，可以类比均值、方差等统计量；而估计值（estimate）是具体的值，例如对一系列具体的样本计算出的一个估计量的具体值。可以认为估计量是一个函数、一种方式，而估计值是该函数在具体情形下的值。

作为对问题的回答，估计量的设计是任意的，什么都可以作为估计量，但解答自然有效果好坏之分。对于一个估计量，主要可以从总体上*估计得准不准*和*稳不稳定*两个角度来评估其表现。"准不准" 涉及估计量的无偏性（Unbiasedness），详见 @sec:mvue_unbiasedness；"稳不稳定" 涉及估计量的方差，详见 @sec:mvue_mini_var。

== Unbiasedness <sec:mvue_unbiasedness>

无偏就是估计量的期望和真值一致，该性质确保了在样本量足够大时估计是保证准确的：

$
EE(hat(theta)) = theta quad "or" quad "bias"(hat(theta), theta) = EE(hat(theta)) - theta = 0
$

== Minimum Variance <sec:mvue_mini_var>

估计量的方差可以看作在不同样本组合下计算出的估计值的方差，该方差越小估计量越可靠。具体地，同样样本量下，两个无偏估计量，方差小的那个计算出的估计值靠近真值的概率更大。

这也可以看作是收敛速度的一种体现。

== Cost Functions and Optimality Criterion

#Cre("TODO") Cost functions.

均方误差（mean square error，MSE）定义为：

$
"mse"(hat(theta)) &= EE[(hat(theta) - theta)^2] \
&= EE{[(hat(theta) - EE(hat(theta))) + (EE(hat(theta)) - theta)]^2} \
// &= EE[(hat(theta) - EE(hat(theta)))^2] + (EE(hat(theta)) - theta)^2 \
// &= "var"(hat(theta)) + "bias"(hat(theta))^2 \
&= underbrace(EE[(hat(theta) - EE(hat(theta)))^2], "var"(hat(theta))) + underbrace((EE(hat(theta)) - theta), "bias"(hat(theta)))^2
$

以上拆分表明 MSE 可以视作对偏差和方差的综合评估指标。顺便，无偏估计量的偏差为零，故其 MSE 就等于其方差。

#Cre("TODO") Lec1 P11 Example 2 例子，说明 MSE 估计量不是一直能找到（可能依赖未知量）

== MVUE

对于无偏估计量 $hat(theta)$ 即 $"mse"(hat(theta)) = "var"(hat(theta))$，若任何其他估计量 $tilde(theta)$ 的方差都不低于它，即：

$
forall tilde(theta), quad "var"(hat(theta)) <= "var"(tilde(theta))
$

则 $hat(theta)$ 为所有 $theta$ 上的最小方差无偏估计量（MVU）。

MVU 不保证存在，例如很多情况下不同估计量的方差随着真值变化而变化，故在真值的不同区间可能存在不同的 MVU。
