#import "../generic.typ": *

#import "@preview/cetz:0.4.2"
#import "@preview/fletcher:0.5.8" as fletcher: diagram, node, edge

= Cramér-Rao Lower Bound

#Cre("TODO")

设估计量 $hat(theta) = g(bold(x))$，其中 $bold(x) = [x[0], x[1], dots, x[N - 1]]^T$。若其为无偏估计量，则其期望应等于真值：

$
forall theta, quad EE(hat(theta)) = integral g(bold(x)) p(bold(x); theta) dif bold(x) = theta
$

== CRLB Theorem

对于满足一些条件的一类参数（待估计量）$theta$，接下来介绍的定理可以用于寻找 $theta$ 所有无偏估计量 $hat(theta)$ 的方差所能达到的最小值，即下界。

这也意味着，如果我们可以找到一个无偏估计量，且其方差等于这个下界，那么它就将是我们喜欢的最小方差无偏估计量（MVUE）。

要引出定理，首先来准备一些东西。

=== Score Function and Regularity Conditions

定义 score function 为对数似然函数的梯度，可以衡量其陡峭程度：

$
s(bold(x); theta) = (partial ln p(bold(x); theta))/(partial theta)
$

若 $s(bold(x); theta)$ 存在且有限，并且有（#Cre("TODO") 何意味）：

$
integral (partial p(bold(x); theta))/(partial theta) dif bold(x) = partial/(partial theta) integral p(bold(x); theta) dif bold(x)
$

则概率分布 $p(bold(x); theta)$ 符合如下正则条件（regularity condition）：

$
forall theta, quad EE[s(bold(x); theta)] = EE[(partial ln p(bold(x); theta))/(partial theta)] = 0
$ <equ:crb_regularity_cond>

除非该 PDF 的非零定义域取决于 $theta$（#Cre("TODO") 何意味）。

== Fisher Information

对于 score function 的方差，我们可以证明其为 Fisher information $I(theta)$，定义为：

$
I(theta) = -EE[(partial^2 ln p(bold(x); theta))/(partial theta^2)] = -EE[((partial ln p(bold(x); theta))/(partial theta))^2]
$

证明如下：

#Cre("TODO") P14

== CRLB Theorem

如果参数（待估计量）$theta$ 的分布 $p(bold(x); theta)$ 满足 @equ:crb_regularity_cond 的正则条件，则该参数的任何无偏估计量 $hat(theta)$ 的方差存在下界：

$
"var"(hat(theta)) >= 1 / I(theta) = 1 / (-EE[(partial^2 ln p(bold(x); theta))/(partial theta^2)]) = 1 / (-EE[((partial ln p(bold(x); theta))/(partial theta))^2])
$

== Find the MVU Estimator



== CRLB for the Gaussian Models
