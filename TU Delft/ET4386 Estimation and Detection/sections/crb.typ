#import "../generic.typ": *

#import "@preview/cetz:0.4.2"
#import "@preview/fletcher:0.5.8" as fletcher: diagram, node, edge

= Cramér-Rao Lower Bound

#Cre("TODO")

设估计量 $hat(theta) = g(bold(x))$，其中 $bold(x) = [x[0], x[1], dots, x[N - 1]]^T$。若其为，则其期望应等于真值：

$
forall theta, quad EE(hat(theta)) = integral g(bold(x)) p(bold(x); theta) dif bold(x) = theta
$

== Score Function and Regularity Conditions

定义 score function 为对数似然函数的梯度，用以衡量其陡峭程度：

$
s(bold(x); theta) = (partial ln p(bold(x); theta))/(partial theta)
$

// 其期望：

// $
// EE[s(bold(x); theta)] = EE[(partial ln p(bold(x); theta))/(partial theta)]
// $

若 $s(bold(x); theta)$ 存在且有限，并且有（#Cre("TODO") 何意味）：

$
integral (partial p(bold(x); theta))/(partial theta) dif bold(x) = partial/(partial theta) integral p(bold(x); theta) dif bold(x)
$

则概率分布 $p(bold(x); theta)$ 符合如下正则条件（regularity condition），即 score function 的均值为零：

$
forall theta, quad EE[s(bold(x); theta)] = EE[(partial ln p(bold(x); theta))/(partial theta)] = 0
$

除非该 PDF 的非零定义域取决于 $theta$（#Cre("TODO") 何意味）。

它们的用处是，若满足正则条件，我们就可以估计出估计量方差的下界（lower bounds），从而有希望得到 MVU。

== Fisher Information

对于 score function 的方差，我们可以证明其为 Fisher information，即：

$
I(theta) = -EE[(partial^2 ln p(bold(x); theta))/(partial theta^2)] = -EE[((partial ln p(bold(x); theta))/(partial theta))^2]
$

证明如下：

#Cre("TODO") P14

== CRLB Theorem

== CRLB for the Gaussian Models
