#import "../generic.typ": *

#import "@preview/cetz:0.4.2"
#import "@preview/fletcher:0.5.8" as fletcher: diagram, node, edge

= Estimation Theory

== Minimum Variance Unbiased Estimator (MVUE)

=== Philosophy

// Supposing to fetch a set of *random samples* $X = {X_1, X_2, dots, X_N}, forall 1<=n<=N$ from the *probabiliy distribution* $p_(X_n) (x_n; theta)$.

// The subscript $X_n$ in $p$ indicates whose distribution it is; $x_n$ is the virtual variable; the variables after the semicolon such as $theta$ represent the *parameters* of the distribution.

从概率分布

$
p_(X_n) (x_n; theta)
$

中采样得一系列样本 $X = {X_1, X_2, dots, X_N}, forall 1<=n<=N$。其中，$p$ 下标的 $X_n$ 代表是谁的概率分布，$x_n$ 是作为自变量的形式变量，分号后的参量如 $theta$ 表示分布的参数。

有时我们可能会简化这个概率分布的写法，比如省略随机变量和参数的标记写为 $p(x_n)$ 等，但只是简化的写法，不代表省略的信息不存在。

// Then what we want is to:
// + Recover the unknown parameters $theta$ from the measurements $X$;
// + Provide a performance measure of the estimated parameters $hat(theta)$;
// + Discuss the statistical *optimality*.

接下来我们希望做的是：
+ 从观测数据 $X$ 中恢复未知参数 $theta$；
+ 提供一个评估所选参数估计量 $hat(theta)$ 的性能指标；
+ 讨论所选估计量的统计*最优性*。

这分别是*解决问题*、*解决得如何*以及*还有没有可能解决得更好*三件事。

=== Estimates and Estimators

对估计问题的解答方式是设计一个估计量（estimator）。这是一种统计量，可以类比均值、方差等统计量；而估计值（estimate）是具体的值，例如对一系列具体的样本计算出的一个估计量的具体值。可以认为估计量是一个函数、一种方式，而估计值是该函数在具体情形下的值。

作为对问题的回答，估计量的设计是任意的，什么都可以作为估计量，但解答自然有效果好坏之分。对于一个估计量，主要可以从总体上*估计得准不准*和*稳不稳定*两个角度来评估其表现。"准不准" 涉及估计量的无偏性（Unbiasedness），详见 @sec:mvue_unbiasedness；"稳不稳定" 涉及估计量的方差，详见 @sec:mvue_mini_var。

// 具体点说，对于一组给定的、具体的观测数据，它可能来自对观测模型的采样，但是在观测的瞬间就已经是确定的数值。于是用这组观测数据配合选取的估计量，计算得到的估计值也将是一个确定的数值。前面所说的 “准不准” 和 “稳定不稳定” 是无法在这样的一次估计中体现的，而是对多次估计结果统计特性的总结。无偏性保证多次观测平均结果逼近真值，而小方差则确保单次估计都更接近真值——形象一些就是那个经典的打靶例子。

=== Unbiasedness <sec:mvue_unbiasedness>

无偏就是估计量的期望和真值一致，该性质确保了在样本量足够大时估计值是准确的：

$
EE(hat(theta)) = theta quad "or" quad "bias"(hat(theta), theta) = EE(hat(theta)) - theta = 0
$

=== Minimum Variance <sec:mvue_mini_var>

估计量的方差可以看作在不同样本组合下计算出的估计值的方差，该方差越小估计量越可靠。具体地，同样样本量下，两个无偏估计量，方差小的那个计算出的估计值靠近真值的概率更大。

这也可以看作是收敛速度的一种体现，可以说方差小的统计量更高效（efficient）。

=== Cost Functions and Optimality Criterion

#Cre("TODO") Cost functions

为了评估估计量的好坏，需要规定好具体的指标。常用的指标如均方误差（mean square error，MSE）定义为：

$
"mse"(hat(theta)) &= EE[(hat(theta) - theta)^2] \
&= EE{[(hat(theta) - EE(hat(theta))) + (EE(hat(theta)) - theta)]^2} \
// &= EE[(hat(theta) - EE(hat(theta)))^2] + (EE(hat(theta)) - theta)^2 \
// &= "var"(hat(theta)) + "bias"(hat(theta))^2 \
&= underbrace(EE[(hat(theta) - EE(hat(theta)))^2], "var"(hat(theta))) + underbrace((EE(hat(theta)) - theta), "bias"(hat(theta)))^2
$

以上拆分表明 MSE 可以视作对偏差和方差的综合评估指标。顺便，*无偏估计量的偏差为零，故其 MSE 就等于其方差*。

对于一个估计量，通过最小化理论 MSE 推导得出的 “最优” 估计量称为 MSE 估计量（MSE estimator）。看起来很完美，但有时这个估计量的表达式会依赖未知参数，导致实际上无法实现（not realizable）。举个例子，对于测量模型：

$
x[n] = theta + w[n], quad n = 0, dots, N - 1
$

欲估计参数 $theta$ 我们采用形式如下的估计量：

$
hat(theta) = a / N sum_(n=0)^(N-1) x[n]
$

其中 $a$ 是待定系数，现在我们希望通过调节 $a$ 得到一个 MSE 最小的估计量。简单计算其期望和方差：

$
EE[hat(theta)] = a theta, quad "var"[hat(theta)] = (a^2 sigma^2) / N
$

计算偏差，平方后同方差相加得 MSE 的表达式：

$
"mse"(hat(theta)) = (a^2 sigma^2) / N + (a - 1)^2 theta^2
$

通过令其对 $a$ 的导数 $(dif "mse"(hat(theta)))/(dif a)$ 为 $0$ 求解最优的系数 $a$，过程略，最终得到：

$
a^* = theta^2 / (theta^2 + sigma^2 \/ N)
$

该结果包含待估计参数的真值 $theta$，故实际情况下无法实现。

=== Minimum Variance Unbiased Estimator (MVUE)

由前，我们有时无法稳定地通过最小化 MSE 的方法获得一个好的估计量，但我们总得规定一个 “最优”。大部分情况下，我们会考虑使用*最小方差无偏*（minimum variance unbiased，MVU）性质作为最优估计量的标准。具体地，我们通常首先希望保证估计是无偏的，从而确保准确性，然后再追求估计的稳定性，即小方差。

对于无偏估计量 $hat(theta)$ 有 $"mse"(hat(theta)) = "var"(hat(theta))$，若任何其他估计量 $tilde(theta)$ 的方差都不低于它，即：

$
forall tilde(theta), quad "var"(hat(theta)) <= "var"(tilde(theta))
$

则 $hat(theta)$ 为所有 $theta$ 上的最小方差无偏估计量。

这样的*全局 MVUE 并不保证存在*，例如很多情况下估计量的方差其实随着真值变化而变化，这导致真值处于不同区间时 MVUE 可能是不同的。

此外，即便 MVUE 存在，也没有一个固定的方法去找到它，通常需要具体情况具体分析，具体将在后面的章节中讨论。
