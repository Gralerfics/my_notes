#import "../generic.typ": *

#import "@preview/cetz:0.4.2"
#import "@preview/fletcher:0.5.8" as fletcher: diagram, node, edge

= Least-Squares Estimates

// 统计学估计参数theta，频率学派和贝叶斯学派，最小二乘的例子
// 非线性的情况，非线性最小二乘求解算法，nls.m（？）
// P24，TODO，RLS、SLS、……

== Notations and Symbols

按课件的偏好，我们给出统一的估计问题的模型。首先是*信号模型*：

$
y = F theta + L epsilon, quad epsilon ~ cal(N)(0, I_m)
$

其中 $theta$ 是待估计的参数，$L epsilon$ 是噪声项。规范一下维度，$F in RR^(m times n)" "(m>=n)$，$y in RR^m$， $theta in RR^n$，$epsilon in RR^m$，而 $L in RR^(m times m)$ 是非奇异方阵。*更清晰一些*：

$
[y]_(m times 1) = [F]_(m times n) [theta]_(n times 1) + [L]_(m times m) [epsilon]_(m times 1), quad epsilon ~ cal(N)([0]_(m times 1), I_m)
$

我们令 $epsilon$ 是服从*单位高斯分布*的（均值为 $0$，协方差为 $I$），再用一个分开的 $L$ 来允许其拓展到噪声协方差不为单位矩阵的其他情况。我们也可以直接计算 $L epsilon$ 的均值：

$
E[L epsilon] = L dot E[epsilon] = L dot 0 = 0
$

和协方差：

$
E[(L epsilon - E[L epsilon])(L epsilon - E[L epsilon])^T] &= E[(L epsilon)(L epsilon)^T] \
&=E[L epsilon epsilon^T L^T] \
&=L E[epsilon epsilon^T] L^T \
&=L L^T
$

顺便，根据高斯分布的性质，$L epsilon$ 仍然服从高斯分布（$L$ 是线性变换）。于是*前述模型等效于*把 $L epsilon$ 合并成一个均值为零、协方差为 $L L^T$ 的新 $epsilon$：

$
y = F theta + epsilon, quad epsilon ~ cal(N)(0,L L^T)
$ <equ:lse_note_equv_model>

== The Linear Least-Squares Problem

我们先只考虑单位高斯分布的 $epsilon$，即基于上述信号模型中 $L = I_m$ 的情况（可省去），提出最基础的线性*最小二乘问题*：

$
min_theta epsilon^T epsilon = min_theta norm(epsilon)^2_2, quad s.t." "y = F theta + epsilon, quad epsilon ~ cal(N)(0, I_m)
$

顺便，"s.t." 指 "subject to"。

// #blockquote([
//     #Cre("TODO")1. 这里我们考虑 deterministic 的 $theta$，...
//     #Cre("TODO")2. $epsilon$ 含义的变化，...
// ])

对该问题的最小方差无偏估计量为：

$
hat(theta) = (F^T F)^(-1) F^T y
$

#Cre("TODO")proof of minimum variance unbiased estimator

== The Weighted Linear Least-Squares Problem

当 $L$ 不为 $I_m$ 时，问题为：

$
min_theta epsilon^T epsilon = min_theta norm(epsilon)^2_2, quad s.t." "y = F theta + L epsilon, quad epsilon ~ cal(N)(0, I_m)
$ <equ:lse_wlls_basic>

将条件两侧同乘 $L^(-1)$ 可得到：

$
min_theta epsilon^T epsilon = min_theta norm(epsilon)^2_2, quad s.t." "L^(-1) y = L^(-1) F theta + epsilon, quad epsilon ~ cal(N)(0, I_m)
$

形式同上节，将解中的 $y$ 和 $F$ 换为 $L^(-1) y$ 和 $L^(-1) F$ 得到该问题的解：

$
hat(theta) &= (F^T (L^T)^(-1) L^(-1) F)^(-1) F^T (L^T)^(-1) L^(-1) y \
&= (F^T (L L^T)^(-1) F)^(-1) F^T (L L^T)^(-1) y \
&:= (F^T W F)^(-1) F^T W y
$

其中，记 $W = (L L^T)^(-1)$。

#blockquote([
    由前 @equ:lse_note_equv_model 的思路将 $L epsilon$ 合起来。令 $epsilon' = L epsilon$，即有 $epsilon = L^(-1) epsilon'$，则 @equ:lse_wlls_basic 可化为：

    $
    min_theta (epsilon')^T W epsilon', quad s.t." "y = F theta + epsilon', quad epsilon‘ ~ cal(N)(0, L L^T)
    $

    不想看撇号，把变量换回来，问题等价于：

    $
    min_theta norm(epsilon)_W^2 := min_theta epsilon^T W epsilon, quad s.t." "y = F theta + epsilon, quad epsilon ~ cal(N)(0, L L^T)
    $

    与这里的写法类似，我们之后可能在范数右下角标记 $W$ 来表示这种形式的二次型。

    从这个角度来看，这里的 $W$ 表达了某种 "权重"（#underline("W")eight）的意味，其决定了 $epsilon$ 不同元素（或交叉项）在优化问题中的相对影响力大小，故该问题类型称为加权线性最小二乘问题。
])

== Nonlinear Least-Squares Problem

顺便考虑一下非线性的情况：

$
min_theta norm(epsilon)_2^2, quad y = f(theta) + L epsilon, quad epsilon ~ cal(N)(0, I_m)
$

方法自然就是在 $hat(theta)$ 处将 $f(theta)$ 线性化。由泰勒展开，在 $theta = hat(theta)$ 附近，$f(theta)$ 可一阶近似为：

$
f(theta) approx f(hat(theta)) + lr((dif f(theta))/(dif theta)|)_(theta = hat(theta)) + dots
$

代入 $y = f(theta) + L epsilon$ 即有：

$
underbrace(y - f(hat(theta)), e) approx underbrace(lr((dif f(theta))/(dif theta)|)_(theta = hat(theta)), F) underbrace((theta - hat(theta)), Delta theta) + L epsilon
$

求解算法自然又是迭代法：

#image("../figures/image.png")

#Cre("TODO")文字；nls.m 例程

== The Stochastic Linear Least-Squares Problem


