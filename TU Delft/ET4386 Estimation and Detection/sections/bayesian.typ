#import "../generic.typ": *

#import "@preview/cetz:0.4.2"
#import "@preview/fletcher:0.5.8" as fletcher: diagram, node, edge

== Bayesian Estimators

// Philosophy

#Cre("TODO")

=== Bayesian Mean Square Error (Bmse)

在贝叶斯框架下，待估计的 $theta$ 也被视作随机变量，问题从估计其具体值转为估计其统计特征，这使得我们可以系统性地融合关于 $theta$ 的先验知识，相关详细说明可见 Notes of Filtering & Identification Section 3.5。

我们依旧尝试最小化 MSE，在贝叶斯框架下我们称其为 Bayes MSE（Bmse）。由于 $theta$ 也成为了随机变量，Bmse 相较 MSE 对 $bold(x)$ 的积分还多了一重对 $theta$ 的积分。结果上来看，MSE 还依赖 $theta$ 作为参数（即 $theta$ 不同会导致结果不同），而 Bmse 直接把 $theta$ 的统计特征，即所有可能性都算进去了。

$
"mse"(hat(theta); theta) &= EE[(hat(theta) - theta)^2] = integral (hat(theta) - theta)^2 p(bold(x); theta) dif bold(x) \
"Bmse"(hat(theta)) &= EE[(hat(theta) - theta)^2] = integral.double (hat(theta) - theta)^2 p(bold(x), theta) dif bold(x) dif theta
$

=== Minimum Mean Square Error Estimator (MMSE)

由贝叶斯定理，有 $p(bold(x), theta) = p(theta|bold(x)) p(bold(x))$，代入 Bmse 有：

$
"Bmse"(hat(theta)) = integral.double (hat(theta) - theta)^2 p(bold(x), theta) dif bold(x) dif theta = integral [integral (hat(theta) - theta)^2 p(theta|bold(x)) dif theta] p(bold(x)) dif bold(x)
$

如此我们将关于已知数据 $bold(x)$ 的部分移到外侧。我们希望最小化 Bmse，又由于对所有 $bold(x)$ 都有 $p(bold(x)) >= 0$，我们只需要对所有 $bold(x)$ 最小化中间那部分即：

$
min_(hat(theta)) integral (hat(theta) - theta)^2 p(theta|bold(x)) dif theta
$

求导并令其为零：

$
partial/(partial hat(theta)) integral (hat(theta) - theta)^2 p(theta|bold(x)) dif theta &= 2 integral (hat(theta) - theta) p(theta|bold(x)) dif theta \
&= 2 hat(theta) - 2 integral theta p(theta|bold(x)) dif theta = 0
$

即有：

$
hat(theta) = integral theta p(theta|bold(x)) dif theta = EE(theta|bold(x))
$

即后验分布的均值（mean of the posterior）。

#Cre("TODO") 总结（P21）里有个 commutes over affine transformations。

#Cre("TODO") P22 起的 Vector process，没怎么懂。

==== Example: Gaussian Measurements and Gaussian Prior

#Cre("TODO") P11。

// x[n] = A + w[n]，噪声 w[n] 是高斯，参数 A 的先验是高斯
// 计算 p(bold(x)|A) 和 p(A)，而 p(A|bold(x)) 也是高斯，算出来（？）
// hat(A) 就是 EE(A|bold(x))
// 最后搞出来其实就是样本均值（MVUE）和先验 mu_A 的线性插值，系数是按方差来的
// 若数据量够大总会慢慢趋于样本均值（压过先验）；若无先验，即先验方差趋于无穷则也趋于样本均值
// 对比 Bmse 和 MSE，在贝叶斯框架下前者是更小了，也就是在该框架下可以认为先验知识让估计更准确了

==== Example: Bivariate and Multivariate Gaussian Process

#Cre("TODO") P15。

// x 和 y 联合高斯，协方差矩阵 C
// 反正算出 EE(y|x) 和 "var"(y|x) 之类

$
EE(y|x) &= EE(y) + ("cov"(y, x))/("var"(x)) (x - EE(x)) \
"var"(y|x) &= "var"(y) - ("cov"(x, y)^2)/("var"(x)) = "var"(y) (1 - ("cov"(x, y)^2)/("var"(x)"var"(y))) \
&= "var"(y) (1 - rho^2)
$

#Cre("TODO") P16。

// x 和 y 换成 k x 1 和 l x 1 的向量 bold(x) 和 bold(y)，联合高斯

$
EE(mat(delim: "[", bold(x); bold(y))) = mat(delim: "[", EE(bold(x)); EE(bold(y))), quad C = mat(delim: "[", C_(x x), C_(x y); C_(y x), C_(y y))
$

$
EE(bold(y)|bold(x)) &= EE(bold(y)) + C_(y x) C_(x x)^(-1) (bold(x) - EE(bold(x))) \
C_(y|x) &= C_(y y) - C_(y x) C_(x x)^(-1) C_(x y)
$

#Cre("TODO") P17。

// 例子，bold(x) = bold(1) A + bold(w)，套公式：bold(x) 是 bold(x)，k = N；A 是 bold(y)，l = 1，二者联合高斯
// 套完代矩阵逆引理化简一下

==== Example: General Linear Gaussian Model

#Cre("TODO") P20。

// bold(x) = H bold(theta) + bold(w), quad bold(w) ~ cal(N)(bold(0), C)
// 还有先验 bold(theta) ~ cal(N)(bold(mu)_theta, C_theta)
// 老一套总之把 EE(bold(theta)|bold(x)) 和 C_(bold(theta)|bold(x)) 算出来
// 还是老一套矩阵逆引理化简，得到和 F&I 中类似形式的结果（？？？怎么好像不大一样？），见 Notes of Filtering & Identification Section 3.5。

=== Bayes Risk

#Cre("TODO") P10，怎么感觉和后面 detection 的 Bayes risk 不在说同一个东西。

// 转而最小化一个关于误差 $epsilon = theta - hat(theta)$ 的代价函数 $cal(C)(epsilon)$ 的期望
// 当 $cal(C)(epsilon) = epsilon^2$ 时，最优估计量就是前面 MMSE 那套的 $EE[theta|bold(x)]$，即后验分布的均值
// "Absolute" error：考虑 $cal(C)(epsilon) = abs(epsilon)$，最后会推出估计量是后验分布的中位数
// "Hit-or-miss" error：考虑 $cal(C)(epsilon) = cases(0\, quad &abs(epsilon) <= delta, 1\, &abs(epsilon) > delta)" ", quad "with" delta -> 0$，最后会算出来估计量在后验最大值的附近，i.e., the mode of the posterior

// 这些或都可称为 MMSE？

=== Maximum a Posteriori (MAP) Estimator

#Cre("TODO") P13。

// 相较于 MMSE 中求后验分布的均值，MAP 估计量直接最大化后验分布，计算更简单。
// 和 MLE 的关系
// 和 MMSE 的关系
// one-to-one function 什么的？

=== Linear MMSE Estimator (LMMSE)

#Cre("TODO") P15。

// 除了高斯假设下还行，贝叶斯估计量基本都比较难算
// LMMSE 限定估计量为线性估计量

$
hat(theta) = sum_(n=0)^(N-1) a_n x[n] + a_N
$

// 在这个限定条件下再去最小化 Bmse
// 求解，第一步，hat(theta) 的形式代入 Bmse，对 a_N 求导置零求出 a_N
// 然后代回 Bmse，再对 bold(a) = [a_0, a_1, dots, a_(N-1)]^T 求导置零，求出最优解

$
hat(theta) = EE(theta) + C_(theta x) C_(x x)^(-1) (bold(x) - EE(bold(x)))
$

// 例子：P18；还有 P19 的高斯马尔可夫（这个情况下和 MMSE 一致）
// 顺便，和数据不相关的参数无法被 LMMSE 估计出来
