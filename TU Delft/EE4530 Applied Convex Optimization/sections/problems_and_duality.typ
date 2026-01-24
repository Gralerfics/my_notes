#import "../generic.typ": *

#import "@preview/cetz:0.4.2"
#import "@preview/fletcher:0.5.8" as fletcher: diagram, node, edge

= Canonical Problems

== Optimization Problems

标准优化问题可写成如下形式：

$
&"minimize"& quad &f_0 (x) \
&"subject to"& &f_i (x) <= 0, &quad &i = 1, dots, m \
&& &h_i (x) = 0, &&i = 1, dots, p
$

其中，$x in RR^n$ 是优化变量，$f_0: RR^n -> RR$ 是目标（或代价）函数，$f_i: RR^n -> RR$ 是不等式约束函数，$h_i: RR^n -> RR$ 是等式约束函数。

#underline[最优值（optimal value）]是：

$
p^* = inf {f_0 (x) | f_i (x) < 0, i = 1, dots, m, h_i (x) = 0, i = 1, dots, p}
$

也就是满足约束的最小的#underline[代价函数值]。若其为 $infinity$ 则说明问题#underline[不可行（infeasible）]，即没有 $x$ 满足约束；若其为 $-infinity$ 则说明问题#underline[无下界（unbounded below）]。

#Cre("TODO")
// feasible x；
// optimal x；$X_"opt"$ 是 set of optimal points；
// locally optimal，若有 R > 0 使得 x is optimal for 同样但加了邻域约束 $norm(z-x)_2 <= R$ 的优化问题。

#Cre("TODO")
// ！！！！！注意隐式约束（implicit constraints），主要是指定义域的约束，所有函数的定义域的交集为 $cal(D)$ 称为问题的定义域；
// 写问题里的是显式约束，没有显式约束的问题是无约束问题（unconstrained）。

== Feasibility Problem

#Cre("TODO") 将标准问题中代价函数设为 $f_0 (x) = 0$ 的一类特殊问题，只检查约束，$p^* = 0$ 则说明 feasible，任意 feasible $x$ 都是最优的，若为 $infinity$ 则说明 infeasible。

== Convex Optimization Problems

凸优化问题的标准形式可表述如下：

$
&"minimize"& quad &f_0 (x), &&&&"is convex" \
&"subject to"& &f_i (x) <= 0, &quad &i = 1, dots, m, &quad &f_i "are convex" \
&& &a^T x = b_i, &&i = 1, dots, p, &&"all affine"
$

注意其中要求 $f_0$ 和 $f_i$ 都为凸函数，且等式约束都是仿射的。经常也#underline[把等式约束整体写成 $A x = b$ 的形式]。

若 $f_0$ 是 quasiconvex 的（其他 $f_i$ 还是要求为凸）则称问题为 quasiconvex 的。

#Cre("TODO") 重要性质：凸优化问题的可行解集是凸集。

#Cre("TODO") local and global optima，凸优化问题中任何局部最优点都是全局最优点。证明思路假设有两个局部最优一个比另一个低，连一条线那么这条线上一定有在高的那个邻域内比它更低的，和它局部最优的假设冲突了。

#Cre("TODO") 对于可微的 $f_0$，$x$ 最优的充要条件是它 feasible 且有：

$
nabla f_0 (x)^T (y - x) >= 0, quad forall "feasible" y
$

也就是在 $x$ 处，所有到可行点的方向和梯度方向的内积都是大于等于零的，即在同一侧超平面中，而梯度方向是下降方向的反方向，所以就是选任意可行点都不会再让代价函数下降了。

// #Cre("TODO") （？？？）如果非零则梯度 $nabla f_0 (x)$ 定义了一个在 $x$ 处支持可行解集的的超平面。

// 例子：1、对于无约束问题，$nabla f_0 (x) = 0$ 是 $x$ 最优的充要条件，意思就是没有任何方向可以下降了；2、对于等式约束 Ax = b 问题，最优的充要条件是存在一个 nu 使得 $nabla f_0 (x) + A^T nu = 0$：
// #image("/assets/image-42.png")
// ，总之因为要满足 Ax=b 的约束，方向d必须属于A的零空间，那样A(x+td)=b才继续满足，然后如果最优了就说明梯度的方向要不可行，也就是和这个零空间正交，也就是属于A的行空间，也就是负梯度可以由A^T的列向量拼出来；3、对于非负象限最小化问题，最优的充要条件是x_i>0时梯度i分量=0，而x_i=0时梯度i分量大于或等于0（负梯度的i分量小于或等于0），就是如果不在边界上，梯度必须已经不能动了，如果在边界上则还可以允许它朝负方向，因为约束限制死了也过不去。

=== Equivalent Convex Problems

#Cre("TODO") 当两个问题的解可以简单地互相推得，则这两个问题 (informally) equivalent。一些常见变换可以保持问题的凸性：

#Cre("TODO") 消除等式约束。把 $A x=b$ 解一下表示成 $x = F z+x_0$，然后不用 $x$ 了，全部换成 $F z+x_0$，并改为优化 $z$，就去掉了仿射的等式约束。其中 $x_0$ 是一个特解即有 $A x_0 = b$，而 $F$ 的列空间是 $A$ 的零空间，于是 $F z$ 所得向量属于 $A$ 的零空间，于是 $A(F z) = 0$，所以 $A x = A(F z + x_0) = A F z + A x_0 = b$ 满足约束。

#Cre("TODO") 引入等式约束。和上面正相反，如果 $f_i$ 里面可以提取出仿射关系那么就可以加新变量，新变量和原变量的仿射关系定义就是一种等式约束，然后全一起优化。

#Cre("TODO") 为线性不等式约束引入 slack 变量。$a_i^T x <= b_i$ 可以变为 $a_i^T + s_i = b_i$ 且 $s_i >= 0$，$s_i$ 是新引入的 slack 变量，一起优化。

#Cre("TODO") 上图形（epigraph form），是标准凸优化问题的一种等价形式，就是将问题转为优化 $f_0$ 的上界 $t$，其他约束不变，加一条 $f_0 (x) - t <= 0$。因为 $t$ 最小化了那么 $f_0$ 也就最小化了。

#Cre("TODO") 在部分变量范畴上最小化。就是拆分最小化问题，如果代价函数某个变量没什么约束，那其实可以把代价函数换成一个最小化/下确界函数，在这里面先把一个变量取可能值时的函数值取个最小的，然后优化问题里再去优化另一个。

=== Linear Program (LP)

#Cre("TODO")

$
&"minimize"& quad &c^T x + d \
&"subject to"& &G x prec.curly.eq h \
&& &A x = b
$

#Cre("TODO") 全都 affine，可行集是一个 polyhedron。

// 例子：1、diet problem；2、代价函数是仿射函数的最大值，那么就可以引入一个 t，然后加约束各子函数都 <= t，然后去优化 f_0 = t；3、多面体的切比雪夫中心，……，好例子。

=== Quadratic Program (QP)

#Cre("TODO")

$
&"minimize"& quad &1/2 x^T P x + q^T x + r, quad P in SS_+^n \
&"subject to"& &G x prec.curly.eq h \
&& &A x = b
$

#Cre("TODO") 相比 LP 就是把代价函数换成二次的

// 例子：1、最小二乘；2、带随机代价的线性规划，也就是虽然是线性规划，但变量是随机变量，且代价函数是均值和方差的某种组合，例如 MSE，于是会有二次项。

=== Quadratically Constrained Quadratic Program (QCQP)

#Cre("TODO")

$
&"minimize"& quad &1/2 x^T P_0 x + q_0^T x + r_0, &quad &P_0 in SS_+^n \
&"subject to"& &1/2 x^T P_i x + q_i^T x + r_i <= 0, &&i = 1, dots, m, quad P_i in SS_+^n \
&& &A x = b
$

#Cre("TODO") 如果 $P_1, dots, P_m in SS_(++)^n$，则可行集是 $m$ 个椭球和一个仿射集的交集。

=== Second-Order Cone Programming (SOCP)

#Cre("TODO")

$
&"minimize"& quad &f^T x \
&"subject to"& &norm(A_i x + b_i)_2 <= c_i^T x + d_i, quad i = 1, dots, m \
&& &F x = g
$

其中 $A_i in RR^(n_i times n)$，$F in RR^(p times n)$，不等式约束为二阶锥（SOC）约束：

$
(A_i x + b_i, c_i^T x + d_i) in "second-order cone in" RR^(n_i + 1)
$

#Cre("TODO") 在 $n_i = 0$ 时退化为 LP，在 $c_i = 0$ 时退化为 QCQP，即比 QCQP 和 LP 更一般性。看起来 SOCP 的目标函数是线性的，QCQP 的凸二次代价函数怎么退化过来？可以引入辅助变量 t，添加约束 “那个二次代价函数 <= t”，当 Q 半正定时可以写成 Q^(1\/2) x 的二范数平方 <= t - q^T x，两边取根号就变成了二阶锥约束的一部分，这个 t - q^T x 还可以在引入一个 u >= 0，令其 >= u 来消去线性项？TODO

=== Robust Linear Programming

#Cre("TODO") 如前，LP 可能有时可以要求以一定概率满足约束条件。

// 1、deterministic approach：把参数 a_i 不确定性刻画成要求在一个椭球内的参数取值都需要被满足，然后可以转 SOCP；
// 2、stochastic approach：假设参数 a_i 是高斯的，然后去约束满足约束的概率，然后也转 SOCP。

=== Semidefinite Program (SDP)

#Cre("TODO")

$
&"minimize"& quad &c^T x \
&"subject to"& &x_1 F_1 + x_2 F_2 + dots + x_n F_n + G prec.curly.eq 0, quad F_i, G in SS^k \
&& &A x = b
$

其中的不等式约束称为线性矩阵不等式（linear matrix inequality，LMI），若有多组 LMI 可以合并为一组，用分块对角矩阵拼起来即可。

// LP and SOCP as SDP。

#Cre("TODO") 首先 LP 可以转为 SDP 问题，即将 $A x prec.curly.eq b$ 转为 $bold("diag")(A x - b) prec.curly.eq 0$，也就是逐元素小于或等于零变成对角矩阵半负定。

#Cre("TODO") 然后 SOCP 可以转为 SDP，即将 $norm(A_i x + b_i)_2 <= c_i^T x + d_i, quad i = 1, dots, m$ 转为：

$
mat(
    delim: "[",
    (c_i^T x + d_i) I, A_i x + b_i;
    (A_i x + b_i)^T, c_i^T x + d_i
) succ.curly.eq 0, quad i = 1, dots, m
$

=== Eigenvalue Minimization

$
&"minimize"& quad &lambda_max (A(x))
$

其中 $A(x) = A_0 + x_1 A_1 + dots + x_n A_n$，$A_i in SS^k$。这可以等价转化为如下 SDP：

$
&"minimize"& quad &t \
&"subject to"& &A(x) prec.curly.eq t I
$

#Cre("TODO") 其实就是用了 $lambda (A - t I) = lambda (A) - t$，用特征值的定义可以简单证明。

= Duality

== Lagrange Dual Problem

#underline[Lagrangian]：#Cre("TODO")

$
L(x, lambda, nu) = f_0 (x) + sum_(i=1)^m lambda_i f_i (x) + sum_(i=0)^p nu_i h_i (x)
$

// lambda_i 对应的是 f_i (x) <= 0，注意方向符号
