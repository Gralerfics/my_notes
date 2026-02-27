#import "../generic.typ": *

#import "@preview/cetz:0.4.2"
#import "@preview/fletcher:0.5.8" as fletcher: diagram, node, edge

= Stability

== Stability Theory

// Lec 1.2

#Cre("TODO") P5，$abs(x)_cal(A) := inf_(y in cal(A)) norm(x - y)$，即 $x$ 到集合 $cal(A)$ 中的点距离的下界（最小距离）。显然，若 $x$ 在集合中则该值为 $0$。

#Cre("TODO") ……

#Cre("TODO") 连续函数 $sigma: RR_(>=0) -> RR_(>=0)$。

#Cre("TODO") $cal(K)$ 函数：零点为零，严格单调增。

#Cre("TODO") $cal(K)_infinity$ 函数：$cal(K)$ 的基础上，发散到无穷远 $lim_(s->infinity) sigma(s) = infinity$。

#Cre("TODO") 函数 $beta: RR_(>=0) times NN -> RR_(>=0)$，加了一个整数变量。

#Cre("TODO") $cal(K L)$ 函数：1、$forall t >= 0$，$beta(dot, t)$ 是 $cal(K)$ 函数；2、$forall s >= 0$，$beta(s, dot)$ 非增，且 $lim_(t->infinity) beta(s, t) = 0$，即若将 $s$ 作为横轴，调节 $t$ 在坐标轴上绘制一系列 $beta(dot, t)$ 的图像，则随着 $t$ 增加，图像是逐渐向横轴靠拢（下压）的，且无穷远处会压缩到零。

#Cre("TODO") 连续函数 $gamma: RR^n -> RR_(>=0)$，向量映射到非负实数。

#Cre("TODO") 正定（$cal(P D)$）：$gamma(x) = 0 <=> x = 0$，即过原点，按函数值域的定义还隐含函数值非负。

#Cre("TODO") ……

#Cre("TODO") P13，Lyapunov 函数（LF）：$V: RR_n -> RR_(>=0)$，对任意 $x in RR^n$ 满足：
$
alpha_1 (abs(x)_cal(A)) <= V(x) <= alpha_2 (abs(x)_cal(A)) \
V(f(x)) <= V(x) - alpha_3 (abs(x)_cal(A))
$
其中 $alpha_1, alpha_2 in cal(K)_infinity$，$alpha_3 in cal(P D)$。

// 所以参考这里的用法，前面的 K、KL 之类似乎是类似复杂度分析中 Theta、Omega 一类的比较函数；K 函数零点处为 0，单调增，括号内写 $abs(x)_cal(A)$ 则代表其和 $abs(x)_cal(A)$ 细节形状可能有差异但大体性质一致……

#Cre("TODO") ……

== MPC Stability Analysis

=== Without State Constraints

#Cre("TODO") 首先问题是要证明无状态约束、有限时域下，MPC 闭环系统的稳定性。

#Cre("TODO") 闭环系统：
$
x^+ = f(x, kappa_N (x))
$
$x^+ = f(x, u)$ 是离散时间动力学模型，$u$ 的地方替换成根据 MPC 控制策略产生的输入 $kappa_N (x)$ 得到的就是 MPC 闭环系统。

#Cre("TODO") 最优控制序列：
$
u_N^0 (x) = [u_N^0 (0; x), u_N^0 (1; x), dots, u_N^0 (N-1; x)] in UU^N
$
上标 $0$ 代表 “最优”；该序列是按照 MPC 策略，horizon 为 $N$，在处于状态 $x$ 时，优化出的控制序列，称为最有控制序列；其中每一项表示为 $u_N^0 (k; x)$ 代表序列中的第 $k$ 步，所处 $x$ 是产生该序列的前提条件，写在分号后面类似 “参数”。

#Cre("TODO") 所以按定义，控制律 $kappa_N (x)$ 就是挑出以 $x$ 为起点的最优控制序列的第一步：
$
kappa_N (x) = u_N^0 (0; x)
$

#Cre("TODO") 要证明稳定性，考虑 Lyapunov 直接法，即构造一个势能函数 $V(x)$，使其正定并沿任意状态路径单调下降，具体见前。那么选一个什么来当 Lyapunov 函数呢？很自然地我们可以选有限时域最优值函数 $V_N^0 (x)$，也就是我们 MPC 里面优化的那个衡量未来 $N$ 步代价的目标函数，因为 MPC 本身就在一步一步优化它，所以问题变成要证明：
$
V_N^0 (x^+) <= V_N^0 (x) - alpha_3 (abs(x))
$
注意，这个代价函数是我们设计的，我们的目的也是通过设计合适的 stage cost $l$、终端代价函数 $V_f$ 和终端集 $XX_f$ 来使 MPC 稳定。所以我们的目的是构造出一个符合上述条件的代价函数，使其能作为有效的 Lyapunov 函数。好了，现在观察这个条件，看着 $V_N^0 (x^+)$ 和 $V_N^0 (x)$ 很难想出来怎么去设计（why？），所以我们可以强化这个条件，把左边的 $V_N^0 (x^+)$ 替换为一个上界，如果我们构造出了满足新条件的函数，它也就满足原先的 Lyapunov 条件。这其实是一个两头权衡的问题，构造难就放缩条件，但放缩过头符合条件的函数就越少，构造就可能变难。

#Cre("TODO") 上界怎么来？借助最优性的定义就有：
$
V_N^0 (x^+) := V_N (x^+, bold(u)_N^0 (x^+)) = min_u V_N (x^+, bold(u)) <= V_N (x^+, tilde(bold(u))), quad forall tilde(bold(u)) in UU^N
$
也就是 $x^+$ 开始的最优值函数定义为所有可能输入序列导致的值函数之间最小的那个，所以任意构造一个可行的 $tilde(bold(u))$ 它对应的值函数 $V_N (x^+, tilde(bold(u)))$ 都不小于 $V_N^0 (x^+)$，于是新条件：
$
V_N (x^+, tilde(bold(u))) <= V_N^0 (x) - alpha_3 (abs(x))
$

#Cre("TODO") 那么我们构造一个什么 $tilde(bold(u))$？考虑一种移位控制序列，令：
$
tilde(bold(u)) = [u_N^0 (1; x), u_N^0 (2; x), dots, u_N^0 (N-1; x), bb(u)] in UU^N
$
即最优控制序列从第二项开始的序列，末尾加一步补全 $N$ 步。

#Cre("TODO") 把目标函数定义代入条件化简：


