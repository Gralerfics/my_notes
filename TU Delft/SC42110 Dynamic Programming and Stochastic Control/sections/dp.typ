#import "../generic.typ": *

#import "@preview/cetz:0.4.2"
#import "@preview/fletcher:0.5.8" as fletcher: diagram, node, edge

= Dynamic Programming

== Markov Decision Process (MDP)

考虑系统：

$
X_(t+1) = f(X_t, U_t, W_(t+1))
$

其中，状态变量 $X_t in XX$，$U_t in UU$ 是输入（动作、决策等），$W_t$ 是随机干扰（式中用下标 $t+1$ 为了表明干扰是 $t$ 时刻输入施加后才显现的）。我们称其为*受控随机过程*（controlled stochastic process）。

给定有限时域 $T in NN$，定义：
+ 阶段代价（running/stage cost）函数 $g: XX times UU -> RR$ 表示在时刻 $t in {0, 1, dots, T-1}$、状态为 $x$ 时采取动作 $u$ 的代价。
+ 终端代价（terminal cost）函数 $G: XX -> RR$ 表示在终端时刻 $t = T$ 处于状态 $x$ 时的代价。

假设初始时系统状态 $X_0 = x in XX$，对于给定状态序列 $(X_t)_(t=0)^T = (X_0 = x, X_1, dots, X_T)$ 和输入序列 $(U_t)_(t=0)^(T-1) = (U_0, U_1, dots, U_(t-1))$，其对应的轨迹总代价为：

$
sum_(t=0)^(T-1) g(X_t, U_t) + G(X_T)
$

由于我们考虑的是随机过程，所以常研究轨迹的期望代价：

$
EE[sum_(t=0)^(T-1) g(X_t, U_t) + G(X_T) mid(|) X_0 = x]
$

我们的目标则是求解能最小化这个代价的控制输入序列：

$
min_((U_t)_(t=0)^(T-1) in UU^T) quad &EE[sum_(t=0)^(T-1) g(X_t, U_t) + G(X_T) mid(|) X_0 = x] \
"s.t." quad &X_(t+1) = f(X_t, U_t, W_(t+1))
$

考虑到过程中巨大的随机性，我们基于状态和输入的历史信息设计（时变）闭环控制器 $mu_t$：

$
U_t = mu_t (X_t, U_(t-1), X_(t-1), dots, U_0, X_0) in UU, quad t in {0, 1, dots, T-1}
$

但这么做 $mu_t$ 的状态空间就太大了，要设计函数在任意时刻 $t$ 将每条可能轨迹样本 $(x_t, u_(t-1), x_(t-1), dots, u_0, x_0) in XX times (UU times XX)^t$ 映射到控制信号 $u_t$，计算量不现实。甚至由于控制策略可以是随机策略，即 $mu_t$ 可能用概率分布表示，要求更复杂，不过这门课不考虑。

所以我们只研究一类受控随机过程，其不同时刻的干扰 $W_(t+1)$ 之间在给定 $X_t$ 和 $U_t$ 下条件独立，这样的过程称为*Markov 决策过程*（Markov decision process，MDP），而对 MDP 有如下定理成立：

#blockquote([
    *定理 3.1*：

    对于一个 MDP，存在一个无记忆（memory-less）的最优策略（或称为 Markov 策略）：

    $
    U_t = mu_t (X_t), quad forall t in {0, 1, dots, T-1}
    $
])

于是就不用管那么大的输入空间了，只要基于当前状态设计输入即可实现最优策略。应用 Markov 策略后，MDP 就退化为一条 Markov 链，因为每步的决策被控制器决定好了，之前优化输入序列的优化问题也就可以写作优化控制函数的问题，记代价函数：

$
J_0^mu (x) := EE[sum_(t=0)^(T-1) g(X_t, mu_t (X_t)) + G(X_T) mid(|) X_0 = x]
$

优化问题为：

$
min_((mu_t : XX -> UU)_(t=0)^(T-1)) quad &J_0^mu (x) \
"s.t." quad &X_(t+1) = f(X_t, mu_t (X_t), W_(t+1))
$

#Cre("TODO")

类似 Markov 链，我们也可以用转移概率来描述 MDP 的动力学。例如对于有限的状态空间 $XX = [n] := {1, dots, n}$ 和输入空间 $UU = [m]$，定义动作 $k in UU$ 的转移概率矩阵：

$
P_k (i, j) = PP(X_(t+1) = j mid(|) X_t = i, U_t = k), quad forall (i, j, k) in XX times XX times UU
$

于是转移概率矩阵簇 ${P_k}_(k in UU)$ 就完整描述了 MDP 的动力学。
