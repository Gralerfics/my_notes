#import "../generic.typ": *

#import "@preview/cetz:0.4.2"
#import "@preview/fletcher:0.5.8" as fletcher: diagram, node, edge

= Fundamentals

== Markov Decision Process (MDP)

马尔可夫决策过程（Markov decision process，MDP）常建模为：$(S, A, P, R, gamma)$，其中 $S$ 为状态空间，$A$ 为动作空间，$P(s' mid(|) s, a)$ 为状态转移概率，$R(s, a, s')$ 为奖励函数，$gamma in [0, 1)$ 为折扣因子。

马尔可夫过程的核心假设是马尔可夫性：

$
P(S_(t+1) mid(|) S_t, A_t, S_(t-1), A_(t-1), ...)
=
P(S_(t+1) mid(|) S_t, A_t)
$

即状态转移只与当前状态和当前动作有关，与更早的历史信息无关。

== Trajectory, Reward and Return

不断采取动作、在状态间转移并获取奖励（reward），走出的一条轨迹实例可以写作：

$
tau = (s_0, a_0, r_1, s_1, a_1, r_2, s_2, dots)
$

具体地，从 $s_0$ 开始，每一步在状态 $s_t$ 时采取动作 $a_t$，导致转移到 $s_(t+1)$ 同时获取该步的奖励 $r_(t+1)$：

$
s_0
stretch(->, size: #200%)^(a_0)_(r_1)
dots
stretch(->, size: #200%)^(a_(t-1))_(r_t)
s_t
stretch(->, size: #200%)^(a_t)_(r_(t+1))
s_(t+1)
stretch(->, size: #200%)^(a_(t+1))_(r_(t+2))
dots
$

从时刻 $t$ 开始往后的总体奖励称为回报（return）：

$
G_t &= r_(t+1) + gamma G_(t+1) \
&= r_(t+1) + gamma r_(t+2) + gamma^2 r_(t+3) + dots \
&= sum_(k=0)^infinity gamma^k r_(t+k+1)
$

其中折扣因子 $gamma$ 用来刻画 “长期奖励贡献递减” 的特性，不考虑太远的奖励，但也不能只看眼前的奖励。

== Policy and Value/Action-Value Functions

策略用 $pi (a mid(|) s)$ 表示，即在状态 $s$ 下采取动作 $a$ 的概率，根据策略采样就能采出一条轨迹 $tau$。
