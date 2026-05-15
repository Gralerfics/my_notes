#import "../generic.typ": *

#import "@preview/codly:1.3.0": *

#show: codly-init.with()

= On-Policy Actor-Critic

考虑与先前学习值函数的方式不同的另一条思路，直接学习策略分布。定义参数化的策略 $pi_theta (A mid(|) S)$，依旧最优化总体期望回报 $J^(pi_theta) := EE_pi [sum_(t=0)^infinity gamma^t R_t]$。先考虑最朴素的方法，即直接求 $J^(pi_theta)$ 的梯度用于梯度上升。

== Policy Gradient Theorem

对于随机轨迹 $Tau = {S_0, A_0, R_0, S_1, dots, R_(n-1), S_n}$，我们之前定义过整条轨迹的折扣回报 @equ:discounted_return，这里我们再定义从某个指定时刻 $t$ 作为起点的折扣回报：

$
G_t := sum_(k=0)^infinity gamma^k R_(t+k)
$

课件上这里用的符号是 $R_t$，但我们这里大写 $R_t$ 已经用于表示奖励的随机变量，为避免混淆采用 $G_t := G_t (Tau)$，样本记为 $g_t := g_t (tau)$。于是总期望回报又可写成 $J^(pi_theta) := EE_(pi_theta) [G_0]$。在继续推导前，我们先将轨迹分布显式写出为：

$
p_Tau^(pi_theta) (Tau) = rho(S_0) product_(t=0)^(n-1) pi_theta (A_t mid(|) S_t) P(S_(t+1) mid(|) S_t, A_t)
$

其对数梯度为：

$
nabla_theta ln p_Tau^(pi_theta) (Tau)
&= nabla_theta [ln rho(S_0) + sum_(t=0)^(n-1) [ln pi_theta (A_t mid(|) S_t) + ln P(S_(t+1) mid(|) S_t, A_t)]] \
&= nabla_theta sum_(t=0)^(n-1) ln pi_theta (A_t mid(|) S_t)
$

接下来推导总期望回报的梯度：

$
nabla_theta J^(pi_theta)

&= nabla_theta EE_(pi_theta) [G_0] \

&= nabla_theta sum_(tau) g_0 (tau) p_Tau^(pi_theta) (tau) \

&= sum_(tau) g_0 (tau) [nabla_theta p_Tau^(pi_theta) (tau)] \

&= sum_(tau) g_0 (tau) [p_Tau^(pi_theta) (tau) nabla_theta ln p_Tau^(pi_theta) (tau)] \

&= sum_(tau) [g_0 (tau) nabla_theta ln p_Tau^(pi_theta) (tau)] p_Tau^(pi_theta) (tau) \

&= EE_(pi_theta) [G_0 nabla_theta ln p_Tau^(pi_theta) (Tau)] \

&= EE_(pi_theta) [G_0 nabla_theta sum_(t=0)^(n-1) ln pi_theta (A_t mid(|) S_t)] \

&= EE_(pi_theta) [sum_(t=0)^(n-1) #Cbl($G_0$) nabla_theta ln pi_theta (A_t mid(|) S_t)]
$

到这一步实际上已经可以用了（用轨迹样本估计这个期望），但我们还可以把估计的方差降低一些。直觉上来说，由于之后的动作和之前的奖励无关，所以 $G_0$ 中 $t$ 时刻之前的部分是可以挪出去的。我们还可以证明挪出去这部分期望为零，可以直接消去，具体地，将 $G_0$ 拆分为：

$
G_0 := sum_(k=0)^infinity gamma^k R_k
&= sum_(k=0)^(t-1) gamma^k R_k + sum_(k=t)^infinity gamma^k R_k \
&= sum_(k=0)^(t-1) gamma^k R_k + gamma^t sum_(k=0)^infinity gamma^k R_(k+t) \
&=: C_t + gamma^t G_t
$

于是之前的结果可拆分整理为：

$
nabla_theta J^(pi_theta)
&= sum_(t=0)^(n-1) EE_(pi_theta) [C_t nabla_theta ln pi_theta (A_t mid(|) S_t)] + EE_(pi_theta) [sum_(t=0)^(n-1) gamma^t G_t nabla_theta ln pi_theta (A_t mid(|) S_t)]
$

记 $H_t = {S_0, A_0, R_0, S_1, dots, S_t}$，对上式第一项每项有：

$
EE_(pi_theta) [C_t nabla_theta ln pi_theta (A_t mid(|) S_t)]

&=^((1)) EE_(H_t) [EE_(A_t ~ pi_theta (dot mid(|) S_t)) [C_t nabla_theta ln pi_theta (A_t mid(|) S_t) mid(|) H_t]] \

&=^((2)) EE_(H_t) [C_t EE_(A_t ~ pi_theta (dot mid(|) S_t)) [nabla_theta ln pi_theta (A_t mid(|) S_t) mid(|) H_t]]
$

上式中第 $(1)$ 步是全概率公式将轨迹分布分解为前半段和 $H_t$ 和 $A_t$（$t$ 时刻之后的随机变量在期望中没有用到所以都可以消去）；第 $(2)$ 步是由于给定 $H_t$ 时 $C_t$ 可视为常数，可提取到期望之外。而其中：

$
EE_(A_t ~ pi_theta (dot mid(|) S_t)) [nabla_theta ln pi_theta (A_t mid(|) S_t) mid(|) H_t]
&= sum_(a_t) nabla_theta ln pi_theta (a_t mid(|) S_t), quad S_t "given by" H_t \
&= sum_(a_t) pi_theta (a_t mid(|) S_t) nabla_theta ln pi_theta (a_t mid(|) S_t) \
&= nabla_theta sum_(a_t) pi_theta (a_t mid(|) S_t) \
&= 0
$

于是总期望回报式中关于 $C_t$ 的项为 $0$，最终有：

$
nabla_theta J^(pi_theta) = EE_(pi_theta) [sum_(t=0)^infinity gamma^t G_t nabla_theta ln pi_theta (A_t mid(|) S_t)]
$

即*策略梯度定理*（Policy Gradient Theorem）。这个式子可以直观上理解一下，$nabla_theta ln pi_theta (A_t mid(|) S_t)$ 是 "可以提高 $S_t$ 下选动作 $A_t$ 概率的参数更新方向"，$gamma^t G_t$ 是采取该动作后所得相对初始时刻的折扣回报，所以总期望回报的参数更新方向就是把所有时刻这些子更新方向按回报加权求和并取期望。#underline[再具体一点，某个动作带来的折扣回报越大，提高该动作概率的参数优化方向就会在总优化方向中占比越大，优化后选取该动作的概率也就越大。]

=== REINFORCE

REINFORCE 算法的核心思路就是直接用采样所得轨迹样本去估计这个期望，从而更新策略。我们要最大化总期望回报，即最小化损失函数 $cal(L)_pi [theta] := -J^(pi_theta)$，梯度为：

$
nabla_theta cal(L)_pi [theta]
&= -EE_(pi_theta) [sum_(t=0)^infinity gamma^t G_t nabla_theta ln pi_theta (A_t mid(|) S_t)] \
&approx nabla_theta [-1/m sum_(i=1)^m sum_(t=0)^(n-1) gamma^t g_t^i ln pi_theta (a_t^i mid(|) s_t^i)]
$

REINFORCE 是典型的 on-policy 策略梯度算法，直接使用当前策略采样并更新，策略用概率分布表示，是*随机策略*（stochastic policy）。其#underline[简单直接，不需要估计价值函数]，缺点则是#underline[样本利用效率低（朴素 on-policy 算法的通病，策略更新后旧样本不再能使用），估计方差较大（见后）]。

之前的 DQN 是与其相对的另一个极端 off-policy value-based 算法，完全依赖动作值函数估计，策略是通过最大化动作值函数得到的*确定策略*（deterministic policy）。由于可以复用，样本利用率高，但是在很难学出一个泛化性优秀的价值函数时就容易不稳定。

== Actor-Critic Algorithms

回到 REINFORCE，我们提到它的估计方差（通常）比较大，主要来源于 $g_t^i = sum_(j=0)^infinity gamma^j r_(t+j)^i$。

// 先直观上举个例子，有两个状态 $s_1$ 和 $s_2$，和两种动作 $a_1$ 和 $a_2$，假如在样本中某个时刻

=== Asynchronous Advantage Actor-Critic (A3C)

=== Continuous Control

== Trust Region Methods

=== Trust-Region Policy Optimization (TRPO)

=== Proximal Policy Optimization (PPO)

=== Group Relative Policy Optimization (GRPO)
