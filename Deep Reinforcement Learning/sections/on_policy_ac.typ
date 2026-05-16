#import "../generic.typ": *

#import "@preview/codly:1.3.0": *

#show: codly-init.with()

= Policy Optimization and On-Policy Actor-Critic

考虑与先前学习值函数的方式不同的另一条思路，直接学习策略分布。定义参数化的策略 $pi_theta (a mid(|) s)$，依旧最优化总体期望回报 $J^(pi_theta) := EE_pi [sum_(t=0)^infinity gamma^t R_t]$。先考虑最朴素的方法，即直接求 $J^(pi_theta)$ 的梯度用于梯度上升。

这一类优化策略的方法可以按特点归纳到一张谱图上，如 @fig:policy_opt_spectrum 所示，本章节先讨论左侧三类方法。顺便，图中 TRPO/PPO 被划分在了 off-policy loss 中，这是不太常见的，虽然 TRPO/PPO 会多次使用部分更新前策略的样本，但由于更新幅度小一般还是算作 on-policy 算法，所以放到该章节中，详见后文。

#figure(
  image("../figures/policy_opt_spectrum.png", width: 90%),
  caption: [策略梯度类算法谱图（来自 TU Delft 课件）]
) <fig:policy_opt_spectrum>

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

记 $H_t = {S_0, A_0, R_0, S_1, dots, R_(t-1), S_t}$，对上式第一项每项有：

$
EE_(pi_theta) [C_t nabla_theta ln pi_theta (A_t mid(|) S_t)]

&=^((1)) EE_(H_t) [EE_(A_t ~ pi_theta (dot mid(|) S_t)) [C_t nabla_theta ln pi_theta (A_t mid(|) S_t) mid(|) H_t]] \

&=^((2)) EE_(H_t) [C_t EE_(A_t ~ pi_theta (dot mid(|) S_t)) [nabla_theta ln pi_theta (A_t mid(|) S_t) mid(|) H_t]]
$ <equ:C_t_out_of_expectation_with_H_t>

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

REINFORCE 算法的核心思路就是直接用采样所得轨迹样本去估计这个期望，从而更新策略。我们要最大化总期望回报，即最小化损失函数 $cal(L)_(pi_theta) [theta] := -J^(pi_theta)$，梯度为：

$
nabla_theta cal(L)_(pi_theta) [theta]
&= -EE_(pi_theta) [sum_(t=0)^infinity gamma^t G_t nabla_theta ln pi_theta (A_t mid(|) S_t)] \
&approx nabla_theta [-1/m sum_(i=1)^m sum_(t=0)^(n-1) gamma^t g_t^i ln pi_theta (a_t^i mid(|) s_t^i)]
$

REINFORCE 是典型的 on-policy 策略梯度算法，直接使用当前策略采样并更新，策略用概率分布表示，是*随机策略*（stochastic policy）。其#underline[简单直接，不需要估计价值函数]，缺点则是#underline[样本利用效率低（朴素 on-policy 算法的通病，策略更新后旧样本不再能使用），估计方差较大（见后）]。

如 @fig:policy_opt_spectrum 所示，之前的 DQN 是与其相对的另一个极端 off-policy value-based 算法，完全依赖动作值函数估计，策略是通过最大化动作值函数得到的*确定策略*（deterministic policy）。由于可以复用，样本利用率高，但是在很难学出一个泛化性优秀的价值函数时就容易不稳定。

== Actor-Critic Algorithms

回到 REINFORCE，我们提到它的估计方差（通常）比较大，主要来源于样本回报中的随机性：

$
g_t^i = sum_(j=0)^infinity gamma^j r_(t+j)^i
$

同样都是无偏估计量，在相同样本数量下方差较小的估计量估计更稳定（每个样本更可信），这是我们希望降低估计量方差的原因。

// 先直观上举个例子，有两个状态 $s_1$ 和 $s_2$，和两种动作 $a_1$ 和 $a_2$，假如在样本中某个时刻

之前我们从 $G_0$ 中删掉了 $C_t$，这是各种可能性中与未来无关的共同部分，通常可以降低方差，现在我们考虑是否还有其他共同部分。如前所述，我们希望衡量的是采取各动作后结果的优劣，从而将其作为权重加到各更新方向上，所以仅体现出相对优劣即可。考虑从 $G_t$ 中减去一个 baseline $B_t$：

$
EE_(pi_theta) [sum_(t=0)^infinity gamma^t (G_t - B_t) nabla_theta ln pi_theta (A_t mid(|) S_t)]
$

当然，首先要保证这么做之后估计还是无偏的，借用 @equ:C_t_out_of_expectation_with_H_t 处的 $H_t$ 定义，我们要求选取的 baseline 只能和 $H_t$（链上 $S_t$ 及之前的随机变量）有关，记为 $B_t := b (H_t)$。和前面证明 $C_t$ 项为零的思路基本一致：

$
&EE_(pi_theta) [sum_(t=0)^infinity gamma^t b (H_t) nabla_theta ln pi_theta (A_t mid(|) S_t)] \

=& sum_(t=0)^infinity gamma^t EE_(pi_theta) [b (H_t) nabla_theta ln pi_theta (A_t mid(|) S_t)] \

=& sum_(t=0)^infinity gamma^t EE_(H_t) [EE_(A_t ~ pi_theta (dot mid(|) S_t))[b (H_t) nabla_theta ln pi_theta (A_t mid(|) S_t) mid(|) H_t]] \

=& sum_(t=0)^infinity gamma^t EE_(H_t) [b (H_t) #Cbl($EE_(A_t ~ pi_theta (dot mid(|) S_t))[nabla_theta ln pi_theta (A_t mid(|) S_t) mid(|) H_t]$)] \

=& 0
$

这里第一步将无限时域求和提到期望外面，严格来说需要满足一些有界性条件，在我们的问题假设下（例如奖励有界、折扣因子等）基本都成立，此处省略。该项为零也就说明在 $G_t$ 上减去任意仅和 $H_t$ 有关的 baseline 都不破坏估计的无偏性。

我们通常#underline[选取值函数估计] $V^(pi_theta) (S_t)$ 当 baseline，它仅和当前状态 $S_t$ 有关，符合条件。选择值函数作为 baseline 是符合直觉的，$G_t - V^(pi_theta) (S_t)$ 可以解读为#underline[在该状态下采取某动作所得回报比 "平均情况下在该状态所得回报" 好多少]，相当于做了一个归一化。这个对 "好多少" 的衡量也被称为*优势*（advantage），对应*优势函数*（advantage function）：

$
A^pi (s, a) := Q^pi (s, a) - V^pi (s)
$

#blockquote([
    *关于 baseline 的选取与对方差的优化*：

    此外，选值函数直观简单，但在 "最小化梯度方差" 这一目标下并非最优。TODO
])

#underline[接下来的问题是 REINFORCE 框架下并没有估计值函数]，所以需要一个分支负责估计值函数。值函数是对当前策略价值的衡量，评估一个策略的好坏，所以该分支称为*批评者*（critic），维护一个#underline[参数化的值函数估计] $v_(phi.alt) (s) approx V^(pi_theta) (s)$；原本负责探索、采样的分支则称为*执行者*（actor），并用估计的优势进行策略更新，此即 Actor-Critic 名称的由来。

=== Value Function Estimation

那么 critic 具体怎么学习 $v_(phi.alt) (s)$？基于采样所得轨迹 $tau_infinity = {s_t, a_t, r_t}_(t=0)^infinity$，采用均方误差损失函数：

$
cal(L) [phi.alt] := sum_(t=0)^infinity (overbrace(underbrace(y_t (tau_infinity), "targets") - underbrace(v_(phi.alt) (s_t), "value"), "advantage"))^2
$

其中目标（targets）的选择有很多种。首先是#underline[多步目标]（n-step targets）：

$
y_t^n (tau_n) := sum_(k=0)^(n-1) gamma^k r_(t+k) + gamma^n v_(phi.alt) (s_(t+n))
$

即前 $n$ 步用样本计算，余下的远未来项用值函数估计，对于任意 $n in NN$，多步目标的期望都等于 $V^pi (s_t)$（#Cre("TODO")：bootstrapping？）。

注意到当 $n=1$ 时，$y_t^1 (tau_infinity) := r_t + gamma v_(phi.alt) (s_(t+1))$ 就是 #underline[TD 目标]，收敛慢（迭代次数多）而方差较小。然后是 #underline[Monte-Carlo 目标]，收敛快（迭代次数少）但方差较大：

$
y_t^"MC" (tau_infinity) := lim_(n->infinity) y_t^n (tau_n)
$

在 Monte-Carlo 和 TD 之间用一个系数 $lambda in [0, 1]$ 作折中，即 $"TD"(lambda)$ 目标：

$
y_t^lambda (tau_infinity) := (1-lambda) sum_(n=0)^infinity lambda^n y_t^(n+1) (tau_infinity)
$

在 $lambda = 0$ 时它就是 TD 目标 $y_t^1 (tau_infinity)$；在 $lambda = 1$ 时它就是 Monte-Carlo 目标 $y_t^"MC" (tau_infinity)$。

=== Advantage Actor-Critic (A2C)

Actor-Critic 只是架构的名称，在我们选取值函数作为 baseline，将其与回报期望组合为 "优势" 进行估计时，实际上已经算是 Advantage Actor-Critic 即 A2C 算法。

朴素的优势估计使用 $g_t^i - v_(phi.alt) (s_t^i)$，对应样本损失函数：

$
cal(L)'_(pi_theta) [theta] := -1/m sum_(i=1)^m sum_(t=0)^(n-1) gamma^t ln pi_theta (a_t^i mid(|) s_t^i) [#Cbl($g_t^i - v_(phi.alt) (s_t^i)$)]
$

需要注意的是这里的 $g_t$ 需要通过整条（后半条）轨迹样本计算得到，这种优势估计是基于完整轨迹的 *Monte-Carlo* 估计。减去 baseline 虽然降低了一部分方差，但整条轨迹包含大量的随机性，样本方差通常很大，这部分方差靠 baseline 无法解决。

为此，我们考虑与 Monte-Carlo 估计相对的，前面 Q-Learning 中提到过的*有限差分*（temporal difference，TD）估计。这是基于 bootstrapping 的一步估计，可能引入偏差，但因为单个样本短，能显著降低方差。具体地，假设 $v_(phi.alt) (s_(t+1)^i) approx V^(pi_theta) (s_(t+1)^i)$，用 TD 估计作为优势估计，采用样本损失函数：

$
cal(L)''_(pi_theta) [theta] := -1/m sum_(i=1)^m sum_(t=0)^(n-1) gamma^t ln pi_theta (a_t^i mid(|) s_t^i) [#Cbl($r_t^i + gamma v_(phi.alt) (s_(t+1)^i) - v_(phi.alt) (s_t^i)$)]
$

由于 bootstrapping，这个估计是有偏的（证明 TODO）。

自然，除了 MC 和 TD 估计，与值函数估计中 $"TD"(lambda)$ 对应地，优势估计中也有*广义优势估计*（generalized advantage estimation，GAE）：

$
hat(A)_t^lambda (tau_m) := (1-lambda)/(lambda-lambda^(m-t+1)) sum_(n=1)^(m-t) lambda^n (y_t^n (tau_m) - v_(phi.alt) (s_t))
$

课件中用的符号是 $G_t^lambda$，我们之前用 $G$ 表示折扣回报了，故换成 $hat(A)$，更直白地说明它是对优势函数的估计。当 $lambda -> 0$ 时，$hat(A)_t^lambda (tau_m) -> y_t^1 (tau_m) - v_(phi.alt) (s_t) = r_t + gamma v_(phi.alt) (s_(t+1)) - v_(phi.alt) (s_t)$ 对应 TD 估计；当 $lambda -> 1$ 时，$hat(A)_t^lambda (tau_m) -> sum_(k=0)^(m-t) gamma^k r_(t+k) - v_(phi.alt) (s_t)$ 对应 MC 估计（#Cre("TODO")：推导/理解（加权））。对 GAE 求期望有：

$
&EE_(pi) [hat(A)_t^lambda (tau_m) mid(|) S_t = s_t, A_t = a_t] \
=& (1-lambda)/(lambda-lambda^(m-t+1)) underbrace(sum_(n=1)^(m-t) lambda^n, (1-lambda^(m-t+1))/(1-lambda) - 1) (underbrace(EE_(pi) [y_t^n (tau_m) mid(|) s_t, a_t], Q^pi (s_t, a_t)", "forall n) - v_(phi.alt) (s_t)) \
=& Q^pi (s_t, a_t) - v_(phi.alt) (s_t) \
=^!& A^pi (s_t, a_t)
$

即 GAE 是优势函数的无偏估计（在 $v_(phi.alt) (s_t)$ 估计无偏的前提下）。

顺便注意，值函数估计和优势估计是分开的两件事，但确实在目标/估计选取上是对应的，性质也相似（收敛速度和方差）。

=== Asynchronous Advantage Actor-Critic (A3C)

A3C 算法的核心在于 asynchronous，即异步。A3C 并行地运行多个环境，每个环境都是一个 A2C 在自己线程里获取和更新全局网络参数，互不干扰。

异步避免了全局等待，可提高 CPU/GPU 利用率，并且由于多个环境并行，自然地产生多样化的数据，起到了一部分经验重放缓存的作用，即削弱非 i.i.d. 样本之间的相关性。

=== Convergence of Actor-Critic Methods

若满足满足以下假设，Actor-Critic 方法所得回报通常可以收敛到一个（局部）最大值：
+ MDP 是有限的、遍历的并且奖励有界。
+ actor 和 critic 所用参数化模型是线性的，并使用了恰当的基函数。
+ #underline[actor 的学习率低于 critic]。
+ critic 所用 $"TD"(lambda)$ 目标的 $lambda$ 足够大。

对于深度 Actor-Critic 方法没有特别的保证，但：
+ actor 和 critic 的相对学习率还是很重要。
+ GAE 或 $"TD"(lambda)$ 必须 propagate future reward fast，即要快速将未来奖励的影响传递回来，对应的就是使用较大的 $lambda$，早点利用未来奖励的样本，而非慢慢等 $v_(phi.alt) (s)$ 等更新。

=== Continuous Control

Actor-Critic 方法可以处理连续动作 $bold(a) in cal(A) subset RR^m$，因为它不需要像 DQN 那样通过离散的最大值函数选取策略；且如果 critic 估计值函数而非动作值函数，不依赖动作，在连续动作场景下就不会导致来不及学习、不稳定的问题。

举个例子，可以让策略网络输出分布参数而非动作，比如用对角高斯分布建模策略：

$
pi_theta (bold(a) mid(|) s) prop exp(-sum_(i=1)^m (a_i - mu_theta (s)_i)^2/(2 sigma_theta (s)_i^2))
$

网络的输出则是 $m$ 个 head 用于 $bold(mu)_theta (s)$ 加上 $m$ 个 head 用于 $bold(sigma)_theta (s)$。

简单情况下，探索采样时可以直接从该策略分布中采样。参数更新时可应用*最大熵正则化*（maximum entropy regularization），#Cre("TODO")，旨在鼓励一定程度的随机性探索，避免过早收敛到确定策略。关于这些内容，之后在 exploration 章节再详细说明。

=== Conclusion

如 @fig:policy_opt_spectrum 所示，A2C/A3C 依旧属于#underline[随机策略]、#underline[on-policy 损失函数]方法，不同于 REINFORCE 的是又引入了 value-based 的内容（critic 用 target networks 进行值函数估计），往 DQN 方向靠近了一点。

== Trust Region Methods

=== Off-Policy Gradients

前面提到过，Actor-Critic 等 on-policy 算法的典型问题在于样本利用率低，只能使用当前策略采样的样本进行更新。为了更好地说明这一点，我们试着直接推导 off-policy 梯度，即尝试使用来自另一个策略 $mu$ 的样本来更新当前策略 $pi_theta$。

首先约定：
+ 行为策略 $mu (a mid(|) s)$ 满足 $(pi_theta (a mid(|) s))/(mu (a mid(|) s)) < infinity, forall a in cal(A), forall s in cal(S)$。
+ $xi_t^(pi_theta) (s)$ 表示采用策略 $pi_theta$ 执行 $t$ 步后状态 $s$ 的分布。

于是之前 REINFORCE 的损失函数梯度可以进一步化为：

$
nabla_theta cal(L)_(pi_theta) [theta]

&= -EE_(pi_theta) [sum_(t=0)^infinity gamma^t G_t #Cgr($nabla_theta ln pi_theta (A_t mid(|) S_t)$)] \

&= -integral xi_t^(pi_theta) (s_t) integral pi_theta (a_t mid(|) s_t) sum_(t=0)^(n-1) gamma^t #Cbl($EE_(pi_theta) [G_t mid(|) s_t, a_t]$) #Cgr($(nabla_theta pi_theta (a_t mid(|) s_t)) / (pi_theta (a_t mid(|) s_t))$) dif a_t dif s_t \

&= -sum_(t=0)^(n-1) integral xi_t^(pi_theta) (s_t) integral cancel(pi_theta (a_t mid(|) s_t)) gamma^t #Cbl($Q^(pi_theta) (s_t, a_t)$) (nabla_theta pi_theta (a_t mid(|) s_t)) / cancel(pi_theta (a_t mid(|) s_t)) dif a_t dif s_t \

&= -sum_(t=0)^(n-1) integral xi_t^(pi_theta) (s_t) integral #Cpu($mu (a_t mid(|) s_t)$) gamma^t Q^(pi_theta) (s_t, a_t) (nabla_theta pi_theta (a_t mid(|) s_t)) / (#Cpu($mu (a_t mid(|) s_t)$)) dif a_t dif s_t \
$

这里把期望展开成积分形式和求和形式只是连续和离散的区别，意思差不多。接下来写成采用策略 $mu$ 的轨迹的期望形式（即 $EE_mu [dot] := EE_(Tau~p_Tau^mu (dot)) [dot]$，与 $EE_pi$ 对应）：

$
nabla_theta cal(L)_(pi_theta) [theta]

&= -sum_(t=0)^(n-1) integral xi_t^mu (s_t) (xi_t^(pi_theta) (s_t))/(xi_t^mu (s_t)) integral mu (a_t mid(|) s_t) gamma^t Q^(pi_theta) (s_t, a_t) (nabla_theta pi_theta (a_t mid(|) s_t)) / (mu (a_t mid(|) s_t)) dif a_t dif s_t \

&= -EE_mu [sum_(t=0)^(n-1) (xi_t^(pi_theta) (S_t))/(xi_t^mu (S_t)) gamma^t Q^(pi_theta) (S_t, A_t) (nabla_theta pi_theta (A_t mid(|) S_t)) / (mu (A_t mid(|) S_t))] \

&= -nabla_theta EE_mu [sum_(t=0)^(n-1) (xi_t^(pi_theta) (S_t))/(xi_t^mu (S_t)) gamma^t #Cbl($Q^(pi_theta) (S_t, A_t)$) (pi_theta (A_t mid(|) S_t)) / (mu (A_t mid(|) S_t))] \

&= -nabla_theta EE_mu [sum_(t=0)^(n-1) (xi_t^(pi_theta) (S_t))/(xi_t^mu (S_t)) gamma^t underbrace(#Cbl($(R_t + gamma V^(pi_theta) (S_(t+1)) - V^(pi_theta) (S_t))$), A^(pi_theta) (S_t, R_t, S_(t+1))) (pi_theta (A_t mid(|) S_t)) / (mu (A_t mid(|) S_t))] \
$

这就是尝试用 $mu$ 的样本去更新 $pi_theta$ 的所谓 off-policy 梯度公式。假如现在就直接用这个式子尝试实现一个 Actor-Critic 架构的 off-policy 算法，$V^(pi_theta) (s)$ 用估计的 $v_(phi.alt) (s)$ 近似；$(xi_t^(pi_theta) (S_t))/(xi_t^mu (S_t))$ 难算或者不知道，如果两种策略相近就可以近似为 $1$。于是损失函数梯度为：

$
nabla_theta cal(L)_mu [theta] &approx -nabla_theta EE_mu [sum_(t=0)^(n-1)  gamma^t (R_t + gamma v_(phi.alt) (S_(t+1)) - v_(phi.alt) (S_t)) (pi_theta (A_t mid(|) S_t)) / (mu (A_t mid(|) S_t))] \
&= -nabla_theta EE_mu [sum_(t=0)^(n-1) gamma^t bold(A)_t (pi_theta (A_t mid(|) S_t)) / (mu (A_t mid(|) S_t))] \
&approx nabla_theta cal(L)_(pi_theta) [theta]
$

这个式子大致可以理解为从优势 $A^(pi_theta)_t$（加上标避免和动作 $A_t$ 混淆）中按系数 $(pi_theta (A_t mid(|) S_t)) / (mu (A_t mid(|) S_t))$ 进行*重要性采样*（importance sampling）。// TODO？实际实现中经常省略 $gamma^t$，直接使用 $-nabla_theta EE_mu [sum_(t=0)^(n-1) (pi_theta (A_t mid(|) S_t)) / (mu (A_t mid(|) S_t)) A^(pi_theta)_t]$，形式简洁。

在此基础上，我们就可以去#underline[复用旧策略 $mu$ 采样的样本]了，从而#underline[提高样本利用率，加快训练速度]。具体地，每次使用当前策略采样获得一系列样本，用这些样本多次学习更新得到新策略，如此往复。但若实际尝试，则会发现#underline[样本重复使用次数（repetitions）多了训练就容易不稳定]。

#blockquote([
    *关于重要性采样*：

    假设有一个函数 $f(X)$，我们希望求其在 $X ~ pi (dot)$ 分布上的期望，即：

    $
    E_pi [f] = integral_x pi(x) f(x) dif x
    $

    但假如现在我们无法直接从 $pi (x)$ 上进行采样，只能从另一个分布 $mu (x)$ 上采样（比如说样本已经是现成的来自 $mu$ 的样本，或者 $pi$ 太复杂等情况），那么我们就需要通过来自 $mu$ 的样本间接地获取 $pi$ 下的期望：

    $
    E_pi [f] = integral_x pi(x) f(x) dif x = integral_x mu(x) pi(x)/mu(x) f(x) dif x = E_mu [pi(x)/mu(x) f(x)]
    $

    这样一来，用 $mu$ 上采集的样本，乘上一个重要性采样比 $pi(x)/mu(x)$ 再求期望就可以得到 $pi$ 上采集样本时的期望。
])

=== Trust-Region Policy Optimization (TRPO)

TODOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
// 不稳定的原因归根结底还是在于刚刚直接直接近似为 $1$ 的 $(xi_t^(pi_theta) (S_t))/(xi_t^mu (S_t))$ 项，实际上这个比值很不稳定。如若对某个状态 $s_t$，新策略会有可能采样到而旧策略不会采样到，那么这个比值就会趋于无穷，反之趋于 $0$，即最好要求：

$
xi_t^mu (s) >0, forall s in {s mid(|) xi_t^(pi_theta) (s) > 0} subset.eq cal(S)
$

若新策略较旧策略更新太多，上述情况出现概率就更大，且其他状态下的比值可能也会离 $1$ 比较远，从而导致训练不稳定。

TRPO 算法的思想很粗暴，就是围绕旧行为策略 $mu = pi_(theta')$建立一个*信任域*（trust region），并要求反复使用这部分旧样本更新所得的新策略 $pi_theta$ 不允许超出这个信任域，于是问题变为带约束的优化问题：

$
min_theta cal(L)_mu [theta] quad "s.t." quad EE_mu [sum_(t=0)^(n-1) D_"KL" [mu (dot mid(|) S_t) || pi_theta (dot mid(|) S_t)]] <= delta
$

信任域衡量的是分布之间的 "距离"，所以直接使用了 KL 散度（Kullback-Leibler divergence）的概念，它衡量的是用一个近似概率分布 $q$ 来描述真实分布 $p$ 时所引入的信息损失：

$
D_"KL" (p || q) := integral_(-infinity)^infinity p(x) log (p(x))/(q(x)) dif x
$

这个带约束非线性优化问题看着就很棘手，实际实现中可以采用 *Taylor 估计*（Taylor approximation），也就是在当前参数附近做二阶泰勒展开（Hessian），把问题转为*二次约束二次规划*（quadratic constraints quadratic programming，QCQP）问题。

每次更新，当前（旧）策略为 $mu$，参数位于 $theta_mu$，在此处有梯度 $bold(g) := lr(nabla_theta cal(L)_mu [theta] |)_(theta=theta_mu)$ 和 Hessian 矩阵 $bold(H)$。于是损失函数和约束分别近似为：

$
cal(L)_mu [theta] approx bold(g)^top (theta - theta_mu), quad
D_"KL" [mu || pi_theta] approx 1/2 (theta - theta_mu)^top bold(H) (theta - theta_mu)
$

求解上述优化问题得解：

$
theta^* = theta_mu + alpha sqrt((2 delta)/(bold(g)^top bold(H)^(-1) bold(g))) bold(H)^(-1) bold(g)
$

这个解称为*自然策略梯度*（natural policy gradient）。此外，由于我们是近似得到的，所以照这样更新还是有可能违反 KL 约束，所以式中的 $alpha$ 就是 TRPO 在此基础上添加的系数，可以通过*线搜索*（line search）去找到尽可能大的符合约束的 $alpha$，具体地，比如说先试 $alpha = 1$，不对就尝试乘一个因子缩小一点再尝试，如此往复直到符合约束。

=== Proximal Policy Optimization (PPO)

在 TRPO 上再进一步就是喜闻乐见，在实践中经常使用的*近端策略优化*（proximal policy optimization，PPO）算法了。

如前所述，TRPO 的约束还是太复杂了，虽然可以二阶近似求解但还是有点麻烦，或许可以更简单一点。

// 回到 $cal(L)_mu [theta]$ 中的策略比 $(pi_theta (A_t mid(|) S_t)) / (mu (A_t mid(|) S_t))$，考虑直接将它卡（"clip"）到 $1 plus.minus epsilon.alt$ 的范围里，同样也能确保旧策略和新策略不要离太远，但比 KL 散度约束要更简单粗暴。

具体地，

#Cre("TODO")

与之前会不稳定的直接 off-policy 梯度更新法比较，PPO 随着样本重复使用次数提高显著加速训练，同时稳定性也不错，不过在动作集规模 $abs(cal(A))$ 较大时容易不稳定。

回到谱图 @fig:policy_opt_spectrum，TRPO/PPO 算法依旧是#underline[随机策略]、#underline[采用 Actor-Critic 架构估计值函数]的算法。此外章节初提到过，TRPO/PPO 虽然涉及样本复用，但我们还是将其划分为 on-policy 算法；图中的划分方式倒也可以理解，因为基本的损失函数 $cal(L)_mu [theta]$ 确实是基于 off-policy 思想推导的，所以称为 "off-policy loss"，只不过 TRPO/PPO 在新旧策略分布上做了一些约束以确保没有 "off" 得太远。

=== Group Relative Policy Optimization (GRPO)

#Cre("TODO")
