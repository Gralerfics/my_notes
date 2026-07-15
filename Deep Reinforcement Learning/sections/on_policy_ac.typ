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

== Parameterized Stochastic Policies

我们采用的策略是*随机策略*（stochastic policy），即参数化的是一个概率分布，而不同于 DQN 中直接通过参数化价值函数确定（deterministic）的策略。原因在于，下节中会推导策略梯度定理，其中#underline[涉及（对数）策略分布的梯度]，如果采用确定策略就意味着策略分布退化为 $delta$ 分布，梯度的质量也就会很一般。下个章节的 DDPG、TD3 等算法会基于确定策略展开。

常用的随机策略分布通常也就那几种，离散动作的情况下多采用 *Softmax 策略*：

$
pi_theta (a mid(|) s) = (exp(z_theta (s, a)))/(sum_(a') exp(z_theta (s, a')))
$

非常常见，分布的结构固定后其中的分布参数 $z_theta$ 等可以再用其他模型（如神经网络）进行参数化。Softmax 策略对应的对数梯度为：

$
nabla_theta ln pi_theta (a mid(|) s)
&= nabla_theta [z_theta (s, a) - ln sum_(a') exp(z_theta (s, a'))] \
&= nabla_theta z_theta (s, a) - (sum_(a') exp(z_theta (s, a')) nabla_theta z_theta (s, a'))/(sum_(a') exp(z_theta (s, a'))) \
&= nabla_theta z_theta (s, a) - sum_(a') pi_theta (a' mid(|) s) nabla_theta z_theta (s, a')
$

也就是动作 $a$ 的策略梯度减去所有动作策略梯度按策略分布加权的平均值。

对于连续动作空间则常采用*高斯策略*，设动作维度为 $d$：

$
pi_theta (a mid(|) s) &= cal(N) (mu_theta (s), Sigma_theta (s)) \
&= 1/((2 pi)^(d\/2) abs(Sigma_theta (s))^(1\/2)) exp(-1/2 (a - mu_theta (s))^top Sigma_theta^(-1) (s) (a - mu_theta (s)))
$

对应的对数梯度推导略（#Cre("TODO")）。

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
&= sum_(k=0)^(t-1) gamma^k R_k + gamma^t underbrace(sum_(k=0)^infinity gamma^k R_(k+t), G_t) \
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
$ <equ:EE_At_Nabla_theta_ln_pi_theta_A_S_H_zero>

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

REINFORCE 是典型的 on-policy 策略梯度算法，直接使用当前策略采样并更新，策略用概率分布表示，是随机策略。其#underline[简单直接，不需要估计价值函数]，缺点则是#underline[样本利用效率低（朴素 on-policy 算法的通病，策略更新后旧样本不再能使用），估计方差较大（见后）]。

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
$ <equ:baseline_v_deri>

这里第一步将无限时域求和提到期望外面，严格来说需要满足一些有界性条件，在我们的问题假设下（例如奖励有界、折扣因子等）基本都成立，此处省略。第二步是#underline[条件期望分解]，$EE_(pi_theta)$ 简写的随机变量是整条轨迹，由于期望内部不涉及 $S_t$、$A_t$ 以后的随机变量，可以将它们舍去，剩下的就可以分解为 $H_t$（到 $S_t$ 为止的轨迹）和 $A_t$ 两层。最后为零的项来自 @equ:EE_At_Nabla_theta_ln_pi_theta_A_S_H_zero。

该式为零也就说明在 $G_t$ 上减去任意仅和 $H_t$ 有关的 baseline 都不破坏估计的无偏性。我们通常#underline[选取值函数估计] $V^(pi_theta) (S_t)$ 当 baseline，它仅和当前状态 $S_t$ 有关，符合条件。选择值函数作为 baseline 是符合直觉的，$G_t - V^(pi_theta) (S_t)$ 可以解读为#underline[在该状态下采取某动作所得回报比 "平均情况下在该状态所得回报" 好多少]，相当于做了一个归一化。这个对 "好多少" 的衡量也被称为*优势*（advantage），对应*优势函数*（advantage function）：

$
A^pi (s, a) := Q^pi (s, a) - V^pi (s)
$

#blockquote([
    *关于 baseline 的选取与对方差的优化*：

    此外，选值函数直观简单，但在 "最小化梯度方差" 这一目标下并非最优。#Cre("TODO") 到底怎么降低的方差。
    // https://www.zhihu.com/question/344367451/answer/813387514
])

#underline[接下来的问题是 REINFORCE 框架下并没有估计值函数]，所以需要一个分支负责估计值函数。值函数是对当前策略价值的衡量，评估一个策略的好坏，所以该分支称为*批评者*（critic），维护一个参数化的*值函数估计* $v_(phi.alt) (s) approx V^(pi_theta) (s)$；原本负责探索、采样的分支则称为*执行者*（actor），并用*对优势的估计*进行策略更新。

此即 Actor-Critic 名称的由来，优势估计和值函数估计也是 Actor-Critic 架构下最关键的两个问题。

=== Value Function Estimation

我们先#underline[讨论值函数估计的具体方法]，即 critic 如何学习 $v_(phi.alt) (s)$。基于采样所得轨迹 $tau_infinity = {s_t, a_t, r_t}_(t=0)^infinity$，采用均方误差损失函数：

$
cal(L) [phi.alt] := sum_(t=0)^infinity (overbrace(underbrace(y_t (tau_infinity), "targets") - underbrace(v_(phi.alt) (s_t), "value"), "advantage"))^2
$

其中目标（targets）的选择有很多种。规范起见，定义轨迹样本从 $l$ 时刻到 $r$ 时刻的切片：

$
tau_(l:r) = {s_l, a_l, r_l, s_(l+1), dots, r_(r-1), s_r}
$

顺便定义 $tau_r := tau_(:r) := tau_(0:r)$ 表示从头开始到 $r$ 结束的片段，$tau_infinity$ 也可以纳入其中。

当手上有一段长度为 $n$ 的轨迹切片 $tau_(t:t+n)$ 时，最自然的想法就是把这些样本全都利用上。前 $n$ 步有样本就用样本计算，余下的远未来项用值函数估计，此即#underline[多步目标]（n-step targets）：

$
y_t^n (tau_(t:t+n)) := sum_(k=0)^(n-1) gamma^k r_(t+k) + gamma^n v_(phi.alt) (s_(t+n))
$ <equ:ac_v_estimation_n_step_target>

课件中括号里写的是 $tau_n$ 而非 $tau_(t:t+n)$，只表示了长度容易混淆，这里我们明确写出切片范围。这个 $tau_(l:r)$ 只是用来标注一下定义 $y_t^n$ 所需的变量，实际起点和长度取决于 $y_t^n$ 的上下标，#underline[简单起见]有时也会写包含所有信息的 $tau_infinity$，或者直接省略括号内的参数。

计算可以得知对于任意 $n in NN$，样本来自 $pi$ 的#underline[多步目标的期望]都等于 $V^pi (s_t)$：

$
EE_pi [y_t^n mid(|) s_t]

&= EE_pi [sum_(k=0)^(n-1) gamma^k R_(t+k) + gamma^n v_(phi.alt) (S_(t+n)) mid(|) S_t = s_t] \

&=^! EE_pi [sum_(k=0)^(n-1) gamma^k R_(t+k) + gamma^n V^pi (S_(t+n)) mid(|) S_t = s_t] \

&= EE_pi [sum_(k=0)^(n-1) gamma^k R_(t+k) + gamma^n EE_pi [sum_(k=0)^infinity gamma^k R_(k+t+n) mid(|) S_(t+n), dots, S_t] mid(|) S_t = s_t] \

&= EE_pi [sum_(k=0)^(n-1) gamma^k R_(t+k) + sum_(k=n)^infinity gamma^k R_(t+k) mid(|) S_t = s_t] \

&= EE_pi [sum_(k=0)^infinity gamma^k R_(k+t) mid(|) S_t = s_t] \

&= V^pi (s_t)
$

注意到当 $n=1$ 时，$y_t^1 (tau_infinity) := r_t + gamma v_(phi.alt) (s_(t+1))$ 就是常用的*时序差分*（temporal difference，TD）目标，特点是收敛慢（迭代次数多）而方差较小；当所用样本足够长时就是 *Monte-Carlo* 目标，特点是收敛快（迭代次数少）但方差较大：

$
y_t^"MC" (tau_infinity) := lim_(n->infinity) y_t^n (tau_infinity)
$

现在问题来了，TD 目标和 MC 目标各有优劣，如何综合二者的特点？直接使用多步目标虽然是一种折中方案，但固定了样本长度，并没有实现综合考虑短期和长期目标的效果。考虑无限长样本 $tau_infinity$，我们可以对所有长度 $[1, infinity)$ 的多步目标进行加权平均，选一个系数 $lambda in [0, 1]$，权重就采用几何级数 $lambda^k$，这称为 *$"TD"(lambda)$ 目标*：

$
y_t^("TD"(lambda)) (tau_infinity) := (1-lambda) sum_(k=0)^infinity lambda^k y_t^(k+1) (tau_infinity)
$

其中 $1-lambda$ 是归一化因子，用以确保 $sum_(k=0)^infinity lambda^k = 1/(1-lambda)$ 乘上它后能将权重和化为 $1$。$lambda$ 较小时级数衰减较快，较长的多步目标权重就低，对应地就偏向短期目标。可以验证，$lambda -> 0$ 时它就是 TD 目标 $y_t^1$，$lambda -> 1$ 时它就是 Monte-Carlo 目标 $y_t^"MC"$（#Cre("TODO")）。

=== Advantage Actor-Critic (A2C)

接下来#underline[讨论优势估计的具体方法]。Actor-Critic 只是架构的名称，在我们选取值函数作为 baseline，将其与回报期望组合为 "优势" 进行估计时，实际上已经算是 *Advantage Actor-Critic* 即 A2C 算法。朴素的优势估计使用 $g_t^i - v_(phi.alt) (s_t^i)$，对应样本损失函数：

$
cal(L)'_(pi_theta) [theta] := -1/m sum_(i=1)^m sum_(t=0)^(n-1) gamma^t ln pi_theta (a_t^i mid(|) s_t^i) [#Cbl($g_t^i - v_(phi.alt) (s_t^i)$)]
$ <equ:mc_advantage_loss>

需要注意的是这里的 $g_t^i$ 需要通过第 $i$ 条轨迹样本的连续片段计算得到，故这种优势估计是基于轨迹片段（episodes）的 *Monte-Carlo* 估计。这种估计是#underline[无偏]的，且因为样本信息量大，训练较快；但整条的轨迹包含大量的随机性，导致样本#underline[方差通常很大]，这部分方差无法靠减去一个 baseline 优化。

为此，顺着经典思路，我们考虑与 Monte-Carlo 估计相对的 *TD* 估计。具体地，假设 $v_(phi.alt) (s_(t+1)^i) approx V^(pi_theta) (s_(t+1)^i)$，即允许用当前估计近似未来目标，就可以用 TD 估计作为优势估计，得到相应的样本损失函数：

$
cal(L)''_(pi_theta) [theta] := -1/m sum_(i=1)^m sum_(t=0)^(n-1) gamma^t ln pi_theta (a_t^i mid(|) s_t^i) [#Cbl($r_t^i + gamma v_(phi.alt) (s_(t+1)^i) - v_(phi.alt) (s_t^i)$)]
$ <equ:td_advantage_loss>

由于引入 bootstrapping，这个估计是#underline[有偏]的，但因为单个样本短，能显著#underline[降低方差]。

回顾上节，值函数估计中我们引入 $"TD"(lambda)$ 目标实现 MC 和 TD 的综合，优势估计中自然也有类似的方案，即*广义优势估计*（generalized advantage estimation，GAE）：

$
hat(A)_t^("GAE"(lambda),n) (tau_n) := (1-lambda)/(lambda-lambda^(n-t+1)) sum_(k=1)^(n-t) lambda^k (y_t^k (tau_(t:t+k)) - v_(phi.alt) (s_t))
$ <equ:gae_slides>

课件中用的符号是 $G_t^lambda$，但我们之前用 $G$ 表示折扣回报了，故换成 $hat(A)$，更直白地说明它是对优势函数的估计。这里我们考虑的样本片段 $tau_n$ 到 $s_n$ 就结束的序列，所以求和索引 $k$ 的上限就到 $n-t$，确保 $y_t^k (tau_(t:t+k))$ 有定义；如果是无限样本长度（$n->infinity$）就都是 $infinity$ 了，可以直接省去上标中的 $n$。

从公式角度#underline[直观看一下这个式子]，其实就是用多步目标计算的优势 $y_t^k (tau_(t:t+k)) - v_(phi.alt) (s_t)$ 用几何级数 $lambda^k$ 加权平均，前面是归一化系数。这和 $"TD"(lambda)$ 思路是完全一样的，只是一个在值函数估计里，一个在优势估计里。

#blockquote([
    *关于另一种常见的 GAE 定义及两种定义的关联性*：

    在其他材料中 GAE 经常定义为：

    $
    hat(A)_t^("GAE"(lambda),n) (tau_n) := sum_(l=0)^(n-t-1) (gamma lambda)^l delta_(t+l), quad delta_t := r_t + gamma v_(phi.alt) (s_(t+1)) - v_(phi.alt) (s_t)
    $ <equ:gae_pop>

    这种写法是将 GAE 写成了 TD 优势 $delta_t$ 的加权和形式，形式上简洁一些。我们希望可以证明 @equ:gae_slides 和 @equ:gae_pop 里的两种定义是一致的，但先说结论：在无限长样本下确实如此，但*样本长度有限时二者并不等价*。
    
    先计算优势：

    $
    &y_t^k - v_(phi.alt) (s_t) \
    
    =& sum_(l=0)^(k-1) gamma^l r_(t+l) + gamma^k v_(phi.alt) (s_(t+k)) - v_(phi.alt) (s_t) \
    
    =& sum_(l=0)^(k-1) gamma^l r_(t+l) + gamma^k v_(phi.alt) (s_(t+k)) - v_(phi.alt) (s_t) 
    + sum_(l=1)^(k-1) gamma^l v_(phi.alt) (s_(t+l)) - sum_(l=1)^(k-1) gamma^l v_(phi.alt) (s_(t+l)) \
    
    =& sum_(l=0)^(k-1) gamma^l r_(t+l) + sum_(l=1)^k gamma^l v_(phi.alt) (s_(t+l)) - sum_(l=0)^(k-1) gamma^l v_(phi.alt) (s_(t+l)) \
    
    =& sum_(l=0)^(k-1) gamma^l r_(t+l) + sum_(l=0)^(k-1) gamma^l dot gamma v_(phi.alt) (s_(t+l+1)) - sum_(l=0)^(k-1) gamma^l v_(phi.alt) (s_(t+l)) \
    
    =& sum_(l=0)^(k-1) gamma^l (r_(t+l) + gamma v_(phi.alt) (s_(t+l+1)) - v_(phi.alt) (s_(t+l))) \

    =& sum_(l=0)^(k-1) gamma^l delta_(t+l)
    $

    代入 @equ:gae_slides 有：

    $
    &(1-lambda)/(lambda-lambda^(n-t+1)) sum_(k=1)^(n-t) lambda^k (#Cbl($y_t^k - v_(phi.alt) (s_t)$)) \

    =& (1-lambda)/(lambda-lambda^(n-t+1)) sum_(k=1)^(n-t) lambda^k sum_(l=0)^(k-1) gamma^l delta_(t+l) \

    =& (1-lambda)/(lambda-lambda^(n-t+1)) sum_(l=0)^(n-t-1) sum_(k=l+1)^(n-t) lambda^(k-l) lambda^l gamma^l delta_(t+l) \

    =& sum_(l=0)^(n-t-1) (lambda gamma)^l delta_(t+l) dot (1-lambda)/(lambda-lambda^(n-t+1)) sum_(k=l+1)^(n-t) lambda^(k-l) \

    =& sum_(l=0)^(n-t-1) (lambda gamma)^l delta_(t+l) dot underbrace((1-lambda)/(lambda-lambda^(n-t+1)) sum_(k=1)^(n-t-l) lambda^(k), C(l)) \
    $

    关键步骤就是交换了一下求和顺序，画一下 $k$-$l$ 平面内求和的三角形区域就行。不过#underline[结果很遗憾]，两种定义形式并非完全一致，而且 $C(l)$ 与 $l$ 相关，二者也不是简单的比例关系。但在 $n->infinity$ 时，由于有 $n-t-l = n-t+1 -> infinity$，系数 $C(l)$ 就等于 $1$ 了，即有：

    $
    hat(A)_t^("GAE"(lambda)) (tau_infinity) = sum_(l=0)^infinity (lambda gamma)^l delta_(t+l)
    $ <equ:gae_inf_n_pop>

    @equ:gae_slides 的定义更直观，而 @equ:gae_pop 实际上通常也是以 $n->infinity$ 的形式出现，之后我们主要以前者为主。
])

按照预期，我们可以验证#underline[在两个极端情况下 GAE 分别应该接近 TD 估计和 MC 估计]。当 $lambda -> 0$ 时：

$
lim_(lambda -> 0) hat(A)_t^("GAE"(lambda),n) (tau_n) &= y_t^1 (tau_n) - v_(phi.alt) (s_t) \
&= r_t + gamma v_(phi.alt) (s_(t+1)) - v_(phi.alt) (s_t)
$

这就是 TD 估计，对应 @equ:td_advantage_loss。当 $lambda -> 1$ 时，应用洛必达法则：

$
lim_(lambda -> 1) hat(A)_t^("GAE"(lambda),n) (tau_n)
&= lim_(lambda -> 1) (1-lambda)/(lambda-lambda^(n-t+1)) sum_(k=1)^(n-t) lambda^k (y_t^k (tau_(t:t+k)) - v_(phi.alt) (s_t)) \
&= lim_(lambda -> 1) (dif/(dif lambda) [(1-lambda) sum_(k=1)^(n-t) lambda^(k-1) (y_t^k (tau_(t:t+k)) - v_(phi.alt) (s_t))])/(dif/(dif lambda) [1-lambda^(n-t)]) \
&= 1/(n-t) sum_(k=1)^(n-t) (y_t^k - v_(phi.alt) (s_t))
$

或者直接从定义角度，$lambda->1$ 时权重都相同，所以就是求等权平均，直接写出系数得到结果。注意它*不等于* MC 估计，*但*当 $n->infinity$ 时它就收敛到 MC 估计，用 @equ:gae_inf_n_pop 推导比较方便：

$
hat(A)_t^("GAE"(1)) (tau_infinity) &= sum_(l=0)^infinity gamma^l delta_(t+l) \
&= sum_(l=0)^infinity gamma^l (r_(t+l) + gamma v_(phi.alt) (s_(t+l+1)) - v_(phi.alt) (s_(t+l))) \
&= sum_(l=0)^infinity gamma^l r_(t+l) + sum_(l=0)^infinity gamma^(l+1) v_(phi.alt) (s_(t+l+1)) - sum_(l=0)^infinity gamma^l v_(phi.alt) (s_(t+l)) \
&= sum_(l=0)^infinity gamma^l r_(t+l) + sum_(l=1)^infinity gamma^l v_(phi.alt) (s_(t+l)) - sum_(l=0)^infinity gamma^l v_(phi.alt) (s_(t+l)) \
&= g_t - v_(phi.alt) (s_t) \
$

此外，#underline[GAE 是优势函数的无偏估计]（在 $v_(phi.alt) (s_t)$ 准确的前提下），可以对 GAE 求期望验证：

$
&EE_(pi) [hat(A)_t^("GAE"(lambda),n) (tau_n) mid(|) S_t = s_t, A_t = a_t] \
=& (1-lambda)/(lambda-lambda^(n-t+1)) underbrace(sum_(k=1)^(n-t) lambda^k, (1-lambda^(n-t+1))/(1-lambda) - 1) (underbrace(EE_(pi) [y_t^k (tau_n) mid(|) s_t, a_t], Q^pi (s_t, a_t)", "forall k) - v_(phi.alt) (s_t)) \
=& Q^pi (s_t, a_t) - v_(phi.alt) (s_t) \
=^!& A^pi (s_t, a_t)
$

总结一下，*无限时域* GAE 的两个极端情况分别退化到 TD 和 MC 估计，和上节值函数估计中的无限时域 $"TD"(lambda)$ 相互对应。

=== Asynchronous Advantage Actor-Critic (A3C)

A3C 算法的核心在于异步（asynchronous）。A3C 并行地运行多个环境，每个环境都是一个 A2C，在自己线程里获取和更新全局网络参数，环境之间互不干扰。

异步避免了全局等待，可提高 CPU/GPU 利用率，并且由于多环境并行，按时间顺序自然地会产生多样化的数据，起到了一部分经验重放缓存的作用，即削弱非 i.i.d. 样本之间的相关性。

=== Convergence of Actor-Critic Methods

若满足满足以下假设，Actor-Critic 方法所得回报通常可以收敛到一个（局部）最大值：
+ MDP 是有限的、遍历的并且奖励有界。
+ actor 和 critic 所用参数化模型是线性的，并使用了恰当的基函数。
+ #underline[actor 的学习率低于 critic]。
+ critic 所用 $"TD"(lambda)$ 目标的 $lambda$ 足够大。

对于深度 Actor-Critic 方法没有特别的保证，但：
+ actor 和 critic 的相对学习率还是很重要。
+ GAE 或 $"TD"(lambda)$ 必须 propagate future reward fast，即要快速将未来奖励的影响传递回来，对应的就是使用较大的 $lambda$（即偏向 MC），早点利用未来奖励的样本，而非慢慢等 $v_(phi.alt) (s)$ 等更新。

=== Continuous Control

Actor-Critic 方法可以处理连续动作 $bold(a) in cal(A) subset RR^m$，因为它不需要像 DQN 那样通过离散的最大值函数选取策略；且如果 critic 估计值函数而非动作值函数，不依赖动作，在连续动作场景下就不会导致来不及学习、不稳定的问题。

举个例子，可以让策略网络输出分布参数而非动作，比如用对角高斯分布建模策略：

$
pi_theta (bold(a) mid(|) s) prop exp(-sum_(i=1)^m (a_i - mu_theta (s)_i)^2/(2 sigma_theta (s)_i^2))
$

网络的输出则是 $m$ 个 head 用于 $bold(mu)_theta (s)$ 加上 $m$ 个 head 用于 $bold(sigma)_theta (s)$。

简单情况下，探索采样时可以直接从该策略分布中采样。参数更新时可应用*最大熵正则化*（maximum entropy regularization），定义熵：

$
cal(H) [pi_theta] := -1/n sum_(t=0)^(n-1) integral pi_theta (bold(a) mid(|) s_t) ln pi_theta (bold(a) mid(|) s_t) dif bold(a)
$ <equ:maximum_entropy_reg_entropy>

其表达了策略 $pi_theta$ 在整条轨迹上的平均随机程度。在损失函数中加入正则项：

$
macron(cal(L)) [theta] := cal(L) [theta] - 1/beta cal(H) [pi_theta], quad beta > 0
$

最小化该损失函数时也会一定程度上最大化这个熵，从而使得最优化出来的策略仍然具有随机性。在 on-policy 算法中，采样所用的行为策略就是目标策略。最大熵保持了目标策略的一定随机性，也就鼓励了探索的随机性，避免过早收敛到确定策略，导致错过一些状态空间。顺便给出正则项的梯度：

$
nabla_theta cal(H) [pi_theta] = -1/n sum_(t=0)^(n-1) integral pi_theta (bold(a) mid(|) s_t) ln pi_theta (bold(a) mid(|) s_t) nabla_theta ln pi_theta (bold(a) mid(|) s_t) dif bold(a)
$

=== Conclusion

如 @fig:policy_opt_spectrum 所示，A2C/A3C 依旧属于#underline[随机策略]、#underline[on-policy 损失函数]方法，不同于 REINFORCE 的是又引入了 value-based 的内容（critic 用 target networks 进行值函数估计），往 DQN 方向靠近了一点。

== Trust Region Methods

=== Off-Policy Gradients

前面提到过，Actor-Critic 等 on-policy 算法的典型问题在于样本利用率低，只能使用当前策略采样的样本进行更新。为了更好地说明这一点，我们试着直接推导 off-policy 梯度，即#underline[尝试使用来自另一个策略 $mu$ 的样本来更新当前策略 $pi_theta$]。

首先约定：
+ 行为策略 $mu (a mid(|) s)$ 满足 $(pi_theta (a mid(|) s))/(mu (a mid(|) s)) < infinity, forall a in cal(A), forall s in cal(S)$。
+ $xi_t^(pi_theta) (s)$ 表示采用策略 $pi_theta$ 执行 $t$ 步后状态 $s$ 的分布。

// #Cre("TODO ！！！！！这里的 pi 有问题，有些地方似乎不是用 pi_theta 而是固定为 pi")

由此之前 REINFORCE 的损失函数梯度可以进一步变化：

$
nabla_theta cal(L)_(pi_theta) [theta]

&= -EE_(pi_theta) [sum_(t=0)^infinity gamma^t G_t #Cgr($nabla_theta ln pi_theta (A_t mid(|) S_t)$)] \

&= -EE_(pi_theta) [sum_(t=0)^infinity gamma^t EE_(pi_theta) [G_t mid(|) S_t, A_t] #Cgr($nabla_theta ln pi_theta (A_t mid(|) S_t)$)] quad "(Tower property)" \

&= -integral xi_t^(pi_theta) (s_t) integral cancel(pi_theta (a_t mid(|) s_t)) sum_(t=0)^(n-1) gamma^t #Cbl($EE_(pi_theta) [G_t mid(|) s_t, a_t]$) #Cgr($(nabla_theta pi_theta (a_t mid(|) s_t)) / cancel(pi_theta (a_t mid(|) s_t))$) dif a_t dif s_t \

&= -sum_(t=0)^(n-1) integral xi_t^(pi_theta) (s_t) integral gamma^t #Cbl($Q^(pi_theta) (s_t, a_t)$) nabla_theta pi_theta (a_t mid(|) s_t) dif a_t dif s_t \

&= -sum_(t=0)^(n-1) integral xi_t^(pi_theta) (s_t) integral #Cpu($mu (a_t mid(|) s_t)$) gamma^t Q^(pi_theta) (s_t, a_t) (nabla_theta pi_theta (a_t mid(|) s_t)) / (#Cpu($mu (a_t mid(|) s_t)$)) dif a_t dif s_t
$

这里把期望展开成积分形式和求和形式只是连续和离散的区别，意思差不多。其中 $Q^(pi_theta)$ 来自和之前值函数类似的时移性质（@equ:vfunc_time_shift_prop）：

$
EE_(pi_theta) [G_t mid(|) s_t, a_t]
&= EE_(pi_theta) [sum_(k=0)^infinity gamma^k R_(k+t) mid(|) S_t = s_t, A_t = a_t] \
&= EE_(pi_theta) [sum_(k=0)^infinity gamma^k R_k mid(|) S_0 = s_t, A_0 = a_t] \
&=: Q^(pi_theta) (s_t, a_t)
$

接下来写成采用策略 $mu$ 的轨迹的期望形式（即 $EE_mu [dot] := EE_(Tau~p_Tau^mu (dot)) [dot]$，与 $EE_pi$ 对应）：

$
nabla_theta cal(L)_(pi_theta) [theta]

&= -sum_(t=0)^(n-1) integral xi_t^mu (s_t) (xi_t^(pi_theta) (s_t))/(xi_t^mu (s_t)) integral mu (a_t mid(|) s_t) gamma^t Q^(pi_theta) (s_t, a_t) (nabla_theta pi_theta (a_t mid(|) s_t)) / (mu (a_t mid(|) s_t)) dif a_t dif s_t \

&= -EE_mu [sum_(t=0)^(n-1) (xi_t^(pi_theta) (S_t))/(xi_t^mu (S_t)) gamma^t Q^(pi_theta) (S_t, A_t) (nabla_theta pi_theta (A_t mid(|) S_t)) / (mu (A_t mid(|) S_t))]
$

接下来要作一些*近似*以便推出可用的形式（#Cre("TODO")），考虑只更新策略分布 $pi_theta (dot)$，其他例如 $xi^(pi_theta) (dot)$、$Q^(pi_theta) (dot)$、$V^(pi_theta) (dot)$ 都视作与参数 $theta$ 无关的常数，或者说 `detach` 阻断其对梯度的贡献，记为 $xi^pi (dot)$ 等。顺便记 $eta_t (S_t) := (xi_t^#Cbl($pi$) (S_t))/(xi_t^mu (S_t))$，有：

$
nabla_theta cal(L)_(pi_theta) [theta]

&= -EE_mu [sum_(t=0)^(n-1) eta_t (S_t) gamma^t Q^#Cbl($pi$) (S_t, A_t) (nabla_theta pi_theta (A_t mid(|) S_t)) / (mu (A_t mid(|) S_t))] \

&=^(*) -EE_mu [sum_(t=0)^(n-1) eta_t (S_t) gamma^t (Q^pi (S_t, A_t) - V^pi (S_t)) (nabla_theta pi_theta (A_t mid(|) S_t)) / (mu (A_t mid(|) S_t))] \

&= -nabla_theta EE_mu [sum_(t=0)^(n-1) eta_t (S_t) gamma^t (Q^pi (S_t, A_t) - V^pi (S_t)) (pi_theta (A_t mid(|) S_t)) / (mu (A_t mid(|) S_t))] \

// &= -nabla_theta EE_mu [sum_(t=0)^(n-1) eta_t (S_t) gamma^t underbrace(#Cbl($(R_t + gamma V^pi (S_(t+1)) - V^pi (S_t))$), A^pi (S_t, R_t, S_(t+1))) (pi_theta (A_t mid(|) S_t)) / (mu (A_t mid(|) S_t))]
$

其中 $*$ 式的导出参考 @equ:baseline_v_deri 中减去 baseline 不影响期望的推导，具体地：

$
&EE_mu [sum_(t=0)^(n-1) eta_t (S_t) gamma^t V^pi (S_t) (nabla_theta pi_theta (A_t mid(|) S_t)) / (mu (A_t mid(|) S_t))] \

=& sum_(t=0)^(n-1) gamma^t EE_mu [eta_t (S_t) V^pi (S_t) (nabla_theta pi_theta (A_t mid(|) S_t)) / (mu (A_t mid(|) S_t))] \

=& sum_(t=0)^(n-1) gamma^t integral_(s_t) xi_t^mu (s_t) integral_(a_t) cancel(mu (a_t mid(|) s_t)) eta_t (s_t) V^pi (s_t) (nabla_theta pi_theta (a_t mid(|) s_t)) / cancel(mu (a_t mid(|) s_t)) dif a_t dif s_t \

=& sum_(t=0)^(n-1) gamma^t integral_(s_t) xi_t^mu (s_t) eta_t (s_t) V^pi (s_t) underbrace(integral_(a_t) nabla_theta pi_theta (a_t mid(|) s_t) dif a_t, 0) dif s_t \

=& 0
$

按定义 $Q^pi (S_t, A_t) - V^pi (S_t)$ 即优势函数 $A^pi (S_t, A_t)$，此外 $eta_t (S_t) := (xi_t^(pi_theta) (S_t))/(xi_t^mu (S_t))$ 要么难算要么未知模型无法计算，就假设两种策略相近，近似为 $1$。最终梯度式化为：

$
nabla_theta cal(L)_(pi_theta) [theta]

approx -nabla_theta EE_mu [sum_(t=0)^(n-1) gamma^t (pi_theta (A_t mid(|) S_t)) / (mu (A_t mid(|) S_t)) A^pi (S_t, A_t)]
$ <equ:off_policy_gradient_grad_advantage>

这个式子可以理解为按策略比 $(pi_theta (A_t mid(|) S_t)) / (mu (A_t mid(|) S_t))$ 进行*重要性采样*（importance sampling），从而通过实际来自 $mu$ 的优势估计 $hat(A)_t$ 计算假如它来自 $pi_theta$ 时会有的样子。

#blockquote([
    *关于重要性采样*：

    假设有一个函数 $f(X)$，我们希望求其在 $X ~ pi (dot)$ 分布上的期望，即：

    $
    EE_pi [f(X)] = integral_x pi(x) f(x) dif x
    $

    但假如现在我们无法直接从 $pi (x)$ 上进行采样，只能从另一个分布 $mu (x)$ 上采样（比如说样本已经是现成的来自 $mu$ 的样本，或者 $pi$ 太复杂等情况），那么我们就需要通过来自 $mu$ 的样本间接地获取 $pi$ 下的期望：

    $
    EE_pi [f(X)] = integral_x pi(x) f(x) dif x = integral_x mu(x) pi(x)/mu(x) f(x) dif x = EE_mu [pi(x)/mu(x) f(x)]
    $

    这样一来，用 $mu$ 上采集的样本，乘上一个重要性采样比 $pi(x)/mu(x)$ 再求期望就可以得到 $pi$ 上采集样本时的期望。
])

如果要用这个 $nabla_theta cal(L)_(pi_theta) [theta]$ 尝试实现一个 Actor-Critic 架构的 off-policy 算法，主要要考虑如何去估计这个优势函数。记优势估计为 $hat(A)_t$（加标记避免和动作 $A_t$ 混淆），一种方式是做类似 Q-Learning 的一步 TD 估计，由：

$
Q^pi (S_t, A_t) = EE_mu [R_t + gamma V^pi (S_(t+1)) mid(|) S_t, A_t]
$

将 $V^pi (s)$ 用估计的 $v_(phi.alt) (s)$ 代替，由此得到：

$
nabla_theta cal(L)_mu [theta] &approx -nabla_theta EE_mu [sum_(t=0)^(n-1)  gamma^t (#Cbl($R_t + gamma v_(phi.alt) (S_(t+1)) - v_(phi.alt) (S_t)$)) (pi_theta (A_t mid(|) S_t)) / (mu (A_t mid(|) S_t))] \
&= -nabla_theta EE_mu [sum_(t=0)^(n-1) gamma^t (pi_theta (A_t mid(|) S_t)) / (mu (A_t mid(|) S_t)) #Cbl($hat(A)_t$)]
$

当然，这里的优势估计除了这里的 TD 估计，也可以考虑换成之前讨论的 MC、GAE 等形式。实际实现中如果是基于零碎的状态转移样本进行训练（例如 TD 估计）还经常省略会 $gamma^t$，直接使用简洁的 $-nabla_theta EE [(pi_theta (A_t mid(|) S_t)) / (mu (A_t mid(|) S_t)) hat(A)_t]$，毕竟无所谓连续样本之间的关系了。

在此基础上，我们就可以去#underline[复用旧策略 $mu$ 采样的样本]了，从而#underline[提高样本利用率，加快训练速度]。*具体地*，每次都使用当前策略采样获得一系列样本，然后用这些样本#underline[多次]学习更新得到新策略，如此往复。但若实际尝试，则会发现#underline[样本重复使用次数（repetitions）多了训练就容易不稳定]。

=== Trust-Region Policy Optimization (TRPO)

表现不佳的原因首先是刚刚直接近似为 $1$ 的 $(xi_t^pi (S_t))/(xi_t^mu (S_t))$ 项以及重要性采样比 $(pi_theta (A_t mid(|) S_t)) / (mu (A_t mid(|) S_t))$ 容易数值不稳定。例如，对策略 $pi_theta$ 的更新同时改变了状态分布 $xi_t^pi$，而如果对某个状态 $s_t$，新策略会有可能采样到而旧策略不会采样到，那么比值就会趋于无穷，反之趋于 $0$ 等。

// 所以最好要求：

// $
// xi_t^mu (s) >0, forall s in {s mid(|) xi_t^(pi_theta) (s) > 0} subset.eq cal(S)
// $

// 此外，重要性采样比也可能出现类似的问题，不一定满足之前关于其有界的约定。

TRPO 算法的思想很直接，就是围绕旧行为策略 $mu = pi_(theta')$建立一个*信赖域*（trust region），并要求反复使用这部分旧样本更新所得的新策略 $pi_theta$ 不允许超出这个信赖域，于是问题变为带约束的优化问题：

$
min_theta cal(L)_mu [theta] quad "s.t." quad EE_mu [sum_(t=0)^(n-1) D_"KL" [mu (dot mid(|) S_t) || pi_theta (dot mid(|) S_t)]] <= delta
$

信赖域衡量的是分布之间的 "距离"，所以直接使用了 KL 散度（Kullback-Leibler divergence）的概念，它衡量的是用一个近似概率分布 $q$ 来描述真实分布 $p$ 时所引入的信息损失：

$
D_"KL" (p || q) := integral_(-infinity)^infinity p(x) log (p(x))/(q(x)) dif x
$

这个带约束非线性优化问题看着就很棘手，实际实现中可以采用 *Taylor 估计*（Taylor approximation），也就是在当前参数附近做二阶泰勒展开（Hessian），把问题转为*二次约束二次规划*（quadratic constraints quadratic programming，QCQP）问题。

每次更新，当前（旧）策略为 $mu$，参数位于 $theta_mu$，在此处有梯度 $bold(g) := lr(nabla_theta cal(L)_mu [theta] |)_(theta=theta_mu)$ 和 Hessian 矩阵 $bold(H)$。于是损失函数和约束分别近似为：

$
cal(L)_mu [theta] approx bold(g)^top (theta - theta_mu), quad
D_"KL" [mu || pi_theta] approx 1/2 (theta - theta_mu)^top bold(H) (theta - theta_mu)
$

求解上述优化问题得：

$
theta^* = theta_mu + alpha sqrt((2 delta)/(bold(g)^top bold(H)^(-1) bold(g))) bold(H)^(-1) bold(g)
$

这个解称为*自然策略梯度*（natural policy gradient）。此外，由于我们是近似得到的，所以照这样更新还是有可能违反 KL 约束，所以式中的 $alpha$ 就是 TRPO 在此基础上添加的系数，可以通过*线搜索*（line search）去找到尽可能大的符合约束的 $alpha$，具体地，比如说先试 $alpha = 1$，不对就尝试乘一个因子缩小一点再尝试，如此往复直到符合约束。

=== Proximal Policy Optimization (PPO)

在 TRPO 上再稍加修改就是经常使用的*近端策略优化*（proximal policy optimization，PPO）算法了。

如前所述，TRPO 的约束还是太复杂了，虽然可以二阶近似求解但还是有点麻烦，或许可以更简单一点。考虑 $cal(L)_mu [theta]$ 中的策略比 $(pi_theta (A_t mid(|) S_t)) / (mu (A_t mid(|) S_t))$，我们尝试要求它只能卡（"clip"）到 $1 plus.minus epsilon.alt$ 的范围里，这样也能确保旧策略和新策略不要离太远，但比 KL 散度约束要更简单粗暴。具体地，采用如下和损失函数：

$
cal(L)_mu^"clip" [theta] := -EE_mu [sum_(t=0)^(n-1) gamma^t min{(pi_theta (A_t mid(|) S_t)) / (mu (A_t mid(|) S_t)) hat(A)_t, "clip"((pi_theta (A_t mid(|) S_t)) / (mu (A_t mid(|) S_t)), 1-epsilon.alt, 1+epsilon.alt) hat(A)_t}]
$

观察式子，它其实是在截断和不截断所得两种 loss 中选择了较小的一个。考虑优势估计的符号，分类讨论一下会发现当策略比超出范围时，实际也#underline[只有两种情况会采用截断后的信号]：
+ $hat(A)_t > 0$ 且 $(pi_theta (A_t mid(|) S_t)) / (mu (A_t mid(|) S_t)) > 1+epsilon.alt$ 时。
+ $hat(A)_t < 0$ 且 $(pi_theta (A_t mid(|) S_t)) / (mu (A_t mid(|) S_t)) < 1-epsilon.alt$ 时。

这两种情况分别是 "动作好还想大幅提高该动作概率" 以及 "动作不好还想大幅度降低该动作概率"，换句话说就是#underline[只在 loss 会变好时截断以阻止过大幅度的策略改进]。那么截断能干什么？截断处就是把 loss 削平了，于是梯度也就为零，从而#underline[让这部分优势不参与参数更新]。

与之前会不稳定的直接 off-policy 梯度更新法比较，PPO 随着样本重复使用次数提高显著加速训练，同时稳定性也不错，不过在动作集规模 $abs(cal(A))$ 较大时容易不稳定。

回到谱图 @fig:policy_opt_spectrum，TRPO/PPO 算法依旧是#underline[随机策略]、#underline[采用 Actor-Critic 架构估计值函数]的算法。此外章节初提到过，TRPO/PPO 虽然涉及样本复用，但我们还是将其划分为 on-policy 算法；图中的划分方式倒也可以理解，因为基本的损失函数 $cal(L)_mu [theta]$ 确实是基于 off-policy 思想推导的，所以称为 "off-policy loss"，只不过 TRPO/PPO 在新旧策略分布上做了一些约束以确保没有 "off" 得太远。

#blockquote([
    *关于算法原理与工程实现的距离*：

    TRPO 到 PPO 表现的提升大多不是这些算法层面的裁剪带来的，实际上主要依靠 PPO 实现中的一系列工程技巧。毕竟从原理层面上来说，PPO 只是为了降低复杂度实现的对 TRPO 的近似，按理说效果应该

    具体可以查看一篇 Implementation Matters in Deep Policy Gradients: A Case Study on PPO and TRPO（https://arxiv.org/pdf/2005.12729），其中分析了 PPO 标准实现中应用的一系列技巧并做了详细的消融实验。如果将这些技巧应用到 TRPO 上，在某些任务上其性能也可能超越 PPO。还有 https://iclr-blog-track.github.io/2022/03/25/ppo-implementation-details/。

    不止 PPO，绝大部分算法实践中的优秀效果都离不开大量工程技巧的应用，完全按算法朴素地实现通常效果不会太好。
])

=== Group Relative Policy Optimization (GRPO)

强化学习算法也可以用于大语言模型（large language model，LLM）的微调。对于一个预训练语言模型，可以将其推理视为一种预测下一个单词的强化学习策略：

$
pi_theta (a_t mid(|) s_0, a_0, dots, a_(t-1))
$

再加上信赖域机制限制策略更新幅度以防其遗忘预训练知识，此外还需要训练一个奖励模型提供奖励函数，毕竟基本没法手动为这么大的状态空间设计奖励。如果采用 PPO 算法，critic 需要学习一个值函数估计 $v_(phi.alt) (s)$，但这东西将会需要和 LLM 规模相当的参数量，不现实，并且状态空间过大、维度过高很难保证泛化性能。

考虑 Group Relative Policy Optimization（GRPO）算法，对于一次询问 $s_0$ 运行生成 $G$ 条不同轨迹 ${s_n^i}_(i=1)^G$，得到 $G$ 个最终奖励 $r(s_n^i)$。用这一整组奖励的均值替代 $v_(phi.alt) (s_0)$ 作为对值函数的估计：

$
nu(s_0) := 1/G sum_(i=1)^G r(s_n^i) approx V^(pi_theta) (s_0)
$

优势估计也改用 $nu(s)$ 和奖励来近似：

$
hat(A)_t^i
:= underbrace(r_t^i, approx 0) + underbrace(gamma, 1) underbrace(v_(phi.alt) (s_(t+1)^i), approx r(s_n^i)) - underbrace(v_(phi.alt) (s_t^i), approx nu(s_0))
approx (r(s_n^i) - nu(s_0))/"std"({r(s_n^i)}_(i=1)^G)
$

除以标准差保证尺度归一化。很粗暴的近似，直接表达最终奖励相对平均奖励的优势。

GRPO 的表现取决于与任务的适配性，和 PPO 没有绝对的高下之分，不适合效果就可能很差。此外，通常较小的 $G$ 稳定性表现更优。
