#import "../generic.typ": *

#import "@preview/codly:1.3.0": *

#show: codly-init.with()

= Off-Policy Actor-Critic

#Cre("TODO")

== Deterministic Policy Gradients (DPG)

之前的策略都是采用随机概率分布建模，动作来自对随即策略的采样，现在考虑*确定策略*（deterministic policy），用一个确定的函数映射 $bold(pi)_theta : cal(S) -> cal(A)$ 建模。考虑连续动作 $bold(a) in cal(A) subset RR^m$，critic 用参数化模型 $Q_(phi.alt) (s_t, bold(a)_t)$ 估计值函数 $Q^pi (s_t, bold(a)_t)$，基于之前 off-policy 策略梯度推导的某个中间形式进行近似有：

$
nabla_theta cal(L)_(pi_theta) [theta]

&approx -nabla_theta sum_(t=0)^(n-1) integral_(s_t) xi_t^mu (s_t) integral_(bold(a)_t) gamma^t underbrace(Q^pi (s_t, bold(a)_t), approx Q_(phi.alt) (s_t, bold(a)_t)) pi_theta (bold(a)_t mid(|) s_t) dif bold(a)_t dif s_t
$

而确定策略就相当于某个动作概率为 $1$ 的随机策略（即有 $pi_theta (bold(pi)_theta (s) mid(|) s) = 1, forall s$），故可简化积分得：

$
nabla_theta cal(L)_(pi_theta) [theta]

&approx -nabla_theta sum_(t=0)^(n-1) integral_(s_t) xi_t^mu (s_t) gamma^t Q_(phi.alt) (s_t, bold(pi)_theta (s_t)) dif s_t \

&= -nabla_theta sum_(t=0)^(n-1) integral_(s_t) xi_t^mu (s_t) gamma^t Q_(phi.alt) (s_t, bold(pi)_theta (s_t)) dif s_t \

&= -nabla_theta sum_(t=0)^(n-1) gamma^t EE_mu [Q_(phi.alt) (S_t, bold(pi)_theta (S_t))] \

&prop -nabla_theta EE_(t ~ p(dot mid(|) gamma) prop gamma^t) [EE_mu [Q_(phi.alt) (S_t, bold(pi)_theta (S_t))]]
$

最后一步是将 $t$ 视为服从分布 $p(t mid(|) gamma) prop gamma^t$ 的一个离散随机变量，从而将式子化为期望形式。原本需要在完整轨迹样本上进行计算，如此变化后就可以在数据集上对 $t$ 采样来进行计算，特别是针对以状态转移形式存储的数据（transitions），可以在每条数据中存储其在所在轨迹中的时间步 $t$ 来实现计算。不过实践中有时由于 $gamma$ 接近 $1$、数据不包含时间戳 $t$ 等原因会省略 $gamma^t$。

=== Deep DPG (DDPG)

类似 DQN 之于 Q-Learning，使用深度神经网络作为 DPG 的参数化模型，并加入一系列如经验重放缓存、目标网络和软目标更新、Batch-norm over inputs（BN）、加入噪声随机探索等提高稳定性和表现的技术，得到 DDPG。

要实现 DPG 需要维护两个网络（不算目标网络），一个是作为确定策略的 actor 网络 $bold(pi)_theta (s)$，一个是作为动作值函数估计的 critic 网络 $Q_(phi.alt) (s, bold(a))$。二者分别优化各自的损失函数，actor 优化：

$
arg min_theta -EE_mu [sum_(t=0)^(n-1) gamma^t Q_(phi.alt) (S_t, bold(pi)_theta (S_t))]
$

同时 critic 优化：

$
arg min_(phi.alt) EE_mu [sum_(t=0)^(n-1) (R_t + gamma Q_(phi.alt') (S_(t+1), bold(pi)_(theta') (S_(t+1))) - Q_(phi.alt) (S_t, S_t))^2]
$

用 PyTorch 实现大致如下：

#codly(zebra-fill: none)
```python
from torch.optim import Adam
from torch.nn.functional import mse_loss

actor_optimizer = Adam(actor.parameters())
critic_optimizer = Adam(critic.parameters())
actor_t, critic_t = deepcopy(actor), deepcopy(critic)

for _ in range(max_updates):
    batch = self.replay_buffer.sample()

    # Critic update
    q_values = critic(batch['states'], batch['actions'])
    next_actions = actor_t(batch['next_states'])
    targets = batch['rewards'] + gamma * (˜batch['terminals'] \
              * critic_t(batch['next_states'], next_actions))
    critic_optimizer.zero_grad()
    mse_loss(q_values, targets.detach()).backward()
    critic_optimizer.step()

    # Actor update
    q_values = critic(batch['states'], actor(batch['states']))
    actor_optimizer.zero_grad()
    (-q_values * batch['discounts']).mean().backward()
    actor_optimizer.step()

    # Target network updates
    actor_t = target_model_updates(actor_t, actor)
    critic_t = target_model_updates(critic_t, critic)
```

其中 `actor_t` 和 `critic_t` 是目标网络（target network），类似之前在 DQN 的示例中的使用，`.detach()` 阻断目标网络的梯度更新，`target_model_updates` 中可以具体实现目标网络的同步方式（如软更新）。这里的数据中用 `batch['discounts']` 存储了各条 transition 对应的 $gamma^t$。

两个网络需要分别使用两个优化器进行优化，因为优化策略的 actor loss 中包含值函数网络 $bold(pi)_theta$，值函数网络的参数不应该受到策略更新的影响，我们不希望在 actor 优化时同时影响两个网络的参数。

=== Twin Delayed DDPG (TD3)

类似之前对 Double Q-Learning 的讨论，DPG 也有一样的过估计问题：

$
EE [max_theta Q (s, bold(pi)_theta (s))] >= max_theta EE [Q(s, bold(pi)_theta (s))]
$

TD3 是表现最好的 DDPG 变体之一，其对 DDPG 作了一些补充。首先是 TD3 中的 "Twin"，引入 Clipped Double Q-Learning（CDQ），考虑维护两个 critic Q 网络，但和 Double Q-Learning 一个选取最优动作一个评估价值（构造目标）的思路不同，TD3 中这两个网络是平等的，使用时取其中较小的值，所以其实是一种保守估计方法，或者说 pessimistic 的。

第二是 regularization with clipped noise，在更新 Q 网络时最小化 loss：

$
cal(L)_i^Q [phi.alt_i] = EE_mu [sum_(t=0)^(n-1) (R_t + gamma min_(j in {1,2}) Q_(phi.alt'_j) (S_(t+1), bold(pi)_(theta') (S_(t+1)) + bold(epsilon.alt)) - Q_(phi.alt_i) (S_t, S_t))^2]
$

其中 $bold(epsilon.alt) ~ "clip"(cal(N) (bold(0), bold(sigma)^2), -bold(c), bold(c)) in RR ^abs(cal(A))$，除了多了一个 Q 网络，区别就在于在策略上加了一个噪声（clip 以防太大），以平均掉一些确定策略下可能的单点尖峰错误的影响。

第三是 TD3 中的 "Delayed"，critic 和 actor 不再同频率更新，而是 critic 更新多次后才更新一次 actor。这是由于 actor 的 loss 依赖 critic 的值函数估计，此举尝试让 critic 估计更稳定后再更新 actor。

在谱图 @fig:policy_opt_spectrum 中，DDPG 和 TD3 属于 off-policy 确定策略算法，critic 估计动作值函数，并且采用 pessimism 估计稳定 Q 值估计。

== Robust Decision Making

#Cre("TODO")

=== Regularization

=== Robust Reinforcement Learning

=== Planning as Inference

=== Stochastic Policy Gradients (SPG) and Deterministic Policy Gradients (DPG)

=== Reparameterization Trick



== Maximum Entropy RL



=== Soft Actor-Critic (SAC)

SAC 可以基于 TD3 的框架实现，具体地：1、使用重放缓冲、pessimism，但不用 clipped noise；2、
