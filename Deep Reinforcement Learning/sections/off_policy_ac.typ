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

要实现 DPG 需要维护两个网络，一个是作为确定策略的 actor 网络 $bold(pi)_theta (s)$，一个是作为动作值函数估计的 critic 网络 $Q_(phi.alt) (s, bold(a))$。二者分别优化各自的损失函数，actor 优化：

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

=== Twin Delayed DDPG (TD3)



== Robust Decision Making

== Maximum Entropy RL

=== Soft Actor-Critic (SAC)
