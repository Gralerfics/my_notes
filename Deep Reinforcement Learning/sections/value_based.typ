#import "../generic.typ": *

#import "@preview/cetz:0.4.2"
#import "@preview/cetz-plot:0.1.3"
#import "@preview/fletcher:0.5.8" as fletcher: diagram, node, edge

#import "@preview/codly:1.3.0": *

#show: codly-init.with()

= Value-Based Learning

== Policy Iteration and Value Iteration

（略）以迭代的方式求解最优 $V^* (s)$，在其基础上贪心制定最优策略。

// 课程跳过这个，毕竟不 deep。如果要写可写到 value iteration 收敛性的证明，涉及 Bellman 算子压缩映射和不动点。

== Tabular Q-Learning and SARSA

既然可以迭代求解 $V^* (s)$，索性直接求 $Q^* (s, a)$，更方便求解策略 $a^* = arg max_a Q^* (s, a)$。考虑状态和动作空间离散的情况，$Q (s, a)$ 实际上是表格式的（tabular），包含共 $abs(cal(S)) times abs(cal(A))$ 个元素。

过程类似值迭代，先初始化：

$
Q (s, a) <- 1, quad forall s, a
$

让智能体在环境中进行探索，用*有限差分*（temporal difference，TD）目标进行更新。具体地，任意 $t$ 时刻发生状态转移时得到一组新数据 $(s_t, a_t, r_t, s_(t+1))$，定义目标：

$
y_t := r_t + gamma max_a Q (s_(t+1), a)
$

此处#underline[*目标*（target）表示对真实最优动作值函数 $Q^* (s_t, a_t)$ 的估计]，用它搭配一个学习率 $alpha$ 对当前 $Q (s_t, a_t)$ 进行更新：

$
Q (s_t, a_t) <- (1 - alpha) Q (s_t, a_t) + alpha y_t
$ <equ:tab_ql_update_interp>

这种形式形似 “插值”，或者另一种更像 “更新” 的等价形式也比较常用：

$
Q (s_t, a_t) <- Q (s_t, a_t) + alpha (y_t - Q (s_t, a_t))
$ <equ:tab_ql_update_residue>

其中 $y_t - Q (s_t, a_t)$ 被称为*有限差分误差*（TD error）。如此迭代，最终将会有 $Q (s, a) -> Q^* (s, a)$，证明暂略。

值得注意的是，智能体探索时所实际执行的策略即*行为策略*（behavior policy）是基于 $Q$ 估计值计算的，而更新式中 $max_a Q (s_(t+1), a)$ 说明更新时所用的策略即*目标策略*（target policy）是取当前估计下的最优策略。故在 Q-Learning 中，行为策略和目标策略不一致，这种类型的算法称为 off-policy 算法。

与 Q-Learning 相对应，自然也有 on-policy 的 SARSA 算法，其每次状态转移采样得到 $(s_t, a_t, r_t, s_(t+1), a_(t+1))$（这五项的首字母也是算法得名原因），选取的目标为：

$
y_t^"SARSA" := r_t + gamma Q (s_(t+1), a_(t+1))
$

该目标直接计算自实际样本，行为策略和目标策略一致，是 on-policy 的。

表格式 Q-Learning 最大的问题在于*维度诅咒*（curse of dimensionality），状态和动作空间较大时 $Q$ 的元素数量指数级上升，让学习几乎成为不可能，也难以处理状态或动作连续的情况。

== Value Approximation

为了表格式算法 cannot scale 的问题，我们可以将 $V$ 或 $Q$ 函数参数化，用相对少的参数去表达函数值：

$
v_theta (s) approx V^pi (s), quad q_theta (s, a) approx Q^pi (s, a)
$

其中 $theta$ 表示参数。参数化的方式有很多，例如 SoB（sum of basis）、神经网络等。

参数化后的值函数自然无法再用表格式算法的方式进行更新，通常我们采用梯度更新的方式对参数进行更新。于是我们现在需要考虑的问题变为：1、用什么*损失函数*（loss function）；2、用什么计算*梯度*（gradients）。

首先我们规定一些符号，每次发生状态转移的相关信息都是一条 TD 样本：

$
(s_t, a_t, r_t, s_(t+1), a_(t+1))
$

采用策略 $pi$ 采样得到的一条或多条轨迹中收集到的 $n$ 条 TD 样本构成数据集，记为：

$
cal(D) = {(s_t, a_t, r_t, s_(t+1), a_(t+1))}_(t=1)^n
$

=== Residual Value Gradients

最常用的损失函数仍然是*均方误差*（mean-squared error，MSE），考虑值函数估计 $v_theta (s)$ 和目标之间在整个数据集 $cal(D)$ 上的均方误差：

$
cal(L)[theta] := EE_cal(D) [(overbrace(underbrace(r_t + gamma v_theta (s_(t+1)), "target" y_t) - v_theta (s_t), "TD-error"))^2]
$ <equ:loss_vfunc_rvg>

如果用动作值函数则为：

$
cal(L)[theta] := EE_cal(D) [(overbrace(underbrace(r_t + gamma q_theta (s_(t+1), a_(t+1)), "target" y_t) - q_theta (s_t, a_t), "TD-error"))^2]
$

以上都是 on-policy 的，而后者实际上就是 SARSA。这一套损失函数设计与梯度更新方式称为*残差梯度*（residual gradients）法，平方项内的部分称为*残差*（residues）。

残差的形式和 @equ:tab_ql_update_residue 中 $y_t - Q(s_t, a_t)$ 也是一致的，就如前文所述，$y_t$ 是对最优动作值函数的近似估计，#underline[残差本质上就是在衡量当前估计的 $Q(s_t, a_t)$ 是否满足最优贝尔曼方程] @equ:q_bellman_optimal，优化该损失函数就是在尝试推动 $q_theta$ 收敛到 $Q^*$。

但要注意，目标 $y_t$ 中用到了 $q_theta (s_(t+1), a_(t+1))$，它也是基于当前参数 $theta$ 的一个#underline[近似值而不是真实的] $Q^* (s_(t+1), a_(t+1))$，这导致它实际无法真的令 $q_theta$ 收敛到 $Q^*$。类似这种 “用估计的参数当目标去优化参数估计” 的行为我们称为*自举*（bootstrapping），也是引发问题的根源。细节上的不同使得这种自举迭代可能让结果收敛（例如值迭代和 Q-Learning），也可能导致发散或产生稳态误差。

接下来具体讨论这个问题。为了更清晰地展示参数的变化，我们以残差梯度法更新表格式 Q-Learning 为例进行说明。它也可以视作一种参数化，只不过参数就是所有表格值 $theta := bold(Q) = {Q_(s,a)} in RR^(abs(cal(S)) times abs(cal(A)))$。写出 Q-Learning 观察到转移 $(s_t, a_t, r_t, s_(t+1))$ 后的一步损失：

$
cal(L)[bold(Q)] := (r_t + gamma Q_(s_(t+1), a^*) - Q_(s_t, a_t))^2, quad a^* := arg max_a Q_(s_(t+1), a)
$

及一步梯度下降更新式：

$
bold(Q) <- bold(Q) - alpha/2 nabla_bold(Q) cal(L)[bold(Q)]
$

由于参数就是表格元素，各元素对 $bold(Q)$ 求导就是对应位置的单位冲激函数：

$
nabla_bold(Q) Q_(s, a) = delta_(s,a), quad forall s in cal(S), forall a in cal(A)
$

其中 $delta_(s,a) = 1$ 当且仅当自变量为 $(s, a)$，课件则采用的类似 $delta (a = dot, s = dot)$ 的写法。于是表格每项的更新为：

$
Q_(s, a) <- Q_(s, a) - alpha (r_t + gamma Q_(s_(t+1), a^*) - Q_(s_t, a_t)) (gamma delta_(s_(t+1), a^*) - delta_(s_t, a_t)), quad forall s in cal(S), forall a in cal(A)
$

由于其中两个 $delta$ 函数的存在，这里的一步更新实际上只会更新两个元素：

$
Q_(s_t, a_t) &<- (1 - alpha) Q_(s_t, a_t) + alpha (r_t + gamma Q_(s_(t+1), a^*)) \
Q_(s_(t+1), a^*) &<- (1 - alpha gamma^2) Q_(s_(t+1), a^*) - alpha gamma (r_t - Q_(s_t, a_t))
$

第一部分 $Q_(s_t, a_t)$ 的更新没有什么问题，$r_t + gamma Q_(s_(t+1), a^*)$ 就是 $y_t$，和之前 @equ:tab_ql_update_interp 中完全一致；但这里还多出了第二部分 $Q_(s_(t+1), a^*)$ 的更新，与前文标准的表格式 Q-Learning 不一致了。

这会导致两个问题，#underline[首先最重要的是多了这一项的残差梯度是*有偏的*（biased），无法再令参数收敛到正确值]，印证之前的猜想。因为利用贝尔曼最优算子我们证明 $Q^*$ 是最优贝尔曼方程的唯一不动点解，且标准表格式 Q-Learning 会收敛到这个解，但现在残差梯度化简后与其存在偏差，自然收敛结果也有偏差。

第二是#underline[违背了*因果性*（causality）]。具体地，更新的本质是树立一个目标对当前参数 $Q_(s_t, a_t)$ 进行调整，但残差梯度同时还去更新了未来状态的参数 $Q_(s_(t+1), a^*)$。这种#underline[参数估计和目标同时在变]的情况自然会导致更新过程不稳定，参数信息传播和学习速度变慢。

=== Value Semi-gradients

为了解决前述问题，标准 Q-Learning 常采用的是*半梯度*（semi-gradients）更新的方式。比如上节例子中残差梯度法同时更新 $Q_(s_t, a_t)$ 和 $Q_(s_(t+1), a^*)$，半梯度法则只用 $Q_(s_t, a_t)$ 这一部分梯度更新，从而和表格式 Q-Learning 保持一致并得以套用其收敛性结论。

回到参数化损失 @equ:loss_vfunc_rvg，问题本质是目标 $y_t$ 依赖当前参数估计 $theta$，那么我们可以通过将目标中的 $v_theta$ 替换为另一组迟滞参数驱动的 $v_(theta')$ 来削减这种耦合：

$
cal(L)[theta] := EE_cal(D) [(underbrace(r_t + gamma v_(theta') (s_(t+1)), "bootstrapped target" y_t) - v_theta (s_t))^2]
$

这样在更新时就有 $nabla_theta v_(theta') (s_(t+1)) = 0$，对应舍去的那一部分梯度。

实际更新时 $theta'$ 的数值是和 $theta$ 相等的（$theta' = theta$），所以说明白点就只是反向传播时切断了 $theta'$ 的这一部分梯度，前向传播时照旧。

以上更新方式称为*半梯度有限差分学习*（semi-gradient TD-learning）。在 PyTorch 具体实现中，给相应项加上 `.detach()` 即可实现只截断某一部分反向传播的操作。该方案比残差梯度法更快且无偏，对#underline[线性模型可以证明收敛]，例如前文表格 Q-Learning 的例子。

但#underline[采用神经网络等模型参数化值函数时就很容易发散]，这是由于模型估计误差、自举性等因素复杂综合导致的。引入 $theta'$ 只能弱化耦合而无法消除，因为终究无法得到真实的 $Q^*$ 充当优化目标。

若要进一步可以考虑*神经拟合 Q 迭代*（neural-fitted Q-iteration，NFQ）算法，这是第一个成功的深度强化学习算法，就是迭代比较慢。相比普通的半梯度学习，它会在每次阶段性收敛后再将 $theta$ 同步到 $theta'$，一步一个脚印显著降低发散概率。

=== Online Q-Learning

接下来讨论一个实例的具体实现细节，考虑用半梯度 TD 学习实现 Q-Learning，神经网络作为参数化模型。这是最基础的深度 Q-Learning，为后文 DQN 做铺垫。

首先有关参数化模型 $q_theta (s, a)$，自然的设计是将 $s$ 和 $a$ 作为输入层，输出层输出 $Q$ 值。但对于离散动作 Q-Learning，更新时需要反复计算 $max_a q_theta (s_(t+1), a)$，如果能一次前向传播计算出所有动作 $a$ 对应的 $q_theta (s_(t+1), a)$ 值（输出一个向量）相比依次计算最后取最值自然是更高效的。所以我们将状态作为输入层，输出层则让每个动作对应一个输出（one output per action）。

引入半梯度法后，模型可以表达为：

$
underbrace(q_theta (s_t, a_t), "value") =^! max_pi EE_pi [sum_(k=0)^infinity gamma^k r_(t+k) mid(|) s_t, a_t] approx underbrace(r(s_t, a_t) + gamma EE[max_(a') q_(theta') (s_(t+1), a')], "target")
$

等号上的叹号表达我们希望二者相等，即希望 $q_theta$ 是 $Q^*$ 的良好估计。用 PyTorch 实现节选如下：

#codly()
```python
from torch.optim import RMSprop
from torch.nn.functional import mse_loss
optimizer = RMSprop(q.parameters())
while True:
    # sample episode from environment
    state, action, reward, term, next_state = environment.sample(q)
    # compute left and right side of Bellman eq.
    value = q(state).gather(dim = -1, index = action)
    target = reward + gamma * (~term * q(next_state).max(dim = -1)[0])
    # gradient descent step on supervised regression loss
    optimizer.zero_grad()
    mse_loss(value, target.detach()).backward()
    optimizer.step()
```

其中，可以观察到半梯度法 `target.detach()` 的使用。此外，代码中的 `term` 指示是否是终止状态（terminal state），确保在终止状态不再加入未来的奖励。

== Stabilization Issues and Techniques

前文讨论过，朴素的半梯度在线 Q-Learning 很容易发散。从机器学习的角度来看，这种不稳定性源于其违背了一系列假设，例如：
+ 回归目标不平稳（non-stationary），即前述自举性的问题；
+ 贪婪策略改变了训练分布；
+ 网络会遗忘旧样本；
+ 状态转移不是独立同分布的；

TODO

=== On-policy and Off-policy Sampling

在强化学习中没有传统监督学习的训练集和测试集划分环节，所有数据来自与环境交互的采样。

=== Experience Replay

TODO

如何利用 Q-Learning off-policy 的特性？

*重放缓存*（experience replay buffer）：

优先重放缓存（prioritized experience replay buffer）

=== Exploration and Neural Networks

TODO

$epsilon.alt$-贪婪探索策略（$epsilon.alt$-greedy exploration policy）

=== Target Networks

TODO

半梯度 Q-Learning

硬目标更新

软目标更新

== Deep Q-Networks (DQN)

TODO
