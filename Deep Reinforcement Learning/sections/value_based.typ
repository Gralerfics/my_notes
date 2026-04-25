#import "../generic.typ": *

#import "@preview/cetz:0.4.2"
#import "@preview/cetz-plot:0.1.3"
#import "@preview/fletcher:0.5.8" as fletcher: diagram, node, edge

= Value-Based Learning

== Policy Iteration and Value Iteration

（略）以迭代的方式求解最优 $V^* (s)$，在其基础上贪心制定最优策略。

// 课程跳过这个，毕竟不 deep。如果要写可写到 value iteration 收敛性的证明，涉及 Bellman 算子压缩映射和不动点。

== Tabular Q-Learning and SARSA

既然可以迭代求解 $V^* (s)$，索性直接求 $Q^* (s, a)$，更方便求解策略 $a^* = arg max_a Q^* (s, a)$。考虑状态和动作空间离散的情况，$Q (s, a)$ 实际上是表格式的（tabular），包含共 $abs(cal(S) times cal(A))$ 个元素。

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
$

这种形式形似 “插值”，或者另一种更像 “更新” 的等价形式也比较常用：

$
Q (s_t, a_t) <- Q (s_t, a_t) + alpha (y_t - Q (s_t, a_t))
$

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
$



以上都是 on-policy 的，而后者实际上就是 SARSA。

接下来我们以表格式 Q-Learning 为例讨论一下 residual gradients 关于因果性的问题。

TODO

=== Value Semi-gradients



// *端到端学习*（end-to-end learning）
