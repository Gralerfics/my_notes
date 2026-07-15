#import "../generic.typ": *

#import "@preview/codly:1.3.0": *

#show: codly-init.with()

= Exploration

*探索*（exploration）是数据的来源，之前我们已经使用过一些探索技巧，例如 *$epsilon.alt$-贪婪策略*，即以 $epsilon.alt$ 的概率（一般会随探索逐渐下调直至一个下限）随机选择动作，其余情况按目标策略采样：

$
pi_theta^(epsilon.alt) (a mid(|) s) := (1 - epsilon.alt) pi_theta (a mid(|) s) + epsilon.alt 1/abs(cal(A))
$

还有基于价值函数的 *Boltzmann 探索*（uninformed Boltzmann exploration）：

$
pi_theta^beta (a mid(|) s) := exp(beta Q_theta (s, a)) / (sum_(a') exp (beta Q_theta (s,a')))
$

当 $beta = 1$ 时，就是用 softmax 函数把价值函数表示成了动作概率分布；当 $beta -> 0$ 时，各动作概率趋于相似，整体策略趋于纯随机探索；当 $beta > 1$ 时，本身价值函数较大的动作权重也更大，策略趋于只选择最大 $Q$ 值动作的确定性策略。在 @equ:maximum_entropy_reg_entropy 还提到过*最大熵正则化*，基于损失函数，保持策略的随机性。

优化的探索策略是相对纯随机探索来说的，那么#underline[怎样的探索策略算是更智能的]？一个自然的选择是#underline[优先选择已经观察到高回报的动作]，即围绕当前最优策略（exploitation policy）进行探索。毕竟过去的经验认为某个动作价值比较高，那么之后也可以更频繁地采用它。前述 $epsilon.alt$-贪婪策略和 Boltzmann 探索策略都属于这一类，在当前最优下引入少量随机性。这类策略也天生存在缺陷，容易过快收敛于短期奖励指示的最优奖励（前期容易达到的奖励），并不断巩固这一策略最终错过可能存在的长期高回报。

另一种是从广泛探索的角度出发，#underline[优先选择回报确定性不高的动作]。以基于价值函数估计的算法为例，价值函数估计的准确性依赖于在这部分状态空间上探索的充分性；回报估计波动大，说明探索不充分。特别是对于 Q-Learning 这样的 off-policy 算法，优先探索还没有被充分探索、评估的区域可以加速价值函数估计的收敛，辅助智能体更好地做出长期判断。不过相较值迭代等算法中可以通过判断 "更新后价值函数是否还改变" 来确定是否已经收敛，这类策略是通过回报样本来衡量不确定性的。这就导致如果系统本身就有随机性，称为*偶然不确定性*（aleatoric uncertainty），这类策略就容易陷进去出不来。具体一点，比如说系统在某处给予回报本来就是随机的（比如投个骰子决定回报），再多的探索也不会降低这里回报的不确定性，而这类探索策略还是会义无反顾地选择它，无法脱身。

既然客观的偶然不确定性无法降低（irreducible），那么就让智能体更主观地评估 "不确定性"，#underline[优先选择还没有尝试过的动作]，就不会陷入循环。动作是否已经尝试过完全是从智能体的视角来说的，未尝试过的动作表示其对应状态空间尚未被完全探索，这种不确定性可称为*认知不确定性*（epistemic uncertainty），不会像偶然不确定那样因系统性质导致无法降低。但这类策略也有问题，一是难以处理连续动作空间的情况，二是会导致过度探索——本来智能的策略就是希望用有限的算力探索过于庞大的状态空间中的低维有效信息，现在这样就有退回到暴力搜索的意思了。

== Uncertainty

前面提到了无法降低的偶然不确定和可以降低的认知不确定性，都属于各种因素导致的不确定性（uncertainty）。除此之外，如果是基于模型的方法，错误的模型也会导致不可降低的模型偏差（model bias）。

均方误差是非常常用的损失函数，假设在数据分布 $cal(D)$ 上计算，$X~rho(dot)$, $Y~cal(N)(dot mid(|) mu(X), sigma^2 (X))$，我们可以将 MSE 分解为几个部分：#Cre("TODO")

$
&EE_(cal(D)) underbrace([EE [(Y - f_(cal(D)) (x))^2 mid(|) cal(D)]], "generalization error") \
=& underbrace(EE [sigma^2 (X)], "aleatoric") + underbrace(EE [VV_(cal(D)) [f_(cal(D)) (X) mid(|) X]], "epistemic") + underbrace(EE [(EE_(cal(D)) [f_(cal(D)) (X) mid(|) X] - mu(X))^2], "model-bias")
$

#Cre("TODO")

== Thompson Sampling

== Optimistic Exploration
