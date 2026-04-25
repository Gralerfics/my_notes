#import "../generic.typ": *

#import "@preview/cetz:0.4.2"
#import "@preview/fletcher:0.5.8" as fletcher: diagram, node, edge

= Fundamentals

== Markov Chains

我们将智能体的行为建模为一个*随机策略*（stochastic policy）$pi (a mid(|) s)$，即在状态为 $s$ 的条件下采取动作 $a$ 的概率分布。*状态*（state）是一种对当前已知信息的描述，例如智能体所处环境、内部状态等；*动作*（action）则是智能体可选的行动，采取后可能导致状态发生变化，转移到另一个状态，即*状态转移*（state transition）。

在强化学习的背景下，我们会对 “在某状态下采取某动作” 的行为发放*奖励*（reward）以引导智能体的学习。奖励的发放也可以选取其他时机，例如采取动作后根据实际到达的新状态发放奖励等，不同约定互相之间通常都可转换，只是定义和符号差别。此处我们约定，#underline[某时刻的奖励仅为该时刻状态和动作的函数]。

随时间发展状态不断转移，智能体和环境将演化出一条由状态、动作和奖励构成的随机序列：

$
Tau_n = {S_0, A_0, R_0, S_1, A_1, R_1, dots, S_(n-1), A_(n-1), R_(n-1), S_n}
$

序列的各个元素为随机变量，这是由于我们的策略 $pi$ 和状态转移模型 $P$ 都将使用随机分布建模。严格来说，由前面的约定，我们采用的奖励函数是关于状态和动作的确定（deterministic）函数，故 $R_t$ 并不是随机变量（或者视为只有一种概率为 $1$ 的取值）。如果需要，奖励函数自然也可以用概率模型建模，用以刻画奖励的随机性；但奖励终究是人为设计的，此举通常没有特别的必要。

从该随机序列进行一次采样即得到一条*轨迹*（trajectory）：

$
tau_n := {s_t, a_t, r_t}_(t=0)^(n-1) union {s_n}
$

用小写字母表示样本以作区分。$n$ 为时域窗口的长度，可以有 $n -> infinity$。$s_n$ 为终端（terminal）状态，可以用 “在 $n$ 时刻后不再有奖励” 来表达，即 $r_t = 0, forall t >= n$。

如果#underline[状态转移概率 $P$ 仅和当前时刻的状态和动作有关]，即满足马尔可夫性：

$
&P(S_(t+1) mid(|) S_t = s_t, A_t = a_t, S_(t-1) = s_(t-1), A_(t-1) = a_(t-1), dots) \
=&P(S_(t+1) mid(|) S_t = s_t, A_t = a_t)
$

则前述序列可视作一条*马尔可夫链*（Markov chain），其演化过程称为*马尔可夫决策过程*（Markov decision process，MDP），记为 $cal(M) = chevron.l cal(S), cal(A), rho, P, r chevron.r$，其中：
- $cal(S)$ 为状态集合；
- $cal(A)$ 为动作集合；
- $rho$ 为初始状态先验分布，即有 $S_0 ~ rho (dot)$；
- $P$ 为状态转移概率，即有 $S_(t+1) ~ P(dot mid(|) S_t = s_t, A_t = a_t)$；
- $r$ 为奖励函数，按我们的约定有 $R_t equiv r_t = r(s_t, a_t)$。

每时刻采取的动作则由随机策略 $pi$ 描述，即有 $A_t ~ pi (dot mid(|) S = s_t)$。

#blockquote([
  *关于随机变量与样本的符号约定*：

  规范一些来说，我们可以描述一个随机变量符合某个分布，而样本是一个确定值，比如状态转移概率 $P$ 严格来说应写成：

  $
  S_(t+1) ~ P(dot mid(|) S_t = s_t, A_t = a_t)
  $

  但它们在含义上都指向同一个量，通常不影响表达，故为了简便起见，没有歧义时我们会#underline[统一用小写的符号来表达随机变量和样本的含义]，例如：
  
  $
  s_(t+1) ~ P(dot mid(|) s_t, a_t)
  $

  其中条件概率中的随机变量也可能省略，例如 $s_t$ 能看出表示 $t$ 时刻状态就不再写 $S_t = s_t$；如果换其他符号看不出时，就可能显式写 $S_(t+1) = s'$ 或 $s_(t+1) = s'$。
])

#figure(
  image("../figures/markov_chain_flow.png", width: 80%),
  caption: [马尔可夫链示意图（来自课件）。]
) <fig:markov_chain_flow>

总而言之，MDP 加上策略 $pi$ 这#underline[两个要素]构成了马尔可夫链，如图 @fig:markov_chain_flow 所示。因此，随机轨迹序列 $Tau$ 的分布参数依赖 MDP 参数和 $pi$，可记为：

$
Tau ~ p_Tau (dot; cal(M), pi)
$

我们用#underline[分号表达 “参数” 的含义]，和条件概率中的竖线不同；有时我们也将参数写到上标之类的地方，理解即可。严格来说，$n$ 也应该写在 MDP 参数里，不过我们主要考虑 $n -> infinity$，通常就不写了。

在我们的讨论范畴内#underline[通常默认 MDP 模型]（部分证明涉及给定模型下不同策略的比较），所以可以只标记 $pi$ 用以区分，类似地还有之后的 $EE_pi$、$V^pi$ 等记号。再加上刚刚约定的关于随机变量和样本符号的规定，之后会出现如下形式表达：

$
tau ~ p_Tau^pi (dot)
$

== Return and Value Functions

一条轨迹的总体奖励称为*回报*（return），通常我们会每时刻乘上*折扣因子*（discount）来刻画 “长期奖励贡献递减” 的特性，不考虑太远的未来，但也不能只看眼前的奖励，再考虑 $n -> infinity$，定义*折扣回报*（discounted return）：

$
J (tau) = sum_(t=0)^infinity gamma^t r_t
$

折扣回报是基于轨迹的，而轨迹是随机的，为了衡量一个状态的价值，我们需要的是回报的期望：

$
V^pi (s) := EE_(tau~p_Tau^pi (dot)) [sum_(t=0)^infinity gamma^t r_t mid(|) s_0 = s]
$ <equ:vfunc_def>

这个期望回报称为*值函数*（value function），表示在使用策略 $pi$ 的前提下从 $s$ 出发获得的平均折扣回报，加了条件 $s_0 = s$ 表示#underline[对给定策略下某个状态平均价值的一种衡量]。应用全概率公式引入初态先验分布 $rho$ 即可得到#underline[对整个策略平均价值的衡量]：

$
J^pi := EE_(tau~p_Tau^pi (dot)) [sum_(t=0)^infinity gamma^t r_t] = EE_(s_0~rho(dot)) [V^pi (s_0)]
$ <equ:j_def_and_jv_relation>

为了更好地研究策略，我们有时还希望#underline[衡量在某状态下采取某动作的平均价值]，故引入*动作值函数*（action-value function）：

$
Q^pi (s, a) := EE_(tau~p_Tau^pi (dot)) [sum_(t=0)^infinity gamma^t r_t mid(|) s_0 = s, a_0 = a]
$ <equ:qfunc_def>

由此可以看出 $J^pi$、$V^pi$ 和 $Q^pi$ #underline[在条件上是递进的关系]，@equ:j_def_and_jv_relation 中描述了前两者的关系，对后两者则有全概率公式：

$
V^pi (s) &= EE_(a~pi(dot mid(|) s)) [Q^pi (s, a)] \
&= sum_a pi(a mid(|) s) Q^pi (s, a)
$ <equ:vq_relation>

从定义出发可以得到逆过来的关系（证明见后）：

$
Q^pi (s, a) &= r(s, a) + gamma EE_(s'~P(dot mid(|) s, a)) [V^pi (s')] \
&= r(s, a) + gamma sum_(s') P(s' mid(|) s, a) V^pi (s')
$ <equ:qv_relation>

#blockquote([
  *关于值函数与定义中初态时间索引无关的证明*：

  值函数的定义中是从 $s_0 = s$ 开始计算折扣回报，那由于马尔可夫性，它和从任意 $s_k = s$ 开始计算其后（前面的不算）的回报理应是一样的，即有：

  $
  EE_(tau~p_Tau^pi (dot)) [
    sum_(t=0)^infinity gamma^t r_(t+k)
    mid(|) s_k = s
  ]
  =
  EE_(tau~p_Tau^pi (dot)) [
    sum_(t=0)^infinity gamma^t r_t
    mid(|) s_0 = s
  ]
  =: V^pi (s)
  $ <equ:vfunc_time_shift_prop>

  只是时间轴平移了 $k$ 步，直觉上很好理解，但还是需要形式上证明一下，后面会用到。为方便展开取前 $n$ 步，左式为：

  $
  &EE_(tau~p_Tau^pi (dot)) [
    sum_(t=0)^n gamma^t r_(t+k)
    mid(|) s_k = s
  ] \

  = &sum_(a_k, s_(k+1), ..., a_(k+n))
  (
    sum_(t=0)^n gamma^t r(s_(t+k), a_(t+k))
  )product_(t=0)^n
  pi(a_(t+k) mid(|) s_(t+k)) \
  &P(s_(t+k+1) mid(|) s_(t+k), a_(t+k)),
  quad s_k = s
  $

  令 $tilde(s)_t = s_(t+k), tilde(a)_t = a_(t+k)$ 则有 $tilde(s)_0 = s_k = s$，于是上式化为：

  $
  sum_(tilde(a)_0, tilde(s)_1, ..., tilde(a)_n)
  (
    sum_(t=0)^n gamma^t r(tilde(s)_t, tilde(a)_t)
  )
  product_(t=0)^n
  pi(tilde(a)_t mid(|) tilde(s)_t)
  P(tilde(s)_(t+1) mid(|) tilde(s)_t, tilde(a)_t)
  $

  而另一方面，直接从 $s_0 = s$ 出发，右式为：

  $
  &EE_(tau~p_Tau^pi (dot)) [
    sum_(t=0)^n gamma^t r_t
    mid(|) s_0 = s
  ] \

  = &sum_(a_0, s_1, ..., a_n)
  (
    sum_(t=0)^n gamma^t r(s_t, a_t)
  )
  product_(t=0)^n
  pi(a_t mid(|) s_t)
  P(s_(t+1) mid(|) s_t, a_t),
  quad s_0 = s
  $

  二者完全相同，只是符号不同。因此对任意 $n$ 都有：

  $
  EE_(tau~p_Tau^pi (dot)) [
    sum_(t=0)^n gamma^t r_(t+k)
    mid(|) s_k = s
  ]
  =
  EE_(tau~p_Tau^pi (dot)) [
    sum_(t=0)^n gamma^t r_t
    mid(|) s_0 = s
  ]
  $

  令 $n -> infinity$，若 $0 <= gamma < 1$ 且奖励有界，则折扣和收敛，则证得 @equ:vfunc_time_shift_prop。
])

#blockquote([
  *关于 “动作值函数-值函数” 递推关系的证明*：

  我们来证明 @equ:qv_relation。首先：

  $
  Q^pi (s, a) &=^"def" EE_(tau~p_Tau^pi (dot)) [sum_(t=0)^infinity gamma^t r_t mid(|) s_0 = s, a_0 = a] \
  &=^((1)) EE_(tau~p_Tau^pi (dot)) [r_0 + gamma sum_(t=0)^infinity gamma^t r_(t+1) mid(|) s_0 = s, a_0 = a] \
  &=^((2)) r(s, a) + gamma EE_(tau~p_Tau^pi (dot)) [sum_(t=0)^infinity gamma^t r_(t+1) mid(|) s_0 = s, a_0 = a]
  $

  式中第 $(1)$ 步拆出第一项 $r_0$ 和系数 $gamma$ 并调整索引（$t <- t+1$），第 $(2)$ 步代入 $r_0 = r(s_0, a_0)$ 并将其和系数从期望中提出来。表达式应用全概率公式可得：

  $
  r(s, a) + gamma sum_(s') P(s' mid(|) s, a) EE_(tau~p_Tau^pi (dot)) [sum_(t=0)^infinity gamma^t r_(t+1) mid(|) s_0 = s, a_0 = a, s_1 = s']
  $

  由马尔可夫性，给定 $s_1 = s'$ 条件的情况下未来与 $(s_0, a_0)$ 无关，故消去相应条件得：

  $
  Q^pi (s, a) = r(s, a) + gamma sum_(s') P(s' mid(|) s, a) EE_(tau~p_Tau^pi (dot)) [sum_(t=0)^infinity gamma^t r_(t+1) mid(|) s_1 = s']
  $

  再利用时移性质 @equ:vfunc_time_shift_prop 将期望替换为 $V^pi (s')$ 即可得 @equ:qv_relation。
])

== Bellman Expectation Equations

我们可以推导出关于值函数的递推表达式：

$
V^pi (s) &= EE_(a~pi (dot mid(|) s), s'~P(dot mid(|) s, a)) [r(s, a) + gamma V^pi (s')] \
&= sum_(a) pi(a mid(|) s) sum_(s') P(s' mid(|) s, a) [r(s, a) + gamma V^pi (s')]
$ <equ:v_bellman_expectation>

该递推式也被称为*贝尔曼期望公式*（Bellman expectation equation），证明见后。动作值函数的版本是：

$
Q^pi (s, a) &= r(s, a) + gamma EE_(s'~P(dot mid(|) s, a), a'~pi (dot mid(|) s')) [Q^pi (s', a')] \
&= r(s, a) + gamma sum_(s') P(s' mid(|) s, a) sum_(a') pi(a' mid(|) s') Q^pi (s', a')
$ <equ:q_bellman_expectation>

#blockquote([
  *关于期望与条件期望的符号约定*：

  在课件中，值函数是这样定义和展开的：

  $
  V^pi (s) := EE_pi [sum_(t=0)^infinity gamma^t r_t mid(|) s_0 = s] = EE [r(s, a) + gamma V^pi (s') mid(|) inline(mat(delim: #none, a~pi (dot mid(|) s); s'~P(dot mid(|) s, a)))]
  $ <equ:vfunc_def_and_expect_equ_slides_ver>

  // 其记号是有些混乱的，但如果处处都写完整比较繁琐，我们希望表意到位即可。于是我们在此澄清相关约定，之后看情况选取合适的形式书写。
  
  首先，形如 $s_0 = s$ 的条件混淆了随机变量和样本，按照我们之前的定义应该写 $S_0 = s$ 用以区分。不过如前注释所述，为简化符号我们#underline[接受这种混合的写法]。
  
  第二，求期望需要明确的是随机变量及其采样自什么分布。值函数定义中条件期望的条件约束了初始状态 $s_0$，求和式中则使用了各时刻的奖励 $r_t$，比起列举所有变量，我们可以直接令这里的随机变量是服从轨迹概率分布的随机序列 $Tau$，即前述 $Tau ~ p_Tau (dot; cal(M), pi)$，其包含了所有的状态、动作和奖励，对应 @equ:v_bellman_expectation 中 $EE_(tau~p_Tau^pi (dot)) [dot]$ 的写法。#underline[此处我们接受简化]，可写为 $EE_pi [dot]$，类似地还有用 $EE_rho [dot]$ 代表 $EE_(s_0~rho(dot)) [dot]$ 等。
  
  第三，一般不会像第二个等号这样将随机变量和及其分布写到条件期望条件的位置上去，而是写在 $EE$ 下标处或者直接省略。有时考虑公式长度以及手写的困难，也#underline[接受这种写在条件处的做法]。

  如此我们就统一了 @equ:v_bellman_expectation 和 @equ:vfunc_def_and_expect_equ_slides_ver 的形式。接下来补充部分其他公式在课件里的写法，方便对照：

  $
  J^pi := EE_pi [sum_(t=0)^infinity gamma^t r_t] &= EE_rho [V^pi (s_0)] = EE [sum_(t=0)^infinity gamma^t r(s_t, a_t) mid(|) inline(mat(delim: #none, s_0~rho(dot)", "a_t~pi (dot mid(|) s_t); s_(t+1)~P(dot mid(|) s_t, a_t)))]
  $

  $
  Q^pi (s, a) := EE_pi [sum_(t=0)^infinity gamma^t r_t mid(|) inline(mat(delim: #none, s_0 = s; a_0 = a))] = r(s, a) + gamma EE [Q^pi (s', a') mid(|) inline(mat(delim: #none, s'~P(dot mid(|) s, a); a'~pi (dot mid(|) s')))]
  $

  此后关于这类符号问题不再赘述，将依情况选择合适的方式书写。
])

#blockquote([
  *关于贝尔曼期望公式的证明*：
  
  我们先证明#underline[值函数的版本]：

  $
  V^pi (s) &=^"def" EE_(tau~p_Tau^pi (dot)) [sum_(t=0)^infinity gamma^t r_t mid(|) s_0 = s] \

  &=^((1)) EE_(tau~p_Tau^pi (dot)) [r_0 + gamma sum_(t=0)^infinity gamma^t r_(t+1) mid(|) s_0 = s] \

  &=^((2)) EE_(tau~p_Tau^pi (dot)) [
    r_0 + gamma
    EE_(tau~p_Tau^pi (dot)) [
      sum_(t=0)^infinity gamma^t r_(t+1)
      mid(|) s_0 = s, a_0, s_1
    ]
    mid(|) s_0 = s
  ] \

  &=^((3)) EE_(tau~p_Tau^pi (dot)) [
    r_0 + gamma
    EE_(tau~p_Tau^pi (dot)) [
      sum_(t=0)^infinity gamma^t r_(t+1)
      mid(|) s_1
    ]
    mid(|) s_0 = s
  ] \

  &=^((4)) EE_(tau~p_Tau^pi (dot)) [
    r_0 + gamma V^pi (s_1)
    mid(|) s_0 = s
  ] \

  &=^((5)) EE_(a_0~pi (dot mid(|) s), s_1~P(dot mid(|) s, a_0)) [
    r(s, a_0) + gamma V^pi (s_1)
  ] \

  &=^((6)) EE_(a~pi (dot mid(|) s), s'~P(dot mid(|) s, a)) [
    r(s, a) + gamma V^pi (s')
  ]
  $

  分别阐释一下上式每个步骤的理由：最开始的 $"def"$ 式是定义，接下来第 $(1)$ 步从求和中拆出第一项并调整了一下索引（$t <- t + 1$）。
  
  第 $(2)$ 步是应用迭代期望公式（tower property of conditional expectation）将里面套上条件期望，引入额外条件 $A_0 = a_0$ 和 $S_1 = s_1$，由于之前我们约定不区分随机变量的符号，此处为了防止类似 $a_0 = a_0$ 的歧义就只写 $a_0$。

  第 $(3)$ 步利用马尔可夫性，已知 $t=1$ 时状态为 $s_1$，则条件分布和之前时刻的 $s_0$、$a_0$ 都无关，从而省去条件。第 $(4)$ 步直接利用 @equ:vfunc_time_shift_prop 中的时移性质。

  第 $(5)$ 步是将轨迹分布 $p_Tau^pi (dot)$ 对 $(a_0, s_1)$ 的边缘分布取出来。具体地，轨迹分布是链上所有随机变量的联合分布，而第 $(4)$ 步所求期望仅依赖 $s_0$、$a_0$、$r_0$ 和 $s_1$，故无需考虑之后的 $a_1$ 等；余下的 $s_0 = s$ 和 $r_0 = r(s_0, a_1)$ 由条件和定义消去，只留下对 $(a_0, s_1)$ 的边缘分布。最后第 $(6)$ 步只是符号等价代换（$s_1 -> s'$，$a_0 -> a$），得到更通用的形式。

  这里的证明思路和 @equ:qv_relation 的证明思路基本一致，大致都是#underline[引入下一时刻的随机变量以形成递推的结构]。
  
  有了前面的证明，对于#underline[动作值函数的版本]只需将 @equ:vq_relation 代入 @equ:qv_relation 即可得到 @equ:q_bellman_expectation。
])

== Optimal Bellman Equations

到现在为止，我们定义的值函数和动作值函数都是基于策略 $pi$ 的，通过调整策略可以得到不同的值函数集合。

我们定义*最优值函数*（optimal value function）和*最优动作值函数*（optimal action-value function）：

$
V^* (s) := max_pi V^pi (s)
$

$
Q^* (s, a) := max_pi Q^pi (s, a)
$

注意#underline[该定义并不能保证存在一种最优策略]使得每个状态 $s$ 或状态动作对 $(s, a)$ 都能取到这个最优值。不过我们就#underline[假设存在最优策略] $pi^*$：

$
exists pi^*, quad V^(pi^*) (s) equiv V^* (s) >= V^pi (s), quad forall pi, forall s in cal(S)
$ <equ:optimal_policy_exists_assumption>

动作值函数同理。实际上这个假设对大部分情况都成立，因为我们用 $pi (a mid(|) s)$ 对策略建模，本身就可以针对每个 $s$ 设计策略，而将所有状态下的最优子策略拼在一起即可得到最优策略。

之前推导了值函数和动作值函数的关系，这里也有#underline[最优值函数和最优动作值函数的关系]：

$
V^* (s) &= max_a Q^* (s, a)
$ <equ:vq_opt_relation>

#blockquote([
  *关于最优值函数和最优动作值函数关系的证明*：

  首先由定义和 @equ:vq_relation 有：

  $
  V^* (s) = max_pi V^pi (s) = max_pi sum_a pi(a mid(|) s) Q^pi (s, a)
  $

  接下来我们先证明其 $<= max_a Q^* (s, a)$。由定义，对任意 $pi$ 有：

  $
  sum_a pi(a mid(|) s) Q^pi (s, a) <= sum_a pi(a mid(|) s) max_(a') Q^* (s, a') = max_a Q^* (s, a)
  $

  然后证明其 $>= max_a Q^* (s, a)$。构造一个策略 $tilde(pi)$ 如下：在状态 $s$ 采取动作 $tilde(a) = arg max_a Q^* (s, a)$，即有 $Q^* (s, tilde(a)) = max_a Q^* (s, a)$，之后按照最优策略执行。该策略在状态 $s$ 的期望回报为：

  $
  V^tilde(pi) (s) = sum_a tilde(pi) (a mid(|) s) Q^tilde(pi) (s, a) = Q^* (s, tilde(a)) = max_a Q^* (s, a)
  $
  
  又由于 $V^* (s)$ 是所有策略中可以取到的最大期望回报，故有：

  $
  max_pi sum_a pi(a mid(|) s) Q^pi (s, a) = V^* (s) >= V^tilde(pi) (s) = max_a Q^* (s, a)
  $

  两个不等式联合得到 @equ:vq_opt_relation。
])

最优值函数和最优动作值函数也有迭代形式，即*最优贝尔曼方程*（optimal Bellman equations），证明见后：

$
V^* (s) &= max_a [r(s, a) + gamma EE_(s'~P(dot mid(|) s, a)) [V^* (s')]] \
&= max_a [r(s, a) + gamma sum_(s') P(s' mid(|) s, a) V^* (s')]
$ <equ:v_bellman_optimal>

$
Q^* (s, a) &= r(s, a) + gamma EE_(s'~P(dot mid(|) s, a)) [max_(a') Q^* (s', a')] \
&= r(s, a) + gamma sum_(s') P(s' mid(|) s, a) max_(a') Q^* (s', a')
$ <equ:q_bellman_optimal>

课件中#underline[考虑状态可微的情况]，则有最优 Q 值（optimal Q-values）：

$
Q^* (s, a) = r(s, a) + gamma integral P(s' mid(|) s, a) max_(a') Q^* (s', a') dif s'
$

对应的最优策略为在任意状态下都选取令最有动作值函数最大的动作，即：

$
pi^* (a mid(|) s) = 1, quad "iff" a = arg max_(a') Q^* (s, a')
$

#blockquote([
  *关于最优贝尔曼方程的证明*：
  
  从定义出发，利用前面证明的一些结论很容易得到最优动作值函数的版本：

  $
  Q^* (s, a) &= max_pi Q^pi (s, a) \
  &= max_pi [r(s, a) + gamma EE_(s'~P(dot mid(|) s, a)) [V^pi (s')]] \
  &=^* r(s, a) + gamma EE_(s'~P(dot mid(|) s, a)) [max_pi V^pi (s')] \
  &= r(s, a) + gamma EE_(s'~P(dot mid(|) s, a)) [V^* (s')] \
  &= r(s, a) + gamma EE_(s'~P(dot mid(|) s, a)) [max_(a') Q^* (s', a')]
  $

  其中需要注意的是 $*$ 式，这一步将 $max_pi$ 放进了期望中，不具备一般性，因为前者要求同一个策略在所有可能后续状态 $s'$ 上表现最佳，放进期望后则允许对每个状态 $s'$ 分别选取最佳策略。故#underline[欲使该步成立，需要满足前述最优策略存在的假设]，即 @equ:optimal_policy_exists_assumption。

  接下来代入 @equ:vq_opt_relation 即可直接得到最优值函数的版本：

  $
  V^* (s) &= max_a Q^* (s, a) \
  &= max_a [r(s, a) + gamma EE_(s'~P(dot mid(|) s, a)) [V^* (s')]]
  $
])

考虑最优策略价值函数 $J^* := max_pi J^pi$：

$
max_pi J^pi &= max_pi EE_rho [V^pi (s_0)] \
&=^* EE_rho [max_pi V^pi (s_0)] = EE_rho [V^* (s_0)] \
&= EE_rho [max_a Q^* (s_0, a)] \
&=^* EE_(s_0~rho(dot)) [max_a (r(s_0, a) + gamma EE_(s'~P(dot mid(|) s_0, a)) [max_(a') Q^* (s', a')])] \
&= EE_(s~rho(dot)) [max_a (r(s, a) + gamma EE_(s'~P(dot mid(|) s, a)) [max_(a') max_pi Q^pi (s', a')])]
$

其中两个 $*$ 式的含义仍表示需要 @equ:optimal_policy_exists_assumption 的最优策略存在；前者是为了将 $max_pi$ 放进期望中，后者是证明最优动作值函数贝尔曼方程时用到了该假设。最后一步进行了一些等效替换。

// == End-to-end Learning

// *端到端学习*（end-to-end learning）
