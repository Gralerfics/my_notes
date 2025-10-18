#import "../generic.typ": *

#import "@preview/cetz:0.4.2"
#import "@preview/fletcher:0.5.8" as fletcher: diagram, node, edge

= Stochastic Modelling Identification

对于随机过程的建模同对确定信号的建模主要存在两方面差别。第一，确定建模中由于已知 $x[n]$ 的具体样本值，故误差依赖于确定样本值定义，而随机建模中我们只掌握 $x[n]$ 的统计特征，不再适合使用先前的 $e[n]$ 定义。第二就是输入信号的差别，由于对随机过程进行建模，所以输入也不再适合使用单位脉冲信号，而是使用单位方差的白噪声（White Noise），见 @fig:signal_model_stochastic。

在考虑这些区别的前提下，我们还要对要建模的随机过程作平稳性假设，即假定随机过程是 WSS 的。

还是类似地，对于随机过程我们可以将 @sec:dmi_ls_method 中的最小二乘误差换成均方误差 $cal(E)_"MS" = E{abs(x[n]-hat(x)[n])^2}$ 来进行优化，但也会遇到同样的非线性问题，难以处理，需要寻找别的方案。

== Auto-regressive Moving Average Processes

我们先定义一种叫 ARMA 过程的随机过程。考虑使用 ARMA 模型的传递函数（@equ:signal_modelling_arma_tf）对方差为 $sigma_v^2$ 的白噪声 $v[n]$ 进行滤波，得到输出 $x[n]$。

#figure(
    caption: [Diagram of ARMA process generation]
)[
    #diagram(
        spacing: (10mm, 8mm),
        node-stroke: 0.8pt,
        edge-stroke: 0.8pt,
        node((0, 0), [$H(z) = B(z) / A(z) = (sum_(k=0)^q b[k] z^(-k)) / (1 + sum_(k=1)^p a[k] z^(-k))$], inset: 8pt),
        edge((-2, 0), (0, 0), "-|>", [White noise $v[n]$], label-pos: 0),
        edge((0, 0), (2, 0), "-|>", [$x[n]$], label-pos: 1),
    )
] <fig:smi_arma_proc_diagram>

这里*假设滤波器是稳定（Stable）的*，那么输出的随机过程 $x[n]$ 将会是 WSS 的#text(fill: red, "（TODO）")。由于白噪声的功率谱为 $P_v (z) = sigma_v^2$，得到 $x[n]$ 的功率谱：

$
P_x (z) = sigma_v^2 (B(z) B^*(1\/z^*))/(A(z) A^*(1\/z^*))
$

在频域即为：

$
P_x (e^(j omega)) = sigma_v^2 abs(B(e^(j omega)))^2 / abs(A(e^(j omega)))^2
$

我们定义用于这种形式的功率谱的过程为 $"ARMA"(p, q)$ 过程。可以注意到由于对称性，其功率谱有 $2p$ 和极点和 $2q$ 个零点#text(fill: red, "（TODO）")。

#blockquote[
    再次澄清，@sec:signal_modelling_arma_model 中提及 ARMA 模型，此处是在说 ARMA 过程。后者是指将白噪声放入 ARMA 模型后输出信号满足的随机过程。
]

=== Yule-Walker Equations <sec:signal_modelling_yule_walker>

在随机建模问题上，我们希望构建的模型输出具有与目标信号相同的统计特征，如考虑输出的自相关 $r_x(k)$ 同目标信号一致。所以接下来我们需要*建立起模型输出的自相关 $r_x (k)$ 与系统参数 $a[dot]$、$b[dot]$ 乃至单位冲激响应 $h[n]$ 之间的统计关系*。

由定义，对于来自 $v[n]$ 的 ARMA 过程 $x[n]$，满足如下方程：

$
x[n] + sum_(l=1)^p a[l] x[n-l] = sum_(l=0)^q b[l] v[n-l]
$

我们可以由此式推得 $x[n]$ 的自相关与 $x[n]$ 同 $v[n]$ 的互相关之间满足同样形式的关系，具体操作是在等式两侧同乘以 $x^*[n-k]$ 并取期望：

$
E{x[n] x^*[n-k] + sum_(l=1)^p a[l] x[n-l] x^*[n-k]} = E{sum_(l=0)^q b[l] v[n-l] x^*[n-k]}
$

即：

$
E{x[n] x^*[n-k]} + sum_(l=1)^p a[l] E{x[n-l] x^*[n-k]} = sum_(l=0)^q b[l] #text(fill: blue, $E{v[n-l] x^*[n-k]}$)
$

在平稳性假设成立的前提下，代入自相关和互相关的定义（见 @sec:fun_rp_statistic）得：

$
r_x (k) + sum_(l=1)^p a[l] r_x (k-l) = sum_(l=0)^(q) b[l] #text(fill: blue, $r_(v x) (k-l)$)
$

互相关项 $r_(v x) (k-l)$ 的存在让式子仍然包含 $v$，我们可以用表示系统属性的单位冲激响应 $h[n]$ 来换掉它，方法是继续代入 $x[n] = v[n] * h[n] = sum_(m=-infinity)^infinity v[m] h[n-m]$：

$
r_(v x) (k-l) &= E{v[k] x^*[l]} \
&= E{v[k] (sum_(m=-infinity)^infinity v[m] h[l-m])^*} \
&= sum_(m=-infinity)^infinity #text(fill: purple, $E{v[k] v^*[m]}$) h^*[l-m]
$

#blockquote[
    这里使用的 $r_(v x) (k-l) = E{v[k] x^*[l]}$ 由于平稳性假设同前面定义的 $r_(v x) (k-l) = E{v[n-l] x^*[n-k]}$ 是一致的，因为$(n-l)-(n-k)=k-l$。这样换一下推导更简洁。
]

由于 $v[n]$ 是独立同分布、方差为 $sigma_v^2$ 的白噪声，故有：

$
#text(fill: purple, $E{v[k] v^*[m]}$) = cases(
    sigma_v^2\, quad &m = k,
    0\, &"otherwise"
)
$

即该项在 $m != k$ 时都为 $0$，代入前式得到：

$
#text(fill: blue, $r_(v x) (k-l)$) = sigma_v^2 h^*[l-k]
$

于是我们得到了不包含 $v$ 的表达式：

$
r_x (k) + sum_(l=1)^p a[l] r_x (k-l) = sigma_v^2 sum_(l=0)^(q) b[l] h^*[l-k]
$

最后考虑现实情况修饰一下，*假设系统是因果（Causal）的*，即 $h[n]$ 在 $n<0$ 时取值皆为 $0$，那么 $h^*[l-k]$ 在 $l < k$ 时就为 $0$，可以修改等式右侧项所含求和的上下限，并记为 $c[k]$：

#emphasis_equbox([
$
c[k] := sum_(l=0)^(q) b[l] h^*[l-k] = sum_(l=k)^(q) b[l] h^*[l-k] = sum_(l=0)^(q-k) b[l+k] h^*[l]
$
])

顺便，$k>q$ 时该项为 $0$，最后我们得到 *Yule-Walker Equations*：

#emphasis_equbox([
$
r_x (k) + sum_(l=1)^p a[l] r_x (k-l) = cases(
    sigma_v^2 c[k]\, quad &0<=k<=q,
    0\, &k>q
)
$
])

这就是本节的目标，即模型输出的自相关 $r_x (k)$ 与系统参数 $a[dot]$、$b[dot]$ 和单位冲激响应 $h[n]$ 之间的统计关系，为之后求解模型参数铺垫。

=== Modified Yule-Walker Equation (MYWE) Method

由前，我们应用单位方差假设，即 $sigma_v^2 = 1$，后续不再专门提及。



== AM Processes

== MA Processes
