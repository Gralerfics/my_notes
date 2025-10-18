#import "../generic.typ": *

#import "@preview/cetz:0.4.2"
#import "@preview/fletcher:0.5.8" as fletcher: diagram, node, edge

= Stochastic Modelling Identification

对于随机过程的建模同对确定信号的建模主要存在两方面差别。第一，确定建模中由于已知 $x[n]$ 的具体样本值，故误差依赖于确定样本值定义，而随机建模中我们只掌握 $x[n]$ 的统计特征，不再适合使用先前的 $e[n]$ 定义。第二就是输入信号的差别，由于对随机过程进行建模，所以输入也不再适合使用单位脉冲信号，而是使用单位方差的白噪声（White Noise），见 @fig:signal_model_stochastic。

在考虑这些区别的前提下，我们还要对要建模的随机过程作平稳性假设，即假定随机过程是 WSS 的。

还是类似地，对于随机过程我们可以将 @sec:dmi_ls_method 中的最小二乘误差换成均方误差 $cal(E)_"MS" = E{abs(x[n]-hat(x)[n])^2}$ 来进行优化，但也会遇到同样的非线性问题，难以处理，需要寻找别的方案。

#blockquote[
    注意这里开始的 $x[n]$、$v[n]$ 等不再是具体的离散信号，而是代表一个随机过程。实际应用中我们可能并不直接知道其统计特征，此时才需要从具体的信号（实现、样本）中估计出统计特征。如果只有一条实现用于估算，则应满足相应的平稳性、遍历性假设才有意义。
    
    方便起见就不将后面所有方括号改成圆括号了。
]

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

在随机建模问题上，我们希望构建的模型输出具有与目标过程相同的统计特征，如考虑输出的自相关 $r_x(k)$ 同目标过程一致。所以接下来我们需要*建立起模型输出的自相关 $r_x (k)$ 与系统参数 $a[dot]$、$b[dot]$ 乃至单位冲激响应 $h[n]$ 之间的统计关系*。

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
c[k] :&= sum_(l=0)^(q) b[l] h^*[l-k] = sum_(l=k)^(q) b[l] h^*[l-k] = sum_(l=0)^(q-k) b[l+k] h^*[l] \ &= b[k] * h^*[-k]
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

注意，之后我们将默认应用单位方差假设，即取 $sigma_v^2 = 1$。
])

这就是本节的目标，即模型输出的自相关 $r_x (k)$ 与系统参数 $a[dot]$、$b[dot]$ 和单位冲激响应 $h[n]$ 之间的统计关系。

#blockquote[
    顺带一提，$k > q$ 时：

    $
    r_x (k) = - sum_(l=1)^p a[l] r_x (k-l)
    $

    可以用滤波器参数和已知的自相关函数值外推自相关函数之后的值。
]

清晰一点，我们也把它的矩阵形式画出来：

$
mat(
    delim: "[",
    column-gap: #1.0em,
    row-gap: #0.5em,
    augment: #(hline: 4, stroke: (dash: (2pt, 2pt))),
    r_x (0), r_x (-1), dots, r_x (-p);
    r_x (1), r_x (0), dots, r_x (-p+1);
    dots.v, dots.v, dots.down, dots.v;
    r_x (q), r_x (q-1), dots, r_x (q-p);
    r_x (q+1), r_x (q), dots, r_x (q-p+1);
    dots.v, dots.v, dots.down, dots.v;
    r_x (q+p), r_x (q+p-1), dots, r_x (q);
    dots.v, dots.v, dots.down, dots.v;
)
mat(
    delim: "[",
    row-gap: #0.5em,
    1;
    a[1];
    a[2];
    dots.v;
    a[p];
)
=
sigma_v^2
mat(
    delim: "[",
    row-gap: #0.5em,
    augment: #(hline: 4, stroke: (dash: (2pt, 2pt))),
    c[0];
    c[1];
    dots.v;
    c[q];
    0;
    dots.v;
    0;
    dots.v;
)
=
mat(
    delim: "[",
    row-gap: #0.5em,
    augment: #(hline: 4, stroke: (dash: (2pt, 2pt))),
    c[0];
    c[1];
    dots.v;
    c[q];
    0;
    dots.v;
    0;
    dots.v;
)
$

=== Modified Yule-Walker Equation (MYWE) Method

Yule-Walker 方程可以用于从自相关函数求解滤波器参数，但由于 $h^*[l]$ 的存在，它仍然是一个较难处理的非线性问题。

再次澄清，在该问题中我们对随机过程建模而不是对具体的信号建模，故视目标过程的统计特征（如自相关函数）是已知的，如果我们不知其自相关函数 $r_v (k)$ 的值，才需要通过统计方法从一些实现（样本）中估算得到 $hat(r)_v (k)$。

回到参数辨识的问题上来，我们可以先仿照 Pade 法的思路，通过分步求解来*近似*最优结果。先用 $q<k<=q+p$ 的部分估计 $a[dot]$，对应的式子为：

$
mat(
    delim: "[",
    r_x (q), r_x (q-1), dots, r_x (q-p+1);
    r_x (q+1), r_x (q), dots, r_x (q-p+2);
    dots.v, dots.v, dots.down, dots.v;
    r_x (q+p-1), r_x (q+p-2), dots, r_x (q);
)
mat(
    delim: "[",
    a[1];
    a[2];
    dots.v;
    a[p];
)
=
-
mat(
    delim: "[",
    r_x (q+1);
    r_x (q+2);
    dots.v;
    r_x (q+p);
)
$

该方程组称为 *Modified Yule-Walker equations*（MYWE），于是该方法称为 MYWE 法。值得注意的是该方程组形式与 Pade 法中 @equ:deterministic_model_identification_pade_a_bar 的形式完全一样，只是将自相关函数换成了 $x[n]$ 的值。该矩阵也是 Toeplitz 的，所以可以使用 Trench 算法等加速求解。

得到 $a[dot]$ 后，*第二步需要求解 $b[dot]$*。若将 $a[dot]$ 代回 Yule-Walker 方程我们可以得到 $c[dot]$ 的值。但 $c[k] := b[k] * h^*[-k]$，$h[k]$ 甚至还依赖 $b[k]$，想求出 $b[dot]$ 非常困难。课件称 "We skip this"，似乎不打算管这部分。参考书的相关内容大致从 190 页开始，主要提到几种方法。

第一，我们已知 $a[dot]$，用其构造一个 AR 滤波器 $A(z)$ 对 $x[n]$ 进行滤波可以得到新过程 $y[n]$：

$
P_x (z) = (B(z) B^*(1\/z^*))/(A(z) A^*(1\/z^*)) quad attach(-->, t: A(z)) quad P_y (z) = B(z) B^*(1\/z^*)
$

这个过程是一个 MA 过程，我们再用之后 @sec:smi_ma_processes 中的方法处理，估计 $b[dot]$。

第二，不显式地进行滤波，不过本质应该和第一种一样。通过 Yule-Walker 方程上半部分求出 $c[dot]$ 后，求正半轴的拉普拉斯变换得到（因为通过 Yule-Walker 方程只能求出其正半轴的值）：

$
[C(z)]_+ = sum_(k=0)^infinity c[k] z^(-k)
$

相应地，虽然我们不知道，但其负半轴的拉普拉斯变换为：

$
[C(z)]_- = sum_(k=-infinity)^(-1) c[k] z^(-k) = sum_(k=1)^infinity c[-k] z^k
$

由定义 $c[k] := b[k] * h^*[-k]$ 又得到 MA 过程的功率谱：

$
C(z) = B(z) H^*(1\/z^*) = B(z) (B^*(1\/z^*))/(A^*(1\/z^*)) \
quad => quad P_y (z) equiv C(z) A^*(1\/z^*) = B(z) B^*(1\/z^*)
$

我们将其拆开写：

$
P_y (z) = C(z) A^*(1\/z^*) = [C(z)]_+ A^*(1\/z^*) + [C(z)]_- A^*(1\/z^*)
$

由于 $a[k]$ 负半轴值为 $0$，则 $A^*(1\/z^*)$ 只包含 $z$ 的正功率，同时 $[C(z)]_+$ 也是如此#text(fill: red, "（TODO，这里书上写的减号？）")，故 $P_y (z)$ 的 causal part 即：

$
[P_y (z)]_+ = [C(z)]_+ A^*(1\/z^*)]_+
$

于是虽然我们不知道 $c[k]$ 的负半轴部分的值，但通过该式，我们可以用已知的 $c[dot]$ 正半轴的值伙同 $a[dot]$ 求出 $[P_y (z)]_+$，再由共轭对称性得到完整的 $P_y (z)$。最后对其进行谱分解（Spectral Factorization）得到系数 $b[dot]$：

$
P_y (z) = B(z) B^*(1\/z^*)
$

#text(fill: red, "（TODO，要不抄过来）")参考书 192 页还有一个清晰的例子。

=== Extended Yule-Walker Equation Method

相应地，在第一步中我们也可以使用类似 Prony 法的方式，将 $k>q$ 的所有式子都纳入考虑，得到的超定方程组称为 *Extended Yule-Walker equations*。

然后我们公式化求一个最小二乘解，过程也与 Prony 法相似，具体不再赘述。

== Auto-regressive Processes

我们又来考虑 all-pole 的情况，由于只剩下一个 $b[0]$，方程可以简化很多：

#emphasis_equbox([
$
r_x (k) + sum_(l=1)^p a[l] r_x (k-l) = abs(b[0])^2 delta(k), quad k>=0
$
])

由于不存在复杂的 $c[k]$ 的问题，可以直接沿用 Prony 法，先用除了第一个以外的式子求 $a[dot]$，可以注意到式子和 @equ:dmi_prony_all_pole_matrix 完全一样；然后再用第一个式子推得 $b[0]$。这被称为 Yule-Walker 法（什么混乱的起名方式）。

再次说明，如果我们不知道目标过程的自相关，就需要从（满足遍历性假设的）样本中去估计，例如：

$
hat(r)_x (k) = 1/N sum_(n=0)^(N-1) x[n] x[n-k]
$

这么一搞其实就和前面的 Auto-correlation 法（@sec:dmi_finite_data_autocorrelation_method）完全等价了，也符合直觉。

== Moving Average Processes <sec:smi_ma_processes>

对于 MA 过程，我们代入 Yule-Walker 方程后得到：

$
r_x (k) = sum_(l=0)^q b[l] b^*[l-k] = b[k] * b^*[-k]
$

#text(fill: red, "（TODO，书 195 页具体说明）")

即有：

$
P_x (z) = B(z) B^*(1\/z^*)
$

总结来说，就是将自相关函数 $r_x (k)$ 作 z 变换后得到功率谱 $P_x (z)$，然后进行谱分解即可得到结果，举个例子：

$
r_x (k) = 17 delta(k) + 4[delta(k-1) + delta(k+1)]
$

z 变换得到：

$
P_x (z) = 17 + 4z^(-1) + 4z = (4 + z^(-1))(4+z)
$

于是有：

$
B(z) = 4 + z^(-1) quad "or" quad B(z) = 1 + 4z^(-1)
$

此外，还有 Durbin's method 等方法，此处不记录。
