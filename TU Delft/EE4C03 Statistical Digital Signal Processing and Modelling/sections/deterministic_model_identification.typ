#import "../generic.typ": *

#import "@preview/cetz:0.4.2"
#import "@preview/fletcher:0.5.8" as fletcher: diagram, node, edge

= Deterministic Modelling Identification

模型建立后我们还需要对其进行参数辨识，即确定参数序列 $a[k]$ 和 $b[k]$ 的值。选取参数的目标是使模型输出 $hat(x)[n]$ 尽可能接近目标信号 $x[n]$。

== Least Squares (LS) Method <sec:dmi_ls_method>

首先讨论 Deterministic Modelling，我们希望模型输出 $hat(x)[n]$ 能够精确地重现目标信号 $x[n]$，即每个采样点的值都要尽可能接近。定义误差信号 $e'[n] = x[n] - hat(x)[n]$，可以通过最小化均方误差 $cal(E)_"LS" = sum_(n=0)^infinity abs(e'[n])^2$ 来确定模型参数，如 @fig:deterministic_model_identification_diagram_ls 所示。

#figure(
    caption: [Diagram for deterministic model identification (intractable)]
)[
    #diagram(
        spacing: (10mm, 8mm),
        node-stroke: 0.8pt,
        edge-stroke: 0.8pt,
        node((2, 0), [$plus$], inset: 4pt),
        node((0, 0), [$H(z) = frac(B(z), A(z))$], inset: 8pt),
        edge((-2, -1), "rrrr,d", "-|>", [$x[n]$], label-pos: 0),
        edge((-2, 0), "rr", "-|>", [$delta[n]$], label-pos: 0),
        edge((0, 0), "rr", "-|>", [$h[n]$], label-pos: 0.35),
        edge((0, 0), "rr", "-|>", [$-$], label-pos: 0.94, label-side: right),
        edge((2, 0), "rr", "-|>", [$e'[n] = x[n] - h[n]$], label-pos: 1.2),
    )
] <fig:deterministic_model_identification_diagram_ls>

该优化问题可通过求解如下方程组解决（对变量的共轭求偏导的原因详见 @sec:fun_optimization）：

$
cases(
    frac(partial cal(E)_"LS", partial a^*[k]) = 0\, quad quad & k = 1\, 2\, ...\, p,
    frac(partial cal(E)_"LS", partial b^*[k]) = 0\, & k = 0\, 1\, ...\, q
)
$

但这个方程组是非线性的，求解起来非常复杂。

// #blockquote[
//     具体地，由 Parseval 定理可知：

//     $ cal(E)_"LS" = 1/(2 pi) integral_(-pi)^pi abs(E'(e^(j omega)))^2 dif omega $

//     #text(fill: red, "（TODO）")
// ]

== Padé Approximation

注意到我们使用的 ARMA 模型具有 $p+q+1$ 个参数，也就是说它拥有 $p+q+1$ 个自由度，所以理论上我们是可以用其完美拟合信号的前 $p+q+1$ 个样本的，我们先考虑这个任务。

我们对传递函数的形式作一个转化：

$ H(z) = B(z) / A(z) quad => quad H(z) A(z) = B(z) quad => quad h[n] * a[n] = b[n] $

将卷积展开得到：

$ h[n] + sum_(k=1)^p a[k] h[n-k] = b[n] $

对于单位脉冲输入 $delta[n]$，系统的输出 $h[n]$ 就是我们估计的结果 $hat(x)[n]$。若要完全拟合前 $p+q+1$ 个样本，则直接代入 $h[n] = x[n], 0<=n<=p+q$，得到：

$
x[n] + sum_(k=1)^p a[k] x[n-k] = cases(
    b[n]\, quad &n=0\, 1\, dots\, q,
    0\, &n=q+1\, q+2\, dots\, q+p
)
$ <equ:deterministic_model_identification_equation_pade>

这就是我们想要的线性方程组了，包含了 $a[dot]$、$b[dot]$ 和已知的常数 $x[n]$。

#blockquote[
    上述转化过程体现在系统框图上就是在两路上同乘 $H(z)$ 的分母 $A(z)$。即令新的目标误差 $E(z) = A(z) E'(z) = A(z) X(z) - B(z)$，如 @fig:deterministic_model_identification_diagram_pade 所示。

    #figure(
        caption: [Diagram for deterministic model identification]
    )[
        #diagram(
            spacing: (10mm, 8mm),
            node-stroke: 0.8pt,
            edge-stroke: 0.8pt,
            node((2, 0), [$plus$], inset: 4pt),
            node((0, 0), [$B(z)$], inset: 8pt),
            node((0, -1), [$A(z)$], inset: 8pt),
            edge((-2, -1), "rr", "-|>", [$x[n]$], label-pos: 0),
            edge((0, -1), "rr,d", "-|>", [$hat(b)[n]$], label-pos: 0.155),
            edge((-2, 0), "rr", "-|>", [$delta[n]$], label-pos: 0),
            edge((0, 0), "rr", "-|>", [$b[n]$], label-pos: 0.35),
            edge((0, 0), "rr", "-|>", [$-$], label-pos: 0.94, label-side: right),
            edge((2, 0), "rr", "-|>", [$e[n] = hat(b)[n] - b[n]$], label-pos: 1.2),
        )
    ] <fig:deterministic_model_identification_diagram_pade>

    单位脉冲信号经过 $B(z)$ 得到的输出就是 $b[n]$，将其与用目标信号 $x[n]$ 经过 $A(z)$ 得到的 $hat(b)[n]$（此时可视为是对系数序列 $b[n]$ 的估计）作差，得到新的误差。
    
    经过这样的操作后的列出的方程就是线性的了。
]

接下来是对 @equ:deterministic_model_identification_equation_pade 的求解。该方程组包含相同个数的方程和未知数，非奇异的话恰好可以求解出唯一解。清晰一点，我们把矩阵画出来：

#let pade_matrix = [$
mat(
    delim: "[",
    column-gap: #0.7em,
    row-gap: #0.6em,
    x[0], 0, 0, dots, 0;
    x[1], x[0], 0, dots, 0;
    x[2], x[1], x[0], dots, 0;
    dots.v, dots.v, dots.v, dots.down, dots.v;
    x[p], x[p-1], x[p-2], dots, x[0];
    dots.v, dots.v, dots.v, dots.down, dots.v;
    x[q], x[q-1], x[q-2], dots, x[q-p];
    x[q+1], x[q], x[q-1], dots, x[q-p+1];
    x[q+2], x[q+1], x[q], dots, x[q-p+2];
    dots.v, dots.v, dots.v, dots.down, dots.v;
    x[q+p], x[q+p-1], x[q+p-2], dots, x[q];
)
mat(
    delim: "[",
    row-gap: #0.6em,
    1;
    a[1];
    a[2];
    dots.v;
    a[p];
)
=
mat(
    delim: "[",
    row-gap: #0.6em,
    b[0];
    b[1];
    b[2];
    dots.v;
    b[p];
    dots.v;
    b[q];
    0;
    0;
    dots.v;
    0;
)
$ <equ:deterministic_model_identification_matrix_pade>]

#align(center)[
    #cetz.canvas({
        import cetz.draw: *

        content(
            (0mm, 0mm),
            [#pade_matrix]
        )

        set-style(
            stroke: 0.8pt,
            mark: (transform-shape: false)
        )

        let sub_box = (x, y, w, h, color, label, lx, ly, lc) => {
            set-style(stroke: (paint: color, thickness: 0.9pt, dash: (3pt, 2pt)))
            rect((x, y), (x + w, y + h))
            if (label != none) {
                content(
                    (x + w / 2 + lx, y + h / 2 + ly),
                    text(fill: lc, label)
                )
            }
        }

        sub_box(-6.25, -0.97, 9.485, 4.65, blue, $bold(X)_0$, 0, 2.6, blue)
        sub_box(-6.25, -1.06, 1.55, -2.61, purple, $bold(x)_(q+1)$, 0, -1.68, purple)
        sub_box(-4.6, -1.06, 7.835, -2.61, purple, $bold(X)_q$, 0, -1.68, purple)
        
        sub_box(3.7, -1.64, 0.79, 2.62, green, $dash(bold(a))$, 0, -1.6, green)
        
        sub_box(5.5, -0.97, 0.78, 4.65, gray, $bold(b)$, 0, 2.6, gray)
    })
]

在这里，我们先用下半部分（后 $p$ 行）求解 $dash(bold(a))$（即 $a[dot]$）：

$
mat(delim: "[", bold(x)_(q+1), bold(X)_q) bold(a) = bold(0)
quad &<=> quad
mat(delim: "[", bold(x)_(q+1), bold(X)_q) mat(delim: "[", 1; dash(bold(a))) = bold(0) \
quad => quad
bold(X)_q dash(bold(a)) = - bold(x)_(q+1)
quad &=> quad
dash(bold(a)) = - bold(X)^(-1)_q bold(x)_(q+1)
$ <equ:deterministic_model_identification_pade_a_bar>

注意到 $bold(X)_q$ 是一个非对称 Toeplitz 矩阵，存在一些更高效的专用方法用于求解其逆矩阵，如 Trench 算法。接下来，代入上半部分（前 $q+1$ 行）即可得 $b[dot]$：

$ bold(b) = bold(X)_0 mat(delim: "[", 1; dash(bold(a))) $ <equ:dmi_b_X0a>

Padé 法很直接，但显然也存在一系列问题：
+ 不保证所得系统是稳定的；
+ 只约束了模型输出 $hat(x)[n]$ 和目标信号 $x[n]$ 的前 $p + q + 1$ 个样本相同，此后的匹配效果可能不佳；
+ $bold(X)_q$ 可能奇异且无解。

#blockquote[
    关于 $bold(X)_q$ 奇异且无解的情况，可以认为是模型中 $a[0] = 1$ 的默认假设存在问题。如果修改模型，令 $a[0] = 0$，则方程虽奇异但非无解，而是解不唯一。

    如果不是全极点模型，则如此得到的新传递函数在分子和分母中同时具有因子 $z$，可以约去。这本质上是在零处发生了零极点对消，结果上看等效于模型的阶数下降了，即模型阶数存在冗余。
]

== Prony's Method

Padé 法将所有的自由度都用在了序列的前 $p+q+1$ 项上，而 Prony 法的想法很简单，就是降低前面这段序列的拟合要求，从而获得一个从信号整体上看更好的拟合。

=== Prony Normal Equations <equ:dmi_prony_normal_equations>

我们先从比较形式化的角度来推导 Prony 法。具体地，还是借 Padé 法的想法将问题转化为线性的，如 @fig:deterministic_model_identification_diagram_pade 所示。写出整个信号误差的表达式，而不只是前 $p+q+1$ 项：

$
e[n] = cases(
    x[n] + sum_(k=1)^p a[k] x[n-k] - b[n]\, quad &n=0\, 1\, dots\, q,
    x[n] + sum_(k=1)^p a[k] x[n-k]\, &n>q
)
$ <equ:deterministic_model_identification_error_prony>

Prony 法中，我们先通过最小化均方误差的方式求解 $a[dot]$：

$
epsilon_(p, q) = sum_(n=q+1)^infinity abs(e[n])^2 = sum_(n=q+1)^infinity abs(x[n] + sum_(k=1)^p a[k] x[n-k])^2
$ <equ:dmi_epsilon_pq>

这里只考虑 $n>q$ 部分的误差是为了令该部分只与 $a[dot]$ 有关，这是考虑到分步求解 $a[dot]$ 和 $b[dot]$ 的需要，确实可能会牺牲一部分定义的准确性，但相对无限长的 $x[n]$ 影响并不大。接下来我们公式化求偏导令其为零计算最优值：

$
(partial epsilon_(p, q))/(partial a^*[k]) = sum_(n=q+1)^infinity (partial [e[n] e^*[n]])/(partial a^*[k]) = sum_(n=q+1)^infinity e[n] (partial e^*[n])/(partial a^*[k]) = 0, quad k=1, 2, dots, p
$ <equ:dmi_prony_normal_equ_derivation_startpoint>

由定义 @equ:deterministic_model_identification_error_prony 我们知道 $(partial e^*[n])/(partial a^*[k]) = x^*[n-k]$，代入得：

$
sum_(n=q+1)^infinity e[n] x^*[n-k] = 0, quad k=1, 2, dots, p
$ <equ:dmi_prony_orthogonality_principle>

这个式子表达了最小误差和信号间的正交关系，称为 Orthogonality principle。我们继续代入定义 @equ:deterministic_model_identification_error_prony （注意字母 $k$ 用掉了，换用 $l$，不要混淆）得：

$
sum_(n=q+1)^infinity (x[n] + sum_(l=1)^p a[l] x[n-l]) x^*[n-k] = 0, quad k=1, 2, dots, p
$

移项并重排求和符号顺序可以得到：

$
sum_(l=1)^p a[l] 

(sum_(n=q+1)^infinity x^*[n-k] x[n-l]) = -sum_(n=q+1)^infinity x^*[n-k] x[n], quad k=1, 2, dots, p
$

为了简化表达，我们记：

#emphasis_equbox([
$
r_x (k, l) := sum_(n=q+1)^infinity x^*[n-k] x[n-l]
$ <equ:deterministic_model_identification_prony_ne_rx>
])

我们可以顺便观察到 $r_x (k, l) = r_x^* (l, k)$。将其代入原式得到：

#emphasis_equbox([
$
sum_(l=1)^p a[l] r_x (k, l) = -r_x (k, 0), quad k = 1, 2, dots, p
$ <equ:dmi_prony_normal_equ>
])

这被称为 *Prony normal equations*，写成矩阵形式是：

#emphasis_equbox([
$
mat(
    delim: "[",
    r_x (1, 1), r_x (1, 2), dots, r_x (1, p);
    r_x (2, 1), r_x (2, 2), dots, r_x (2, p);
    dots.v, dots.v, dots.down, dots.v;
    r_x (p, 1), r_x (p, 2), dots, r_x (p, p);
)
mat(
    delim: "[",
    a[1];
    a[2];
    dots.v;
    a[p];
)
=
-mat(
    delim: "[",
    r_x (1, 0);
    r_x (2, 0);
    dots.v;
    r_x (p, 0);
)
$
])

记为：

#emphasis_equbox([
$
bold(R)_x dash(bold(a)) = - bold(r)_x
$ <equ:dmi_prony_normal_equ_matrix>
])

可以发现 $bold(R)_x$ 是一个 Hermitian 矩阵。求得 $a[dot]$ 后我们就可以代回 @equ:deterministic_model_identification_error_prony，令 $n=1, 2, dots, q$ 时误差为 $0$，得到 $b[k]$。

#blockquote[
    但我们需要注意一件事，原本最小化 $e[n]$ 的问题仍然是一个联合非线性最小二乘问题，$a[dot]$ 和 $b[dot]$ 是耦合的。所以我们分两步依此求解 $a[dot]$ 和 $b[dot]$ 的方法实际上是一种简化，并*不能保证全局最小化原始误差*。
]

// #blockquote[
//     *注意*，虽然 $r_x (k, l)$ 的定义非常像信号的样本自相关函数，但在这里本质上只是记号上的代换，不要弄错主次。

//     我们不是在估计自相关函数，也没有对产生信号的随机过程作任何平稳性、遍历性的假设。我们的逻辑是这样的：我们用 Prony 法自然地导出了等式，其恰好形式上可以用自相关矩阵来代换，而不是原理上要求我们需要用自相关来计算。所以我们不去纠结其公式形式上同 @equ:fun_rp_autocor_matrix 的细节差异，没有太大的意义。
    
//     我们只需要认识到该式从定义上确实具有和自相关类似的物理意义，即信号延迟后与自身的相似程度，主要目的是拓展我们的理解。

//     #text(fill: red, "（TODO）") // 如果信号真平稳，说明这个能估计信号背后的分布；如果信号自己不争气，则只能拟合这个具体的信号。
// ]

=== An Equivalent Perspective from Pseudoinverse <sec:dmi_equivalent_perspective_from_pseudoinverse>

上节的推导中我们自然地应用最小二乘法将问题视为优化问题处理，实际上我们也可以直接令所有的 $hat(x)[n]=x[n]$，得到一个超定方程组：

#let prony_matrix = [$
mat(
    delim: "[",
    column-gap: #0.7em,
    row-gap: #0.6em,
    x[0], 0, 0, dots, 0;
    dots.v, dots.v, dots.v, dots.down, dots.v;
    x[q], x[q-1], x[q-2], dots, x[q-p];
    x[q+1], x[q], x[q-1], dots, x[q-p+1];
    dots.v, dots.v, dots.v, dots.down, dots.v;
    x[q+p], x[q+p-1], x[q+p-2], dots, x[q];
    dots.v, dots.v, dots.v, dots.down, dots.v;
)
mat(
    delim: "[",
    row-gap: #0.6em,
    1;
    a[1];
    a[2];
    dots.v;
    a[p];
)
=
mat(
    delim: "[",
    row-gap: #0.6em,
    b[0];
    dots.v;
    b[q];
    0;
    dots.v;
    0;
    dots.v;
)
$ <equ:deterministic_model_identification_matrix_prony>]

#align(center)[
    #cetz.canvas({
        import cetz.draw: *

        content(
            (0mm, 0mm),
            [#prony_matrix]
        )

        set-style(
            stroke: 0.8pt,
            mark: (transform-shape: false)
        )

        let sub_box = (x, y, w, h, color, label, lx, ly, lc) => {
            set-style(stroke: (paint: color, thickness: 0.9pt, dash: (3pt, 2pt)))
            rect((x, y), (x + w, y + h))
            if (label != none) {
                content(
                    (x + w / 2 + lx, y + h / 2 + ly),
                    text(fill: lc, label)
                )
            }
        }

        sub_box(-6.25, 0.4, 9.485, 1.95, blue, $bold(X)_0$, 0, 1.25, blue)
        sub_box(-6.25, 0.3, 1.55, -2.65, purple, $bold(x)_(q+1)$, 0, -1.7, purple)
        sub_box(-4.6, 0.3, 7.835, -2.65, purple, $bold(X)_q$, 0, -1.7, purple)
        
        sub_box(3.7, -1.64, 0.79, 2.62, green, $dash(bold(a))$, 0, -1.6, green)
    })
]

接下来，我们可以将原本最小化均方误差的过程 "内化" 到使用伪逆求该方程组最小二乘解的过程当中。两种思路本质上是等效的，前者逻辑更顺畅，而后者有助于从线性空间的角度去理解问题。使用伪逆求解得到：

$ dash(bold(a)) = - text(fill: #purple, bold(X)^+_q) bold(x)_(q+1) = - text(fill: #purple, (bold(X)_q^H bold(X)_q)^(-1) bold(X)_q^H) bold(x)_(q+1) $

即最优的系数 $dash(bold(a))$ 将是如下方程组的解：

$ (bold(X)_q^H bold(X)_q) dash(bold(a)) = - bold(X)_q^H bold(x)_(q+1) $

作如下代换后我们再次得到 @equ:dmi_prony_normal_equ_matrix 的 Prony normal equations：

$ bold(R)_x = bold(X)_q^H bold(X)_q, quad bold(r)_x = bold(X)_q^H bold(x)_(q+1) $

可以计算验证 $bold(R)_x$ 与上节中的定义是一致的：

$
bold(R)_x
&=
bold(X)_q^H bold(X)_q \
&=
inline(
    mat(
        delim: "[",
        x^*[q], x^*[q+1], x^*[q+2], dots;
        x^*[q-1], x^*[q], x^*[q+1], dots;
        dots.v, dots.down, dots.v, dots.down;
        x^*[q-p+1], x^*[q-p+2], x^*[q-p+3], dots;
    )
    mat(
        delim: "[",
        x[q], x[q-1], dots, x[q-p+1];
        x[q+1], x[q], dots, x[q-p+2];
        x[q+2], x[q+1], dots, x[q-p+3];
        dots.v, dots.v, dots.down, dots.v;
    )
) \
&=
display(sum_(n=q+1)^infinity) inline(mat(
    delim: "[",
    column-gap: #1.0em,
    row-gap: #0.5em,
    x^*[n-1] x[n-1], x^*[n-1] x[n-2], dots, x^*[n-1] x[n-p];
    x^*[n-2] x[n-1], x^*[n-2] x[n-2], dots, x^*[n-2] x[n-p];
    dots.v, dots.v, dots.down, dots.v;
    x^*[n-p] x[n-1], x^*[n-p] x[n-2], dots, x^*[n-p] x[n-p];
)) \
&=
inline(
    mat(
        delim: "[",
        column-gap: #1.0em,
        row-gap: #0.5em,
        r_x (1, 1), r_x (1, 2), dots, r_x (1, p);
        r_x (2, 1), r_x (2, 2), dots, r_x (2, p);
        dots.v, dots.v, dots.down, dots.v;
        r_x (p, 1), r_x (p, 2), dots, r_x (p, p);
    )
)
$

该角度还提供了其他有效信息：关于 $A^H A$ 这种形式的矩阵，对任意向量 $bold(a)$ 有 $bold(a)^H (A^H A) bold(a) = (bold(A a))^H (bold(A a)) = norm(bold(A a))^2 >= 0$，这说明 $A^H A$ 是半正定矩阵。

由此，前述 Hermitian 矩阵 $bold(R)_x = bold(X)_q^H bold(X)_q$ 同时也是（半）正定矩阵，而该性质将决定 $A(z)$ 是（临界）稳定的，由此弥补了 Padé 法的一项缺陷#text(fill: red, "（TODO，是吗？是否还需要是 Toeplitz？）")。

此外，若 $bold(R)_x$ 为正定矩阵，则其特征值都是正数，即行列式不为零，矩阵可逆，解存在；若 $bold(R)_x$ 为包含零特征值的半正定矩阵，则奇异，但这实际上说明模型的阶数冗余，可以降低一些再尝试。

=== The Minimum Error and Augmented Normal Equations

由于我们求的是最小二乘解，所以最终拟合信号和真实信号间还是会存在一个最小误差，由 $e[n]$ 定义（@equ:deterministic_model_identification_error_prony）和 $epsilon_(p, q)$ 定义（@equ:dmi_epsilon_pq）继续推导：

$
epsilon_(p, q)
=
sum_(n=q+1)^infinity abs(e[n])^2
&=
sum_(n=q+1)^infinity e[n] (x[n] + sum_(k=1)^p a[k] x[n-k])^* \
&=
sum_(n=q+1)^infinity e[n] x^*[n] + sum_(n=q+1)^infinity e[n] (sum_(k=1)^p a[k] x[n-k])^* \
&=
sum_(n=q+1)^infinity e[n] x^*[n] + sum_(k=1)^p a^*[k] (sum_(n=q+1)^infinity e[n] x^*[n-k])
$

代入解最优时成立的 Orthogonality principle（ @equ:dmi_prony_orthogonality_principle）和 $e[n]$ 定义（@equ:deterministic_model_identification_error_prony）得：

$
epsilon_(p, q)
=
sum_(n=q+1)^infinity e[n] x^*[n]
&=
sum_(n=q+1)^infinity (x[n] + sum_(k=1)^p a[k] x[n-k]) x^*[n] \
&=
(sum_(n=q+1)^infinity x[n] x^*[n]) + sum_(k=1)^p a[k] (sum_(n=q+1)^infinity x[n-k] x^*[n])
$

用自相关序列 $r_x (k, l)$（@equ:deterministic_model_identification_prony_ne_rx）可简化为：

$ epsilon_(p, q) = r_x (0, 0) + sum_(k=1)^p a[k] r_x (0, k) $

化成这种形式后我们可以将 $epsilon_(p, q)$ 统一到方程 $bold(R)_x dash(bold(a)) = - bold(r)_x$ 中去（把常量移到矩阵最左侧一列了）：

#emphasis_equbox([
$
mat(
    delim: "[",
    column-gap: #1.0em,
    row-gap: #0.5em,
    augment: #(hline: 1, vline: 1, stroke: (dash: (2pt, 2pt))),
    text(fill: #blue, r_x (0, 0)), text(fill: #blue, r_x (0, 1)), text(fill: #blue, r_x (0, 2)), dots, text(fill: #blue, r_x (0, p));
    text(fill: #red, r_x (1, 0)), r_x (1, 1), r_x (1, 2), dots, r_x (1, p);
    text(fill: #red, r_x (2, 0)), r_x (2, 1), r_x (2, 2), dots, r_x (2, p);
    dots.v, dots.v, dots.v, dots.down, dots.v;
    text(fill: #red, r_x (p, 0)), r_x (p, 1), r_x (p, 2), dots, r_x (p, p);
)
mat(
    delim: "[",
    row-gap: #0.5em,
    augment: #(hline: 1, stroke: (dash: (2pt, 2pt))),
    text(fill: #red, 1);
    a[1];
    a[2];
    dots.v;
    a[p];
)
=
mat(
    delim: "[",
    row-gap: #0.5em,
    augment: #(hline: 1, stroke: (dash: (2pt, 2pt))),
    text(fill: #blue, epsilon_(p, q));
    0;
    0;
    dots.v;
    0;
)
$
])

或（$bold(u)_1$ 为首元素为 $1$ 其他为 $0$ 的单位向量）：

#emphasis_equbox([
$
bold(R)_x bold(a) = epsilon_(p, q) bold(u)_1
$
])

这样的形式称为 *Augmented normal equations*。

#blockquote[
    我们也可以用矩阵形式推导，更简洁一些。误差序列构成的向量为：

    $ bold(e) = bold(X)_q dash(bold(a)) + bold(x)_(q+1) $

    由最小二乘解的性质，我们取最优解即这个误差被最小化时，它一定是与 $bold(X)_q$ 的所有列正交的，即有：

    $ bold(X)_q^H bold(e) = 0 quad <=> quad bold(X)_q^H (bold(X)_q dash(bold(a)) + bold(x)_(q+1)) = 0 quad <=> quad bold(R)_x dash(bold(a)) = - bold(r)_x $

    对于阶数为 $p, q$ 的滤波器我们关心的最小误差为：

    $ epsilon_(p, q) = norm(bold(e))^2 = bold(e)^H bold(e) = (bold(X)_q dash(bold(a)) + bold(x)_(q+1))^H bold(e) = bold(x)^H_(q+1) bold(e) $

    最后一步是由于 $bold(X)_q^H bold(e) = 0$。继续代入 $bold(e)$ 最终得到同样的式子：

    $
    epsilon_(p, q) &= bold(x)^H_(q+1) bold(e) = bold(x)^H_(q+1) (bold(X)_q dash(bold(a)) + bold(x)_(q+1)) \
    &= bold(x)^H_(q+1) bold(x)_(q+1) + (bold(x)^H_(q+1) bold(X)_q) dash(bold(a)) \
    &= r_x (0, 0) + [r_x (0, 1), r_x (0, 2), dots, r_x (0, p)] dash(bold(a)) \
    &= r_x (0, 0) + sum_(k=1)^p a[k] r_x (0, k)
    $
]

== Special Case: All-pole Modelling

// 回忆 @equ:deterministic_model_identification_matrix_prony，我们通过 $bold(X)_q dash(bold(a)) = - bold(x)_(q+1)$ 求解系数，而在 All-pole 模型中 $q = 0$，我们代入得到一个看起来很干净的方程组：

// $ H(z) = Y(z) / X(z) = b[0] / (1 + sum_(k=1)^p a[k] z^(-k)) $

// $
// mat(
//     delim: "[",
//     x[0], 0, 0, dots, 0;
//     x[1], x[0], 0, dots, 0;
//     x[2], x[1], x[0], dots, 0;
//     dots.v, dots.v, dots.v, dots.down, dots.v;
//     x[p-1], x[p-2], x[p-3], dots, x[0];
//     x[p-2], x[p-3], x[p-4], dots, x[1];
//     dots.v, dots.v, dots.v, dots.down, dots.v;
// )
// mat(
//     delim: "[",
//     a[1];
//     a[2];
//     a[3];
//     dots.v;
//     a[p];
// )
// =
// -mat(
//     delim: "[",
//     x[1];
//     x[2];
//     x[3];
//     dots.v;
//     x[p];
//     x[p+1];
//     dots.v;
// )
// $

我们来研究全极点模型（@equ:ar_model_transfer_function），它在一些物理过程中很常见。

=== All-pole Normal Equations

首先，我们照例将求解 $a[dot]$ 的过程视作优化问题。参考 @equ:dmi_epsilon_pq 的误差定义，$q = 0$ 时我们有：

$
epsilon_(p, 0) = sum_(n=1)^infinity abs(e[n])^2
$

其中 $e[n]$ 定义仍来自 @equ:deterministic_model_identification_error_prony，但由于 $n=0$ 时 $n-k<0$，故 $a[dot]$ 的系数 $x[n-k]$ 值为 $0$，只剩 $x[0] - b[0]$：

$
e[n] = cases(
    x[0] - b[0]\, quad &n=0,
    x[n] + sum_(k=1)^p a[k] x[n-k]\, &n>0
)
$

这时候我们要作个妖，注意 $e[0] = x[0] - b[0]$ 对于 $a[dot]$ 可视作常数，故求解 $a[dot]$ 时最小化 $epsilon_(p, 0)$ 和最小化我们定义的一个新误差 $epsilon_(p)$ 是等价的：

$
epsilon_(p) = sum_(n=0)^infinity abs(e[n])^2
$

区别就是把 $e[0]$ 也放进去了。幸运的是，这样更换误差函数之后，我们依然会导出 Prony normal equations 的形式（见 @equ:dmi_prony_normal_equ），*但其中 $r_x (k, l)$ 的定义变化*为：

$
r_x (k, l) := sum_(n=0)^infinity x^*[n-k] x[n-l]
$ <equ:dmi_allpole_rxkl_new>

区别是求和下限从 $q + 1 = 1$ 换成了 $0$，*这是我们更换误差函数的影响*，具体的原因需要从 @equ:dmi_prony_normal_equ_derivation_startpoint 处开始用新的误差定义重新推导。

#blockquote[
    好吧其实很简单，我们换 $epsilon_(p)$ 后令其对 $a[dot]$ 偏导为零：

    $
    (partial #text(fill: red, $epsilon_(p)$))/(partial a^*[k]) = sum_(n=#text(fill: red, $0$))^infinity (partial [e[n] e^*[n]])/(partial a^*[k]) = sum_(n=#text(fill: red, $0$))^infinity e[n] (partial e^*[n])/(partial a^*[k]) = 0, quad k=1, 2, dots, p
    $

    由上面的 $e[n]$ 定义，我们还是有 $(partial e^*[n])/(partial a^*[k]) = x^*[n-k]$，因为 $n=0$ 时 $x^*[n-k]=0$，恰好和 $(partial (x[0] - b[0]))/(partial a^*[k]) = 0$ 一致，所以可以统一代入得：

    $
    sum_(n=#text(fill: red, $0$))^infinity e[n] x^*[n-k] = 0, quad k=1, 2, dots, p
    $

    继续代入 $e[n]$ 定义，依旧由于 $x^*[n-k]$ 项的存在，$n=0$ 的情况可以合并进来，直接得：

    $
    sum_(n=#text(fill: red, $0$))^infinity (x[n] + sum_(l=1)^p a[l] x[n-l]) #text(fill: blue, $x^*[n-k]$) = 0, quad k=1, 2, dots, p
    $

    移项并重排求和符号顺序可以得到：

    $
    sum_(l=1)^p a[l] 

    (sum_(n=#text(fill: red, $0$))^infinity x^*[n-k] x[n-l]) = -sum_(n=q+1)^infinity x^*[n-k] x[n], quad k=dots
    $

    该式同 Prony normal equations 的形式一致，区别只是 $r_x (k, l)$ 的定义需要如 @equ:dmi_allpole_rxkl_new 中所示，变为从 $0$ 开始求和。
]

到这里我们已经得到了全极点模型情况下求解 $a[dot]$ 的 "normal equations"，但观察一下还可以发现，由于 $x[n]$ 在 $n < 0$ 时值皆为 $0$，代入 @equ:dmi_allpole_rxkl_new 有：

$
r_x (k+1, l+1) &= sum_(n=0)^infinity x^*[n-(k+1)] x[n-(l+1)] \
&= sum_(n=0)^infinity x^*[n-1-k] x[n-1-l] \
&= sum_(n=-1)^infinity x^*[n-k] x[n-l] \
&= x^*[-1-k] x[-1-l] + sum_(n=0)^infinity x^*[n-k] x[n-l] \
&= r_x (k, l), quad (forall k, l>=0)
$

由此，我们可以令：

$
r_x (k-l) := r_x (k, l) = sum_(n=0)^infinity x^*[n-k] x[n-l]
$

即：

#emphasis_equbox([
$
r_x (k) = sum_(n=0)^infinity x^*[n-k] x[n]
$
])

// #blockquote[
//     *再次注意*，虽然这个东西很像平稳性的定义，矩阵、向量的记号也沿用了自相关定义的符号，但仍要强调这只是形式上恰好。
    
//     具体地，我们发现 $bold(R)_x$ 有 $x[n]$ 的自相关矩阵的含义，而在 All-pole 模型下发现这个矩阵变成 Toeplitz 的了，我们宁愿说这个矩阵不再代表信号的自相关矩阵，也不能说 All-pole 模型对信号的性质施加约束导致信号分布变成平稳的了（模型怎么会影响信号呢）。

//     但这确实说明了一件事，即虽然 Prony 等参数辨识方法可以应用在任意信号上，但如果要其物理意义符合得好，则应该假设信号是 ...

//     #text(fill: red, "（TODO）")也不太对，deterministic 信号哪来的平稳性？还是只是纯记号。
// ]

观察可知 $r_x (k)$ 是共轭对称的，即 $r_x (k) = r_x^* (-k)$。代入可以得到适用于 All-pole 模型的更简洁的方程组：

#emphasis_equbox([
$
sum_(l=1)^p a[l] r_x (k - l) = -r_x (k), quad k=1, 2, dots, p
$
])

或：

#emphasis_equbox([
$
mat(
    delim: "[",
    r_x (0), r_x^* (1), dots, r_x^* (p-1);
    r_x (1), r_x (0), dots, r_x^* (p-2);
    dots.v, dots.v, dots.down, dots.v;
    r_x (p-1), r_x (p-2), dots, r_x (0);
)
mat(
    delim: "[",
    a[1];
    a[2];
    dots.v;
    a[p];
)
=
-mat(
    delim: "[",
    r_x (1);
    r_x (2);
    dots.v;
    r_x (p);
)
$ <equ:dmi_prony_all_pole_matrix>
])

这被称为 *All-pole normal equations*。由于矩阵 $bold(R)_x$ 是共轭对称而且 Toeplitz 的，这使得我们可以使用 Levinson-Durbin 算法对它进行高效的求解。

=== Issues on the Numerator Selection

按常规方法，求得 $a[dot]$ 后我们会通过 @equ:dmi_b_X0a 得 $b[0] = x[0]$。但在 @equ:dmi_prony_normal_equations 的最后我们提到 Prony 法中分步求解的方法并不保证全局最优。我们在这里修改 $b[0]$ 的取值并不一定会破坏结果的最优性，因为结果本来就不是最优的；相反，我们甚至有可能通过换一种 $b[dot]$ 的取值方式达到更好的效果。

而所谓更好的效果也不一定是误差均方值的降低，也可能综合其他因素的考量。这里的全极点模型就是一个例子，如果原信号由于噪声或其他干扰导致 $x[0]$ 的值实际上不是那么可信，为了防止整个模型参数都被这一个 $b[0] = x[0]$ 带偏，我们更倾向于令拟合信号 $hat(x)[n]$（在我们的模型中就等于单位脉冲响应 $h[n]$）与目标信号 $x[n]$ 的能量相等：

$
r_(hat(x)) (0) = r_h (0) = r_x (0)
$

推导可以得到应取 $b[0] = sqrt(epsilon_(p))$。

#text(fill: red, "（TODO）")关于如何推导得出该 $b[0]$ 的取值还没搞懂，参考书中称在其 5.2.3 节会讲。

== Finite Data Records for All-pole Cases

前面对 Prony 法的分析都基于一个假设，即 $x[n]$ 定义在整个正时间域上，从 $0$ 到 $infinity$。而现在我们需要考虑当我们只拥有 $[0, N]$ 上样本的情况。

下面要介绍的两种思路通常用于全极点模型，所以我们也只默认讨论全极点模型。

=== Auto-correlation Method <sec:dmi_finite_data_autocorrelation_method>

第一种方法，我们考虑对 $x[n]$ 加矩形窗，或者说视 $x[n]$ 在 $[0, N]$ 以外的部分值为 $0$，然后仍然直接应用 Prony 法求解。

#text(fill: red, "（TODO）")

=== Covariance Method

实际上，Auto-correlation Method 把定义域以外的部分设为 $0$ 本质上是改变了 $x[n]$ 的形态，因为 $0$ 也是正常的信号值。而在一些情况下，这样做并不能达到最好的效果。

第二种 Covariance Method 的结果通常更加准确，其不对信号本身作假设，而是在优化过程中不考虑定义域以外的样本，更加自然。

#text(fill: red, "（TODO）")

// === Generalizations to the Covariance Method

// #text(fill: red, "（TODO）")
