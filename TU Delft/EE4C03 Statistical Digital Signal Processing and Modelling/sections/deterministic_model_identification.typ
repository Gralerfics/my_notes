#import "../generic.typ": *

#import "@preview/cetz:0.4.2"
#import "@preview/fletcher:0.5.8" as fletcher: diagram, node, edge

= Deterministic Modelling Identification

显然，模型建立后我们还需要对模型参数进行辨识，即确定参数序列 $a[n]$ 和 $b[n]$ 的值。而选取参数的目标是使模型输出 $hat(x)[n]$ 尽可能接近目标信号 $x[n]$。

== Least Squares (LS) Method

首先讨论 Deterministic Modelling，我们希望模型输出 $hat(x)[n]$ 能够精确地重现目标信号 $x[n]$，即每个采样点的值都要尽可能接近。定义误差信号 $e'[n] = x[n] - hat(x)[n]$，可以通过最小化均方误差 $cal(E)_"LS" = sum_(n=0)^infinity norm(e'[n])^2$ 来确定模型参数，如 @fig:deterministic_model_identification_diagram_ls_intractable 所示。

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
] <fig:deterministic_model_identification_diagram_ls_intractable>

该优化问题可通过求解如下方程组解决：

$
cases(
    frac(partial cal(E)_"LS", partial a^*[k]) = 0\, k = 1\, 2\, ...\, p,
    frac(partial cal(E)_"LS", partial b^*[k]) = 0\, k = 0\, 1\, ...\, q
)
$

但 $H(z)$ 包含分式，形式复杂，直接求解上述方程组是非线性问题。所以我们进行一些修改，在两路上同乘 $H(z)$ 的分母 $A(z)$。即令新的目标误差 $E(z) = A(z) E'(z) = A(z) X(z) - B(z)$，得到新方案如 @fig:deterministic_model_identification_diagram_ls 所示。

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
] <fig:deterministic_model_identification_diagram_ls>

具体地，单位脉冲信号经过 $B(z)$ 得到的输出就是 $b[n]$，将其与用目标信号 $x[n]$ 经过 $A(z)$ 得到的 $hat(b)[n]$（此时可视为是对系数序列 $b[n]$ 的估计）作差，得到新的误差。经过这样操作后的方程求解就是线性的了，Padé Approximation、Prony's Method 和 Shank's Method 都是基于该思路的技术。

== From a More Formalised Perspective

这一节从形式化的角度来理解，对应上一节中系统框图的变化。

首先再次明确目标，我们已知 $x[n]$，希望求解最优的 $a[n]$ 和 $b[n]$，使得模型输出 $hat(x)[n]$ 尽可能接近 $x[n]$；同时，为了方便求解，我们希望问题是线性的。也就是说，我们希望得到一个直接包含 $a[n]$、$b[n]$ 和 $x[n]$ 的线性方程组。

首先依然考虑 @fig:deterministic_model_identification_diagram_ls_intractable 中的系统，显然最理想情况下我们希望误差直接为零，即 $e'[n] = 0$，也就是使 $h[n] = x[n]$。但 $h[n]$ 和参数 $a[n]$、$b[n]$ 的关系是非线性的，没法通过 $h[n] = x[n]$ 列出一个好处理的方程，所以我们在此之前先作一个转化：

$ H(z) = B(z) / A(z) quad => quad H(z) A(z) = B(z) quad => quad h[n] * a[n] = b[n] $

这个两侧同乘 $A(z)$ 的操作本质上就是上一节中将 @fig:deterministic_model_identification_diagram_ls_intractable 中的系统转化为 @fig:deterministic_model_identification_diagram_ls 中的系统的过程。将卷积展开得到：

$ h[n] + sum_(k=1)^p a[k] h[n-k] = b[n] $

前面说了，我们希望使误差直接为零，所以带入 $h[n] = x[n]$，得到：

$ x[n] + sum_(k=1)^p a[k] x[n-k] = b[n] $ <equ:deterministic_model_identification_equation>

这就是我们想要的线性方程组了，包含了 $a[n]$、$b[n]$ 和 $x[n]$。但该方程组中未知数的个数是 $p + q + 1$，而方程的个数是无穷多个（实际上取决于信号的长度，系数均为零的等式就无意义了），所以该方程组是一个超定（Overdetermined）方程组。

#blockquote[
    上一节中应用最小二乘法将问题视为优化问题处理，而这里我们更直接地令误差为零，得到一个超定方程组。

    对后者我们可以通过伪逆来求最优解，这本质上还是最小二乘法，但将显性的求最小化优化问题的过程，转化为隐性的、公式化的将理想解投影到可行解空间的过程，让前面的思路更加清晰一些。
]

当然，我们也可以直接从 @fig:deterministic_model_identification_diagram_ls 中的系统出发列出方程组。写出 $hat(b)[n]$ 的表达式：

$ hat(b)[n] = x[n] + sum_(k=1)^p a[k] x[n-k] $

直接令 $e[n] = 0$ 即 $hat(b)[n] = b[n]$，得到和 @equ:deterministic_model_identification_equation 一样的方程组，这里属于废话了。

== Padé Approximation

接下来就是对 @equ:deterministic_model_identification_equation 的求解了。对该方程组的求解方法中 Padé Approximation 是最简单的一种，只取 $n = 0, 1, ..., p + q$，得到前 $p + q + 1$ 个方程，和未知数个数相同，非奇异的话刚好可以求解出唯一解。

清晰一点，我们把矩阵画出来：

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

        sub_box(-6.25, -1.3, 9.485, 4.65, blue, $bold(X)_0$, 0, 2.6, blue)
        sub_box(-6.25, -1.4, 1.55, -1.95, purple, $bold(x)_(q+1)$, 0, -1.35, purple)
        sub_box(-4.6, -1.4, 7.835, -1.95, purple, $bold(X)_q$, 0, -1.35, purple)
        
        sub_box(3.7, -1.64, 0.79, 2.62, green, $dash(bold(a))$, 0, -1.6, green)
    })
]

在这里，我们先用下半部分即后 $p$ 行求解 $dash(bold(a))$（即 $a[n]$）：

$ bold(X)_q dash(bold(a)) = - bold(x)_(q+1) quad => quad dash(bold(a)) = - bold(X)^(-1)_q bold(x)_(q+1) $ <equ:deterministic_model_identification_pade_a_bar>

然后代入上半部分即前 $q+1$ 行计算 $b[n]$：

$ bold(b) = bold(X)_0 bold(a) = bold(X)_0 vec(1, dash(bold(a))) $

很直接，但显然也存在一系列问题：
+ 不保证所得系统是稳定的，极点（即 $1 + sum_(k=1)^p a[k] z^(-k) = 0$ 的根）可能在单位圆外；
+ 只约束了模型输出 $hat(x)[n]$ 和目标信号 $x[n]$ 的前 $p + q + 1$ 个样本相同，此后的匹配效果可能不佳；
+ 若 $bold(X)_q$ 奇异则可能无解。

== Prony's Method

=== Extensions of the Padé Approximation

Prony 的想法很简单，就是在 Padé 基础上，将 $x[n]$ 的所有样本点都纳入考虑。还是把矩阵画出来：

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

和 @equ:deterministic_model_identification_matrix_pade 类似，只是矩阵的行数向下无限延伸（直到 $x[n]$ 不再有值的地方）。$bold(x)_(q+1)$ 和 $bold(X)_q$ 也随之延伸，理论上若 $x[n]$ 是无限长的序列，则此二者也为无限维向量和矩阵，不过实际情况下我们需要处理的 $x[n]$ 基本都是有限的序列。由此我们通过引入整个 $x[n]$ 的信息，解决了 Padé 法的第二个问题。

#blockquote[
    由于模型的描述能力是有限的，在多数情况下使用 Prony 对整个序列拟合效果的提升往往也伴随着前 $p + q + 1$ 个样本相对 Padé 法拟合效果的下降。当然，实际上如果我们更关注序列某些特定部分的拟合效果的话，可以只将这些样本对应的等式纳入矩阵求解。
]

考虑到我们引入了更多的等式，使 $bold(X)_q dash(bold(a)) = - bold(x)_(q+1)$ 变成了一个超定方程组，需要借助伪逆来得到最优的解：

$ dash(bold(a)) = - text(fill: #purple, bold(X)^+_q) bold(x)_(q+1) = - text(fill: #purple, (bold(X)_q^H bold(X)_q)^(-1) bold(X)_q^H) bold(x)_(q+1) $

我们再记 $bold(R)_x = bold(X)_q^H bold(X)_q$ 和 $bold(r)_x = bold(X)_q^H bold(x)_(q+1)$，则有：

$ dash(bold(a)) = - (text(fill: #blue, bold(X)_q^H bold(X)_q))^(-1) text(fill: #purple, bold(X)_q^H bold(x)_(q+1)) = - text(fill: #blue, bold(R)_x^(-1)) text(fill: #purple, bold(r)_x) $

注意到对 $A^H A$ 这种形式的矩阵，对任意向量 $bold(a)$ 有 $bold(a)^H (A^H A) bold(a) = (bold(A a))^H (bold(A a)) = norm(bold(A a))^2 >= 0$，即 $A^H A$ 是半正定矩阵。由此，前述 $bold(R)_x = bold(X)_q^H bold(X)_q$ 也是（半）正定矩阵，#text(fill: red, "（TODO）")该性质将决定 $A(z)$ 是（临界）稳定的，由此解决 Padé 法的第一个问题。

此外，若 $bold(R)_x$ 为正定矩阵，则其特征值都是正数，即行列式不为零，得矩阵可逆，解存在，由此解决 Padé 法的第三个问题；#text(fill: red, "（TODO）")若 $bold(R)_x$ 为包含零特征值的半正定矩阵，则奇异，不可逆，但这实际上说明模型的阶数冗余了，可以降低一些再尝试。

=== Prony Normal Equations

显然，上节中我们记 $bold(R)_x = bold(X)_q^H bold(X)_q$ 和 $bold(r)_x = bold(X)_q^H bold(x)_(q+1)$ 是有原因的。首先，我们画出 $bold(R)_x$ 的计算过程：

$
bold(R)_x = // bold(X)_q^H bold(X)_q =
mat(
    delim: "[",
    x^*[q], dots, x^*[q+p-1], dots;
    x^*[q-1], dots, x^*[q+p-2], dots;
    dots.v, dots.down, dots.v, dots.down;
    x^*[q-p+1], dots, x^*[q], dots;
)
mat(
    delim: "[",
    x[q], x[q-1], dots, x[q-p+1];
    dots.v, dots.v, dots.down, dots.v;
    x[q+p-1], x[q+p-2], dots, x[q];
    dots.v, dots.v, dots.down, dots.v;
)
$ <equ:deterministic_model_identification_matrix_prony>

#text(fill: red, "（TODO）")（Rx、rx和相关矩阵）
