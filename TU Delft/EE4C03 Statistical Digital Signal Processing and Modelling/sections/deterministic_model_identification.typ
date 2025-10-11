#import "@preview/ilm:1.4.1": *

#import "@preview/fletcher:0.5.8" as fletcher: diagram, node, edge

#let section_deterministic_model_identification() = {[

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

具体地，单位脉冲信号经过 $B(z)$ 得到的输出就是 $b[n]$，将其与用目标信号 $x[n]$ 经过 $A(z)$ 得到的 $hat(b)[n]$（此时可视为是对系数序列 $b[n]$ 的估计）作差，得到新的误差。经过这样操作后的方程求解就是线性的了，Pade Approximation、Prony's Method 和 Shank's Method 都是基于该思路的技术。

== From a More Formalised Perspective

// 上一节将问题直接视为优化问题处理，这里我们更直接地令误差为零，大概率直接得到一个超定（Overdetermined）的方程组，然后再通过伪逆来求最优解。这本质上还是最小二乘法，但将显性的求最小化优化问题的过程，转化为隐性的将理想解投影到可行参数空间的过程，让前面的思路更加直接一些。或者也可以说，

这一节从形式化的角度来理解，对应上一节中系统框图的变化。

首先再次明确目标，我们已知 $x[n]$，希望求解最优的 $a[n]$ 和 $b[n]$，使得模型输出 $hat(x)[n]$ 尽可能接近 $x[n]$；同时，为了方便求解，我们希望问题是线性的。也就是说，我们希望得到一个直接包含 $a[n]$、$b[n]$ 和 $x[n]$ 的线性方程组。

首先依然考虑 @fig:deterministic_model_identification_diagram_ls_intractable 中的系统，显然最理想情况下我们希望误差直接为零，即 $e'[n] = 0$，也就是使 $h[n] = x[n]$。但 $h[n]$ 和参数 $a[n]$、$b[n]$ 的关系是非线性的，没法通过 $h[n] = x[n]$ 列出一个好处理的方程，所以我们在此之前先作一个转化：

$ H(z) = B(z) / A(z) quad => quad H(z) A(z) = B(z) quad => quad h[n] * a[n] = b[n] $

这个两侧同乘 $A(z)$ 的操作本质上就是上一节中将 @fig:deterministic_model_identification_diagram_ls_intractable 中的系统转化为 @fig:deterministic_model_identification_diagram_ls 中的系统的过程。将卷积展开得到：

$ h[n] + sum_(k=1)^p a[k] h[n-k] = b[n] $

前面说了，我们希望使误差直接为零，所以带入 $h[n] = x[n]$，得到：

$ x[n] + sum_(k=1)^p a[k] x[n-k] = b[n] $

这就是我们想要的线性方程组了，包含了 $a[n]$、$b[n]$ 和 $x[n]$。但该方程组中未知数的个数是 $p + q + 1$，而方程的个数是无穷多个（实际上取决于信号的长度，系数均为零的等式就无意义了），所以该方程组是一个超定（Overdetermined）方程组。

#blockquote[
    上一节中应用最小二乘法将问题视为优化问题处理，而这里我们更直接地令误差为零，得到一个超定方程组。

    对后者我们可以通过伪逆来求最优解，这本质上还是最小二乘法，但将显性的求最小化优化问题的过程，转化为隐性的、公式化的将理想解投影到可行解空间的过程，让前面的思路更加直接一些。
]

// 我们写出 $hat(b)[n]$ 的表达式：

// $ hat(b)[n] =  $

// 我们直接令 $e[n] = 0$，列出方程组：

== Pade Approximation

]}
