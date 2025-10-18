#import "../generic.typ": *

#import "@preview/fletcher:0.5.8" as fletcher: diagram, node, edge

= Signal Modelling

== Motivation

// In fact, modelling is a process of compressing. A parameterized model uses finite parameters to represent complex systems. For signals, we can define a signal model as a mathematical function that can generate signals similar to the target signal.

实际上，建模（Modelling）可以视为对信号的压缩（Compressing）过程。一个参数化的模型可以使用比信号样本数更少的参数数量来表示复杂的信号，实现更高效的存储和传输。

压缩得到的参数也可以视为对信号本质特征的描述，例如蕴含其背后的物理规律等，这就允许我们利用模型对信号的未知部分进行预测（Prediction），或称外推（Extrapolation）。

现在，我们的目标是对给定的数字信号 $x[n]$ 进行建模，即找到一个模型 $H(z)$ 使其输出信号 $hat(x)[n]$ 能够尽可能接近目标信号 $x[n]$。

== Auto-regressive and Moving Average (ARMA) Model <sec:signal_modelling_arma_model>

我们可以根据实际情况使用不同种类的模型来对信号建模，这里以时间序列分析常用的自回归滑动平均模型（Auto-regression and Moving Average Model）即 $"ARMA"(p, q)$ 为例。其传递函数定义如下：

$ Y(z) / X(z) = H(z) = (sum_(k=0)^q b[k] z^(-k)) / (1 + sum_(k=1)^p a[k] z^(-k)) = B(z) / A(z) $ <equ:signal_modelling_arma_tf>

设输入和输出信号分别为 $x[n]$ 和 $y[n]$，$X(z)$ 和 $Y(z)$ 为其对应的拉普拉斯变换函数。这里形式上虽然 $a[n]$ 的索引从 $1$ 开始，但实际上可以取 $a[0] = 1$ 以得到更统一的形式。于是由定义我们有 $Y(z) A(z) = X(z) B(z)$，变换到时域即：

$ a[n] * y[n] = b[n] * x[n] $

展开得到经典的线性常系数差分方程（Linear Constant Coefficient Difference Equation，LCCDE）的形式：

$ y[n] + sum_(k=1)^p a[k] y[n-k] = sum_(k=0)^q b[k] x[n-k] $

写详细一点为：

$ y[n] + a[1] y[n-1] + ... + a[p] y[n-p] = b[0] x[n] + b[1] x[n-1] + ... + b[q] x[n-q] $

显然，这个系统是典型的线性移不变（Linear Shift-Invariant，LSI）系统，具有良好的性质。

=== Autoregressive (AR) Model (All-Pole Model)

若只有自回归的部分，即 $"AR"(p)="ARMA"(p, 0)$：

$ H(z) = Y(z) / X(z) = b[0] / (1 + sum_(k=1)^p a[k] z^(-k)) $ <equ:ar_model_transfer_function>

在时域即：

$ y[n] + sum_(k=1)^p a[k] y[n-k] = b[0] x[n] $

写详细一点为：

$ a[0] y[n] + a[1] y[n-1] + ... + a[p] y[n-p] = b[0] x[n] $

该模型认为当前时刻的输出 $y[n]$ 是前 $p$ 个时刻的输出 $y[n-1], ..., y[n-p]$ 以及当前时刻的输入 $x[n]$ 的线性组合，所以称为自回归模型。由于该模型无零点，故也被称为全极点模型（All-Pole Model）。

=== Moving Average (MA) Model (All-Zero Model)

若只有滑动平均的部分，即 $"MA"(q)="ARMA"(0, q)$：

$ H(z) = Y(z) / X(z) = sum_(k=0)^q b[k] z^(-k) $

在时域即：

$ y[n] = sum_(k=0)^q b[k] x[n-k] $

写详细一点为：

$ y[n] = b[0] x[n] + b[1] x[n-1] + ... + b[q] x[n-q] $

该模型认为当前时刻的输出 $y[n]$ 是前 $q$ 个时刻的输入 $x[n], x[n-1], ..., x[n-q]$ 的线性组合，所以称为滑动平均模型。

== Signal Models

作为一个离散时间系统，$H(z)$ 并不能直接表示信号，而是需要接受一个输入以获得输出。我们将输出信号作为模型对目标信号的估计 $hat(x)[n]$，将输入信号 $x[n]$ 封装为模型的一部分。由于该模型无极点，故也被称为全零点模型（All-Zero Model）。

=== Deterministic Modelling

#figure(
    caption: [Signal model with deterministic input]
)[
    // #cetz.canvas({
    //     import cetz.draw: *

    //     set-style(
    //         stroke: 0.8pt,
    //         mark: (transform-shape: false)
    //     )

    //     content(
    //         (0, 0),
    //         [$H(z)$],
    //         frame: "rect",
    //         padding: (x: 8pt, y: 6pt),
    //         name: "sys",
    //     )

    //     line((-3, 0), "sys.west", mark: (end: "stealth"))
    //     content((-3, 0), [Impulse $delta[n]$], anchor: "east", padding: 5pt)

    //     line("sys.east", (3, 0), mark: (end: "stealth"))
    //     content((3, 0), [$hat(x)[n] = h[n]$], anchor: "west", padding: 5pt)
    // })

    #diagram(
        spacing: (10mm, 8mm),
        node-stroke: 0.8pt,
        edge-stroke: 0.8pt,
        node((0, 0), [$H(z)$], inset: 8pt),
        edge((-3, 0), "r,r,r", "-|>", [Impulse $delta[n]$], label-pos: 0),
        // edge((-3, 0), (0, 0), "-|>", [Impulse $delta[n]$], label-pos: 0),
        edge((0, 0), (3, 0), "-|>", [$hat(x)[n] = h[n]$], label-pos: 1),
    )
] <fig:signal_model_deterministic>

我们可以选择一个已知的、确定的信号，将其固定为系统的输入，稳定得到我们想要的输出信号 $hat(x)[n]$，使其值尽可能接近目标信号 $x[n]$，用于对确定的信号进行建模，称为 Deterministic Modelling。

这个输入信号可以根据实际情况进行选择，符合目标信号的特征的输入信号有时可以减轻模型拟合的负担。在这里我们可以选择使用最简单的单位脉冲信号 $delta[n]$ 作为输入信号，这使得系统的输出信号 $hat(x)[n]$ 即为系统的单位冲激响应 $h[n]$。

=== Stochastic Modelling

#figure(
    caption: [Signal model with stochastic input]
)[
    #diagram(
        spacing: (10mm, 8mm),
        node-stroke: 0.8pt,
        edge-stroke: 0.8pt,
        node((0, 0), [$H(z)$], inset: 8pt),
        edge((-3, 0), (0, 0), "-|>", [White noise $v[n]$], label-pos: 0),
        edge((0, 0), (3, 0), "-|>", [$hat(x)[n]$], label-pos: 1),
    )
] <fig:signal_model_stochastic>

我们还可以选择使用一个已知分布的随机噪声作为输入，得到输出信号 $hat(x)[n]$ 使其统计特征（例如均值、自相关函数）与目标信号 $x[n]$ 一致，从而对随机过程进行建模，称为 Stochastic Modelling。

我们可以选择使用均值为 $0$、方差为 $sigma_v^2$ 的白噪声 $v[n]$ 作为输入信号。这样做的依据是其自相关函数为 $r_v [k] = sigma_v^2 delta[k]$，对其进行傅里叶变换得其功率谱密度为常数 $P_v (omega) = sigma_v^2$，即在所有频率上均有相同的能量分布。

这样的特性保证了我们可以通过对其进行滤波得到任意频率特性的输出信号 $hat(x)[n]$，同时任意频率成分能量均匀，使得输出信号的统计特征与输入信号无关。
