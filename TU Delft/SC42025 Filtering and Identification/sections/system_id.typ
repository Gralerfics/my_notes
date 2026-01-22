#import "../generic.typ": *

#import "@preview/cetz:0.4.2"
#import "@preview/fletcher:0.5.8" as fletcher: diagram, node, edge

#import "@preview/ctheorems:1.1.3": *

= System Identification Cycle

#figure(
    caption: [System identification workflow]
)[
    #resize_box([
        #diagram(
            spacing: (12mm, 10mm),
            node-stroke: 0.8pt,
            edge-stroke: 0.8pt,

            // ----------------------------
            // Main flow (top row)
            // ----------------------------
            node((-12, 0), [Start], shape: "rect", inset: 6pt),

            node((-7, 0), [Experiment
    Design], shape: "rect", inset: 8pt),

            node((-2, 0), [Experiment], shape: "rect", inset: 8pt),

            node((3, 0), [Data pre-
    processing], shape: "rect", inset: 8pt),

            node((8, 0), [Fit model
    to data], shape: "rect", inset: 8pt),

            node((13, 0), [Model
    validation], shape: "rect", inset: 8pt),

            node((18, 0), [ok?], shape: "rect", inset: 8pt),

            node((23, 0), [End], shape: "rect", inset: 6pt),

            // top row arrows
            edge((-12, 0), (-7, 0), "-|>"),
            edge((-7, 0), (-2, 0), "-|>"),
            edge((-2, 0), (3, 0), "-|>"),
            edge((3, 0), (8, 0), "-|>"),
            edge((8, 0), (13, 0), "-|>"),
            edge((13, 0), (18, 0), "-|>"),
            edge((18, 0), (23, 0), "-|>", [Yes], label-pos: 0.6),

            // ----------------------------
            // Side box: Model structure selection (below middle)
            // ----------------------------
            node((5.5, -3.2), [Model
    structure
    selection], shape: "rect", inset: 8pt),

            // Upward arrow from selection to "Fit model to data"
            edge((5.5, -3.2), (8, 0), "-|>"),

            // ----------------------------
            // Feedback arrows to "Model structure selection"
            // (from Data pre-processing, Fit model, Validation)
            // ----------------------------
            // from Data pre-processing down to selection
            edge((3, 0), (3, -3.2), "-|>"),
            // from Fit model to data down to selection (as in figure, a short vertical)
            edge((8, 0), (8, -3.2), "-|>"),
            // from Model validation down to selection
            edge((13, 0), (13, -3.2), "-|>"),

            // connect those verticals into the selection box (horizontal feeds)
            edge((3, -3.2), (5.5, -3.2), "-|>"),
            edge((8, -3.2), (5.5, -3.2), "-|>"),
            edge((13, -3.2), (5.5, -3.2), "-|>"),

            // ----------------------------
            // Global "No" loop back (bottom frame-like loop)
            // ----------------------------
            // From decision diamond down
            edge((18, 0), (18, -6.2), "-|>", [No], label-pos: 0.4),

            // Bottom horizontal back to far left
            edge((18, -6.2), (-9.5, -6.2), "-|>"),

            // Go up to "Experiment Design" (left vertical)
            edge((-9.5, -6.2), (-9.5, 0), "-|>"),

            // "taps" up into the three boxes (as in the picture)
            // Up arrows into Experiment Design / Experiment / Data pre-processing
            edge((-7, -6.2), (-7, 0), "-|>"),
            edge((-2, -6.2), (-2, 0), "-|>"),
            edge((3, -6.2), (3, 0), "-|>"),

            // Make the bottom rail explicit between those tap points
            edge((-9.5, -6.2), (-7, -6.2), "-"),
            edge((-7, -6.2), (-2, -6.2), "-"),
            edge((-2, -6.2), (3, -6.2), "-"),
            edge((3, -6.2), (18, -6.2), "-"),
        )
    ])
]

#Cre("TODO Lec6")

= Prediction Error Methods

#Cre("TODO")

注意之后多处用 $G(q)$ 等指代 $G(q|theta)$ 等，只是为了简便起见，省略参数。

== Models

=== Transfer Function Model

#Cre("TODO") time-shift operator

$
x_(k+1) = q x_k \
x_(k-1) = q^(-1) x_k
$

#Cre("TODO") transfer function

$
y_k = G(q) u_k + H(q) e_k
$

#Cre("TODO") transfer function family

#Cre("TODO") simulation, use ${u_0, dots, u_N}$ to simulate ${y_0, dots, y_N}$

#Cre("TODO") prediction, use ${(u_0, y_0), dots, (u_(k-1), y_(k-1)), u_k}$ to predict $y_k$

// prediction 用 input 和 outputs; simulation 仅用 inputs (?)

$
y_k &= G(q) u_k + H(q) e_k \
H^(-1) (q) y_k &= H^(-1) (q) G(q) u_k + e_k
$

and

$
y_k &= hat(y)_k + e_k \
hat(y)_k &= y_k - e_k \
hat(y)_k &= y_k - e_k + H^(-1) (q) G(q) u_k - H^(-1) (q) G(q) u_k \
hat(y)_k &= H^(-1) (q) G(q) u_k + y_k - underbrace((H^(-1) (q) G(q) u_k + e_k), H^(-1) (q) y_k) \
hat(y)_k &= underbrace(H^(-1) (q) G(q), hat(G)(q)) u_k + underbrace([1 - H^(-1) (q)], hat(H)(q)) y_k
$

#Cre("TODO") 把参数标上，称其为 parameterized predictor:

#emphasis_equbox([
$
hat(y)_k (theta) &= underbrace(H^(-1) (q) G(q), hat(G)(q|theta)) u_k + underbrace([1 - H^(-1) (q)], hat(H)(q|theta)) y_k
$ <equ:id_tfmodel_param_predictor>
])

#blockquote([
    关于预测 $hat(y)_k$ 的这个公式，里面却使用了 $y_k$ 的问题。它的前面算到最后没问题的话是会有 $q$ 算子的，实际用到的是 $k-1$ 或其他更前的项目而非 $y_k$，所以没问题。
])

=== State Space Model (Innovation form)

#Cre("TODO") 见前 innovation form。

$
x_(k+1) &= A x_k + B u_k + K e_k \
y_k &= C x_k + D u_k + e_k
$ <equ:id_state_space_model_form>

#Cre("TODO") convert state-space model to transfer function

$
x_(k+1) &= A x_k + B u_k + K e_k \
q x_k &= A x_k + B u_k + K e_k \
(q I - A) x_k &= B u_k + K e_k \
x_k &= (q I - A)^(-1) B u_k + (q I - A)^(-1) K e_k
$

so 

$
y_k &= C x_k + D u_k \
y_k &= underbrace([C (q I - A)^(-1) B + D], G(q)) u_k + underbrace([C (q I - A)^(-1) K + I], H(q)) e_k
$

// 得到 G(q) 和 H(q)；二者是特征值xxx的不同线性变换。

#Cre("TODO") simulation

#Cre("TODO") prediction …… 最后是（？）parameterized predictor:

$
hat(x)_0 (theta), hat(x)_(k+1) = dash(A)(theta) hat(x)_k + dash(B)(theta) u_k + K(theta) y_k \
hat(y)_k (theta) = C(theta) hat(x)_k + D(theta) u_k
$

=== Summary

#Cre("TODO") 总结两种模型分别的 ground truth 和 predictor 形式，以及之间的关联性。

== Prediction-Error System Identification Components

=== Overview

#Cre("TODO") 1. measurements; 2. parameterized predictor $hat(y)_k (theta)$; 3. minimize $sum_(k=0)^(N-1) 1/N norm(y_k - hat(y)_k (theta))_2^2$. 多个组成：

=== Parameterization

#Cre("TODO") 状态空间方程改写，引入参数 theta；或者用传递函数模型，引入 theta

==== Example: SISO ARMAX

#Cre("TODO") 例子，SISO ARMAX，传递函数，对应状态空间模型（能观正规形）；其 parameterized predictor，以及参数集合；再如上改写

按 ARMAX 模型，设真实模型为：#Cre("TODO 这里 b_0 为 0 是 ARMAX 要求还是正规型要求")

$
y_k = (#Cre($0$) + b_1 q^(-1) + dots + b_n q^(-n))/(1 + a_1 q^(-1) + dots + a_n q^(-n)) u_k + (1 + c_1 q^(-1) + dots + c_n q^(-n))/(1 + a_1 q^(-1) + dots + a_n q^(-n)) e_k
$

将其转为能观正规型（状态空间模型）：

$
&A = mat(
    delim: "[",
    -a_1, 1, 0, dots, 0;
    -a_2, 0, 1, dots, 0;
    dots.v, dots.v, dots.v, dots.down, dots.v;
    -a_(n-1), 0, 0, dots, 1;
    -a_n, 0, 0, dots, 0;
),&
&B = mat(
    delim: "[",
    b_1;
    b_2;
    dots.v;
    b_(n-1);
    b_n;
),
K = mat(
    delim: "[",
    c_1 - a_1;
    c_2 - a_2;
    dots.v;
    c_(n-1) - a_(n-1);
    c_n - a_n;
), \
&C = mat(
    delim: "[",
    1, 0, 0, dots, 0;
),&
&D = 0.
$

如前所述，预测误差法最终是要最小化测量值 $y_k$ 与参数化预测子 $hat(y)_k (theta)$ 之间的误差。根据 @equ:id_tfmodel_param_predictor 等得到预测子：

$
hat(y)_k (theta) &= H^(-1) (q) G(q) u_k + [1 - H^(-1) (q)] y_k \
&= B(q) / C(q) u_k + (C(q) - A(q)) / C(q) y_k \
&= (b_1 q^(-1) + dots + b_n q^(-n))/(1 + c_1 q^(-1) + dots + c_n q^(-n)) u_k + ((c_1 - a_1) q^(-1) + dots + (c_n - a_n) q^(-n))/(1 + c_1 q^(-1) + dots + c_n q^(-n)) y_k
$

这个参数化的预测子共包含 $4n$ 个参数，列举如下：#Cre("TODO 最后的 hat y 是怎么个说法")

$
theta = {a_1, dots, a_n, b_1, dots, b_n, c_1, dots, c_n, hat(y)_(-1), dots, hat(y)_(-n)}
$

或者使用等价的状态空间模型描述该预测子：#Cre("TODO")

$
&dash(A) (theta) = mat(
    delim: "[",
    -c_1, 1, 0, dots, 0;
    -c_2, 0, 1, dots, 0;
    dots.v, dots.v, dots.v, dots.down, dots.v;
    -c_(n-1), 0, 0, dots, 1;
    -c_n, 0, 0, dots, 0;
),&
&dash(B) (theta) = mat(
    delim: "[",
    b_1;
    b_2;
    dots.v;
    b_(n-1);
    b_n;
),&
&K (theta) = mat(
    delim: "[",
    c_1 - a_1;
    c_2 - a_2;
    dots.v;
    c_(n-1) - a_(n-1);
    c_n - a_n;
), \
& C = mat(
    delim: "[",
    1, 0, 0, dots, 0;
),&
&D = 0, \
&hat(x)_n (theta)
$

==== Example: SISO Box-Jenkins

#Cre("TODO") 例子，SISO Box-Jenkins 类似来一遍

==== Notes

首先，对于使用能观正规型进行参数化的预测子，极点可能对 $dash(A)$ 和 $K$ 的变化很敏感。解决方案是：#Cre("TODO") Alternative parameterization: Output normal form (book page 220 and 259)

还有一种极端一点的方式是采用完全参数化（full parameterization），把矩阵每一项都用参数填满，模型灵活度极高，但参数量过大。

=== Optimization Criterion

// 模型，参数分布（贝叶斯），似然函数，求最大值；非凸，可视化；……

=== Numerical Optimization

// Nonlinear least squares problem

==== Gradient Descent

// ...

stable initial guess

constructing a stable conjugate pole pair $p_1 = abs(z) e^(i alpha), p_2 = abs(z) e^(-i alpha)$;
+ sample $0 < abs(z) < 1$ from uniform distribution
+ ……
+ $(z-p_1)(z-p_2)$

$
(q^(-2))/(1-2abs(z) cos(alpha) q^(-1) + abs(z)^2 q^(-2))
$

Slides P6, 7, and $J_N = 1/N E_N^T E_N$, so the gradient ……

efficient way to compute Phi_N and E_N? (theta)

1. run the predictor to get ${hat(y)_k (theta^((i)))}$, then obtain ${e_k (theta^((i)))}$, i.e. $E_N (theta^((i)))$;
2. 按定义求 Phi_N（measurements 是常数可以去掉）;
3. 考虑状态空间模型，求导过程 P10;
4. ……

总之只是避免重复计算那回事。（但如果参数化后求导里面还有参数那每次都得再算了，最好不要这么参数化）

TODO 例子 ……

$
dash(A) = mat(theta_1, 1; theta_2, 0), dash(B) = vec(theta_3, theta_4), C = mat(1, 0), D = 0, K = vec(theta_5, theta_6), hat(x)_0 = vec(theta_7, theta_8)
$

$
(partial hat(y)_1)/(partial theta_1) = theta_7, (partial hat(y)_2)/(partial theta_1) = 2 theta_1 theta_7 + theta_8 + theta_3 u_0 + theta_5 y_0
$

==== (Regularized) Gauss-Newton

用逆 Hessian 替换步长 $mu$（几何意义）。

Hessian 计算复杂，求导后舍去小项近似；以防奇异无法求逆，加一个小量。

=== Analysis of the Estimate

例子，1. 对 theta 线性; 2. 对 theta 非线性（几何解释，方差里的 Hessian）; 3. 模型和真值不符。（在说什么估计量？前面的梯度下降？）

计算偏差和方差；一致估计

MSE 的 bias-variance 分解

= Subspace Methods

子空间法是完全不同于预测误差法的另一套线性系统辨识方法。接下来我们先阐述具体方法，之后再说明子空间法背后的思路。

在下面的分析中，我们考虑前面 @equ:id_state_space_model_form 中的状态空间模型形式，其中包括了状态转移模型、输入信号以及噪声对模型输出信号的影响。我们从只考虑状态转移开始，逐步加入对特定输入信号、任意输入信号乃至噪声的处理能力。

== Discussion: Identification of an Autonomous System

为了梳理方法的整体脉络，先从简化的模型开始。考虑一个确定的自动系统（deterministic autonomous system），即不受输入序列 ${u_k}$ 和噪声的影响，只有内部状态在随时间转移的系统：

$
x_(k+1) &= A x_k \
y_k &= C x_k
$

此时很容易沿时间展开，将输出写成关于 $A$、$C$ 和 $x_0$ 的表达式：

$
y_k = C A^k x_0
$ <equ:id_autosys_output_x0>

回到问题上，我们是已知输出序列 ${y_k}$，希望去估计系统参数和状态序列 ${x_k}$。先不考虑系统参数，如果要估计状态序列 ${x_k}$，我们一定希望*先找到未知状态序列 ${x_k}$ 和已知输出数据之间的联系*。观察前面的 @equ:id_autosys_output_x0 会发现它已经提供了这一关系，我们将一段长度为 $s$ 的输出序列 ${y_0, y_1, dots, y_(s-1)}$ 用该式代换一下：

#let y_s_mat = $
    mat(
        delim: "[",
        y_0;
        y_1;
        y_2;
        dots.v;
        y_(s-1);
    )
$;

#let Os_mat = $
    mat(
        delim: "[",
        C;
        C A;
        C A^2;
        dots.v;
        C A^(s-1);
    )
$;

#let Ts_mat = $
    mat(
        delim: "[",
        D, 0, 0, dots, 0;
        C B, D, 0, dots, 0;
        C A B, C B, D, dots, 0;
        dots.v, dots.v, dots.down, dots.down, dots.v;
        C A^(s-2) B, C A^(s-3) B, dots, C B, D;
    )
$;

#let u_s_mat = $
    mat(
        delim: "[",
        u_0;
        u_1;
        u_2;
        dots.v;
        u_(s-1);
    )
$;

$
#y_s_mat = #Os_mat x_0 =: cal(O)_s x_0
$

// TODO 前面只取长度为 $s$ 的一段窗口感觉就是考虑矩阵过大

这里是一段输出序列同 $x_0$ 之间的关系，那么后面的 $x_1$ 和 $x_2$ 到 $x_(N-1)$ 呢？只要时移这段输出序列窗口即可：

$
mat(
    delim: "[",
    y_k;
    y_(k+1);
    y_(k+2);
    dots.v;
    y_(k+s-1);
)
=
mat(
    delim: "[",
    C A^k;
    C A^(k+1);
    C A^(k+2);
    dots.v;
    C A^(k+s-1);
) x_0
=
mat(
    delim: "[",
    C;
    C A;
    C A^2;
    dots.v;
    C A^(s-1);
) A^k x_0
=
cal(O)_s x_k
$ <equ:id_autosys_shifted_y_x_k>

于是，若我们把未知状态序列写成 $n times N$ 的矩阵（状态维数 $n$，共$N$ 个时刻的数据）：

#let Y_0sN_mat = $
    mat(
        delim: "[",
        y_0, y_1, dots, y_(N-1);
        y_1, y_2, dots, y_N;
        dots.v, dots.v, dots.down, dots.v;
        y_(s-1), y_s, dots, y_(N+s-2);
    )
$;

#let X_0N_mat = $
    mat(
        delim: "[",
        x_0, x_1, dots, x_(N-1);
    )
$;

$
#X_0N_mat
$

则*状态序列与输出数据之间的关系*一样通过水平拼接 @equ:id_autosys_shifted_y_x_k 得到：

$
underbrace(#Y_0sN_mat, Y_(0,s,N)) = cal(O)_s underbrace(#X_0N_mat, X_(0,N))
$

这被称作数据方程，即 *data equation*，是*在已知和未知之间建立起的关系*，有了该式就有依据对系统进行辨识。可以发现该式中输出数据构成了一个 Hankel 矩阵 $Y_(0,s,N)$，各副对角线上元素相等，相当密集。

#blockquote([
    *关于状态变换对数据方程形式无影响的说明*：

    #Cre("TODO")

    #Cre("TODO") 这说明这套子空间法对状态变换不敏感，状态定义具体是什么无所谓，……
])

以上可认为是*关于想到将 $y_k$ 拼接成数据矩阵来操作的这一整套方法的缘由*，不是 "为什么这么拼 $y_k$"，而是 "要考虑到所有的未知 $x_k$ 所以自然地得到了这种 Hankel 矩阵的形式"。

*接下来*就是怎么通过这个关系来进行系统辨识。

#Cre("TODO") 秩的问题，以及需要添加的假设。

#Cre("TODO") 奇异值分解，对应赋值，各种方法都行，因为状态变化无所谓吗（io行为等价的状态空间模型本质都是状态定义不同？）？

// #blockquote([
//     *还有一些东西*：

//     一堆窗口横向再拼，此时拼的就是不同的 $x_k$ 对应的信息。这个过程像*卷积*！时间动力学系统？？？？？？？？？？？？？？？？？？？？？？？？？？？？？？？？？？？？？？？？？？？？？？？？？？？？？？？？？？？？？？？？？？？？？？？？？？？？？？？？？？？

//     还有子空间法后面依赖系统的平移不变性，是用这个平移一格来算 $A$ 的。。。。。。。。。。。。。。。。。。。。。

//     $s$ 的影响？越大越好？？？？？？？？？？？？？？？
// ])

== Discussion: Deterministic SISO Identification with Impulse Response (Ho-Kalman)

#Cre("TODO Deterministic 单输入单输出系统在已知输入为冲激响应下的辨识")

上节中我们忽略了输入序列对系统的影响，

TODO 这个真的有必要吗，直接下节？

== Identification with Deterministic Inputs

#Cre("TODO 换成任意输入；列出式子发现多出和输入相关的项；用正交投影矩阵消除项；")

模型：

$
x_(k+1) &= A x_k + B u_k \
y_k &= C x_k + D u_k
$

利用输入序列 ${u_k}$ 可以迭代展开得到模型的模拟（simulation）输出：

$
y_k = C A^k x_0 + sum_(j=0)^(k-1) C A^(k-1-j) B u_j + D u_k
$

由该式将一段窗口长度为 $s$ 的 $y_k$ 序列写到一个向量中得到：

$
#y_s_mat = underbrace(#Os_mat, cal(O)_s) x_0 + underbrace(#Ts_mat, cal(T)_s) #u_s_mat
$

滑动这个窗口得到新的 $y$ 向量并水平拼接，得到更大的表达式 #Cre("TODO 拼接等价的直觉强化")：

#let U_0sN_mat = $
    mat(
        delim: "[",
        u_0, u_1, dots, u_(N-1);
        u_1, u_2, dots, u_N;
        dots.v, dots.v, dots.down, dots.v;
        u_(s-1), u_s, dots, u_(N+s-2);
    )
$;

$
underbrace(#Y_0sN_mat, Y_(0,s,N)) = cal(O)_s underbrace(#X_0N_mat, X_(0,N)) + cal(T)_s underbrace(#U_0sN_mat, U_(0,s,N))
$

即：

$
Y_(0,s,N) = cal(O)_s X_(0,N) + cal(T)_s U_(0,s,N)
$

后文推导中默认简化符号如下：

$
Y = cal(O)_s X + cal(T)_s U
$ <equ:id_ss_data_equation_with_deterministic_inputs_simplified>

和前面讨论中不同的是，现在我们多出了 $cal(T)_s U$ 这一同输入有关的项，没法再像之前一样直接对 $Y$ 进行奇异值分解然后对应求解 $cal(O)_s X$。

于是自然地，我们希望抹掉这个输入相关项。考虑一个可以用已知的 $U$ 算出来的东西：

$
P_U^perp = I_N - U^T (U U^T)^(-1) U
$

这是一个正交投影矩阵，将向量投影到 $U$ 的行空间的正交补上，即有：

$
U P_U^perp = 0
$

我们给 @equ:id_ss_data_equation_with_deterministic_inputs_simplified 两边乘上它，得到：

$
Y P_U^perp = cal(O)_s X P_U^perp + cal(T)_s (U P_U^perp) = cal(O)_s X P_U^perp
$

于是输入相关项就没了。#Cre("TODO 这里还需要额外的假设，比如确保这个 P 能算")

#blockquote([
    *关于这里正交投影矩阵、行空间和正交补的具体说明*：

    首先对于一个矩阵 $P$，如果它幂等（$P^2 = P$）且共轭对称（$P^H = P$，考虑实数 $P$ 时也就是 $P^T = P$），则可称其为正交投影矩阵（orthogonal projection matrix）。



    TODO 关于，输出中的输入项，即任何输入 u 造成的输出变化，都在 u 张成的线性空间里，或者说本质上是 $y_u = G * u$，一个线性卷积；所以输入可以解释的输出，都属于一个由 u 决定的子空间（分析一下行空间）。

    然后，我们把输出中 u 能解释的部分即 y_u，直接投影消去（用 u 去最大程度拟合 y，然后减去），剩下的就是输入解释不了的部分，或者说输入无关的部分，那么就可以去除输入序列的影响，留下模型真正的状态结构。
    
    （输入中的持续激励？因输入而导致 Y 的列空间中出现的高维满秩结构会淹没正在的状态结构，which 一般是低维结构）。

    子空间法恢复的不是具体的状态（因为状态可以随便用 T 变换），而是状态空间（状态在输出中留下的几何痕迹）（？？？为什么这个不变？）。

    前面去除了输入影响，现在再找点云落在的低维平面（有点像 PCA）。实际过程中可能就是忽略掉过小的奇异值（或者已知模型维度的话就更方便了），来察觉到原系统到底是几维的状态，然后投影上去。所以这玩意应该对噪声是由抵抗力的。

    最接近第一性原理的线性系统辨识方法？
])

#Cre("TODO 需要一个重要结论") $"col"(Y P_U^perp) = "col"(cal(O)_s)$

…… 对 $vec(U, Y)$ 作 LQ 分解：

$
mat(delim: "[", U; Y) = mat(delim: "[", L_11, 0; L_21, L_22) mat(delim: "[", Q_1^T; Q_2^T)
$

#Cre("TODO 其中") $Q_1$ 是 $"row"(U)$ 的一组正交基，$Q_2$ 是 $"row"(U)^perp$ 的一组正交基。（？）

…… 反正一通 ……

最后会知道 $"col"(L_22) = "col"(cal(O)_s)$，即 $L_22$ 的列空间和 $cal(O)_s$ 相同，这意味着对 $L_22$ 进行奇异值分解就可以得到 $cal(O)_s$ 的一组基。

它们或许并不逐元素相同，但因为列空间相同，之间就差一个可逆变换，#Cre("TODO 我们可以辨识出和以前模型不一样 ABCD 但只是差一个状态变换的系统参数，只是状态定义变了，但功能正确，这样就行了")。

奇异值分解后非零奇异值个数就是 $n$，再通过和前面讨论中类似的方法，移位相除得到 $A$，取首行得到 $C$。由此我们通过 $L_22$ 辨识得 $n$、$A$ 和 $C$。

#Cre("TODO 然后？？？") 辨识得 $B$、$D$ 和 $x_0$。

以上即为 MOESP 法？

== TODO

#Cre("TODO 考虑噪声")

== Summary

#Cre("TODO")
