#import "../generic.typ": *

#import "@preview/cetz:0.4.2"
#import "@preview/fletcher:0.5.8" as fletcher: diagram, node, edge

#import "@preview/ctheorems:1.1.3": *

= System Identification Cycle

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

=== State Space Model

#Cre("TODO 有件重要的事，就是说明这个状态空间模型为什么要写成这个形式，有什么好处（指状态里面用 K e_k 描述噪声，输出里用一个 e_k）")

$
x_(k+1) &= A x_k + B u_k + K e_k \
y_k &= C x_k + D u_k + e_k
$

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

#Cre("TODO")

== Discussion: Identification of an Autonomous System

#Cre("TODO 自动系统（即无输入，只有 A 和 C）的辨识")

== Discussion: Deterministic SISO Identification with Impulse Response (Ho-Kalman)

#Cre("TODO Deterministic 单输入单输出系统在已知输入为冲激响应下的辨识")

== Identification with Deterministic Inputs

#Cre("TODO 换成任意输入；列出式子发现多出和输入相关的项；用正交投影矩阵消除项；")

模型：

$
x_(k+1) &= A x_k + B u_k \
y_k &= C x_k + D u_k
$

若已知输入序列 ${u_k}$（#Cre("TODO 初值貌似也是要辨识的")），迭代展开可得模拟输出：

$
y_k = C A^k x_0 + sum_(j=0)^(k-1) C A^(k-1-j) B u_j + D u_k
$

由该式将一段窗口长度为 $s$ 的 $y_k$ 序列写到一个向量中得到：

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
#y_s_mat = underbrace(#Os_mat, cal(O)_s) x_0 + underbrace(#Ts_mat, cal(T)_s) #u_s_mat
$

滑动这个窗口得到新的 $y$ 向量并水平拼接，得到更大的表达式 #Cre("TODO 拼接等价的直觉强化")：

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

这一矩阵是投影矩阵 #Cre("TODO 定义，证明")，投影到 #Cre("TODO 表达") 矩阵 $U$ 的行空间的正交补上 #Cre("TODO 问号")，即有：

$
U P_U^perp = 0
$

我们给 @equ:id_ss_data_equation_with_deterministic_inputs_simplified 两边乘上它，得到：

$
Y P_U^perp = cal(O)_s X P_U^perp + cal(T)_s (U P_U^perp) = cal(O)_s X P_U^perp
$

于是输入相关项就没了。#Cre("TODO 这里还需要额外的假设，比如确保这个 P 能算")

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
