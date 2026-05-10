#import "../generic.typ": *

#import "@preview/cetz:0.4.2"
#import "@preview/fletcher:0.5.8" as fletcher: diagram, node, edge

= Networked Control Systems

== Introduction

== Modelling and Stability of Sampled-Data Systems

$
dots
stretch(harpoon.rt)^(u_k) "ZOH"
stretch(harpoon.rt)^(v(t)) "Plant"
stretch(harpoon.rt)^(xi(t)) "Sensor"
stretch(harpoon.rt)^(y_k "sampled at" s_k) "Controller"
dots
$

假设：
+ 常数采样间隔 $h$，采样点为 $s_k = k dot h$。
+ 事件驱动的控制器和执行器。
+ 全状态可测：$y_k = x_k := xi (k h)$，$xi (t)$ 是连续状态函数。
+ 控制动作 $u_k$ 为 $y_k$ 的函数。
+ 零阶保持（zero-order hold，ZOH），$v(t)$ 为 $u_k$ 经过 ZOH 后的连续输入。

考虑 Plant 处的如下 LTI 系统：

$
dot(xi) (t) &= A xi(t) + B v(t), &quad& forall t in RR^+ \
v(t) &= u_k, && "for" t in [s_k, s_(k+1))
$

状态为 $xi(t)$，输入为 $v(t)$，其解为：

$
xi(t) = e^(A(t-t_0)) xi(t_0) + integral_(t_0)^t e^(A(t-s)) B v(s) dif s
$

就是标准的状态空间模型状态解。因为对于 $t in [s_k, s_(k+1))$ 有 $v(t) = u_k$，取 $t_0 = s_k, t = s_(k+1)$ 代入即为采样所得离散状态：

$
x_(k+1) :&= xi(s_(k+1)) \
&= e^(A h) xi(s_k) + integral_(s_k)^(s_(k+1)) e^(A(s_(k+1)-s)) B u_k dif s \
&= e^(A h) x_k + (integral_0^h e^(A s) B dif s) u_k \
&=: F(h) x_k + G(h) u_k
$

即记 $F(h) = e^(A h)$ 和 $G(h) = integral_0^h e^(A s) B dif s$，也就是连续状态空间模型转为离散模型的公式，注意依赖采样间隔 $h$。

#blockquote([
    *关于矩阵指数函数的积分*：

    若 $A$ 可逆，可以计算一下 $integral_0^h e^(A s) dif s$。应用泰勒展开有：

    $
    integral_0^h e^(A s) dif s
    &= integral_0^h sum_(k=0)^infinity (A s)^k / (k!) dif s \
    &= sum_(k=0)^infinity lr(s (A s)^k / ((k+1)!) |)_0^h \
    &= sum_(k=0)^infinity h (A h)^k / ((k+1)!)
    $

    于是：

    $
    (integral_0^h e^(A s) dif s) A
    = sum_(k=0)^infinity (A h)^(k+1) / ((k+1)!)
    = sum_(k=1)^infinity (A h)^k / (k!)
    = e^(A h) - I
    $

    得到：

    $
    integral_0^h e^(A s) dif s = (e^(A h) - I) A^(-1)
    $
])

#blockquote([
    *关于离散时间线性系统的稳定性*：

    对于线性离散时间系统 $x_(k+1) = macron(A) x_k$，若存在 $c > 0$ 和 $rho in (0, 1)$ 使得：
    
    $
    norm(x_k) <= c rho^(k-k_0) norm(x_0)
    $

    则其原点是一个全局指数稳定（globally exponentially stable，GES）不动点。其还有一个充分必要条件是：$macron(A)$ 的谱半径（spectral radius）小于 $1$，即：

    $
    rho(macron(A)) := max{abs(lambda_1 (macron(A))), abs(lambda_2 (macron(A))), dots, abs(lambda_n (macron(A)))} < 1
    $
])

现在回到先前的离散化系统 $x_(k+1) = F(h) x_k + G(h) u_k$，设计一个静态状态反馈控制律 $u_k = -macron(K) x_k$ 得到闭环系统：

$
x_(k+1) = (F(h) - G(h) macron(K)) x_k
= (e^(A h) - integral_0^h e^(A s) B macron(K) dif s) x_k
$

判断原点是否 GES 的条件即：

$
rho (F(h) - G(h) macron(K)) < 1
$

显然它受采样间隔 $h$ 的影响。通过一些例子还可以发现，随着采样间隔增大系统通常会变得不稳定，但有时继续增大由于频谱混叠（aliasing）反而可能稳定下来。不过这种稳定通常不代表表现好，实际波形在（可见的）采样点之间可能会发生（不可见的）高频震荡。

#blockquote([
    *关于离散时间模型稳定性和采样系统稳定性的关联*：

    TODO
])

== NCS With Delays

$
dots
stretch(harpoon.rt)^("controller-to-actuator delay" tau^"ca ")
stretch(harpoon.rt)^(u_k) "ZOH"
stretch(harpoon.rt)^(v(t)) "Plant"
stretch(harpoon.rt)^(xi(t)) "Sensor"
stretch(harpoon.rt)^(y_k "sampled at" s_k) \
stretch(harpoon.rt)^("sensor-to-controller delay" tau^"sc ")
stretch(harpoon.rt)^(xi^c (t)) "Controller"
stretch(harpoon.rt)^("computational delay" tau^"c")
stretch(harpoon.rt)^(v^c (t))
dots
$

在之前的基础上增改假设：
+ 静态控制器 $u_k = kappa(x_k)$，采样信号 $x_k := xi(s_k)$。
+ #[
    传感器-控制器延迟 $tau^"sc"$：

    $
    xi^c (t) = xi(s_k), quad "for" t in [s_k + tau^"sc", s_(k+1) + tau^"sc")
    $

    即在采样区间后移一段延迟后，控制器得到的（ZOH）延迟状态输入 $xi^c (t)$ 才等于 $xi(s_k)$。
]
+ 计算延迟 $tau^"c"$：控制器输出的 $v^c (t) = kappa(xi^c (t-tau^"c"))$，即 $kappa(xi^c (t))$ 延迟 $tau^"c"$ 后的信号。
+ 控制器-执行器延迟 $tau^"ca"$：$v(t) = v^c (t-tau^"ca")$。
+ 常数总延迟：$tau = tau^"sc" + tau^"c" + tau^"ca"$。
+ 小延迟：$0 <= tau <= h$，即总延迟不超过采样间隔。

照样把系统写出来：

$
dot(xi) (t) &= A xi(t) + B v(t), &quad& forall t in RR^+ \
v(t) &= u_k, && "for" t in [s_k + tau, s_(k+1) + tau)
$

解还是一样，只是 $v(t)$ 有了延迟：

$
xi(t) = e^(A(t-t_0)) xi(t_0) + integral_(t_0)^t e^(A(t-s)) B v(s) dif s
$

同样，考虑 $t in [s_k, s_(k+1))$，代入 $t_0 = s_k, t = s_(k+1)$ 有：

$
x_(k+1) := xi(s_(k+1)) = e^(A h) xi(s_k) + integral_(s_k)^(s_(k+1)) e^(A(s_(k+1)-s)) B #Cbl($v(s)$) dif s
$

和前面不一样的是其中 $v(t)$ 由于延迟不再恒定为 $u_k$：

$
v(t) = cases(
    u_(k-1)"," &quad& t in [s_k, s_k + tau),
    u_k"," && t in [s_k + tau, s_(k+1))
)
$

于是拆分一下代入有：

$
x_(k+1) &= e^(A h) xi(s_k) + integral_(s_k)^(s_k + tau) e^(A(s_(k+1)-s)) B u_(k-1) dif s + integral_(s_k + tau)^(s_(k+1)) e^(A(s_(k+1)-s)) B u_k dif s \
&= e^(A h) x_k + (integral_(h-tau)^h e^(A s) B dif s) u_(k-1) + (integral_0^(h-tau) e^(A s) B dif s) u_k
$

系数分别记为 $F_x (h) := e^(A h)$、$F_u (h, tau) := integral_(h-tau)^h e^(A s) B dif s$ 和 $G_1 (h, tau) := integral_0^(h-tau) e^(A s) B dif s$ 即有：

$
x_(k+1) = F_x (h) x_k + F_u (h, tau) u_(k-1) + G_1 (h, tau) u_k
$

将上时刻输入 $u_(k-1)$ 纳入增广状态 $x_k^e = [x_k^top, u_(k-1)^top]^top$ 中则有：

$
x_(k+1)^e = F(h, tau) x_k^e + G(h, tau) u_k
$

其中：

$
F(h, tau) := mat(delim: "[", F_x (h), F_u (h, tau); 0, 0), quad
G(h, tau) := mat(delim: "[", G_1 (h, tau); I)
$

以上推导须满足假设 $0 <= tau <= h$，若考虑更大的常数延迟如 $tau in [(d-1) h, d h], d > 1$ 则需要将更多的输入纳入状态中。

接下来分析，依旧考虑状态反馈控制器 $u_k = -K x_k^e = -macron(K) x_k - K_u u_(k-1)$，其中 $K = [macron(K), K_u]$，得闭环系统：

$
x_(k+1)^e = (F(h, tau) - G(h, tau) K) x_k^e
$

实际上当且仅当 $K_u = 0$ 时这个控制器才是真的 "static"，即不依赖过去历史且无内部状态。依旧分析 $rho((F(h, tau) - G(h, tau) K))$ 判断其 GES 性，依赖采样间隔和延迟。

== NCS With Time-Varying Sampling Intervals

接下来考虑可能变化的采样间隔，假设：
+ 采样时间现在为采样间隔的累积 $s_(k+1) = sum_(i=0)^k h_i$，即有 $h_k = s_(k+1) - s_k$。
+ 假设采样间隔虽然可能变化，但只能从集合 $h_k in cal(H) := {h^1, dots, h^l}, forall k$ 的有限个候选值中选。
+ 其他假设差不多，不考虑延迟。

和常规 Sampled-Data 系统一样，只是 $F$ 和 $G$ 将与 $h_k$ 相关，即：

$
x_(k+1) = F(h_k) x_k + G(h_k) u_k
$

其中 $F(h_k) = e^(A h_k), quad G(h_k) = integral_0^(h_k) e^(A s) B dif s$，为时变系统。

考虑控制器 $u_k = -macron(K) x_k$，若采样间隔是时变的，这就是一个 switched system 或者说时变系统，其稳定性*不能*靠检查每一步 $k$ 的闭环矩阵 $A_k^"cl" = F(h_k) - G(h_k) macron(K)$ 的特征值来判断，因为实际上决定 $x_k = A_(k-1)^"cl" dots A_0^"cl" x_0$ 稳定性的是这一串矩阵乘积的特征值，不同的交替乘积结果可能不同。

还有就是如果采样间隔 uncertain but constant，即虽然不确定但一旦运行起来就固定采用某个 $h = h^i in cal(H)$ 了，那显然就可以分别判别 $cal(H)$ 中各采样间隔对应闭环矩阵的特征值来判定稳定性了。

== Stability Analysis of Non-Linear Discrete-Time Systems

考虑如下离散时间非线性时变系统：

$
x(k+1) = f(x(k), k), quad forall k
$

状态 $x_k in RR^n$，并且原点是一个不动点，即有：$0 = f(0, k), forall k$。

TODO (uniform) stability；uniformity

TODO (uniform) attractivity；uniformity

TODO local (uniform) asymptotic stability；global

TODO $cal(K)$, $cal(K)_infinity$, $cal(K L)$ functions

TODO (global) uniform asymptotic stability

TODO global exponential stability

TODO LTI systems stability

== Linear Matrix Inequalities (LMIs) for Stability Analysis

TODO positive/negative (semi-)definiteness

TODO 一些关于对称和正定的等价条件和性质

TODO 前述时变采样间隔系统稳定性的线性矩阵不等式（LMI）判据
$
exists P = P^T succ 0, Q = Q^top succ 0, forall h in cal(H), \
(F(h) - G(h) macron(K))^top P (F(h) - G(h) macron(K)) - P prec.eq -Q
$

== A Tutorial on Linear Matrix Inequalities (LMIs) and Semi-Definite Programming (SDP)

== NCS With Packet Losses

== NCS With Communication Constraints (Protocols)

== Multi-Hop Control Systems

== Overview of Advanced (LMI) Techniques for NCS With Multiple Uncertainties

== Event-Triggered Control
