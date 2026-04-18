#import "../generic.typ": *

#import "@preview/cetz:0.4.2"
#import "@preview/fletcher:0.5.8" as fletcher: diagram, node, edge

= Introduction & Mathematical Background

== Model Predictive Control (MPC)

MPC 本身作为一个反馈控制器存在，通过*预测*未来系统行为计算当前*最优*控制输入。

不同于 LQR 在一开始求解最优控制策略后按部就班地执行，MPC 为了动态地抵抗噪声、干扰、模型精度不足的影响，会在执行一小步最优策略后根据所处的新状态重新计算最优策略，循环往复。顺便，由于这种一步一步的离散执行策略，通常只能处理离散时间系统或者将连续时间系统离散化。

具体地，这被称为*滚动域策略*（receding horizon policy），在每一个状态处预测未来行为并求解最优输入序列 $bold(u)^0$，但只取最优输入的第一步作为控制输入。进一步，向前预测无穷远时间的行为是难以做到且意义不大的，所以通常考虑长度为 $N$ 的域（滚动窗口）。

#blockquote([
    *关于每一步所得最优输入序列的关系*：

    TODO // 1、干扰、噪声、模型不准；2、有限的 N。
])

TODO 一些日常例子的建模：Lec 1.1 P17。

== Discrete-time Dynamical Systems

考虑*离散时间状态空间模型*如下：

$
x(k+1) &= f(x(k), u(k)) #h(5em) &x^+ &= f(x, u) \
y(k) &= h(x(k), u(k)) &y &= h(x, u) \
x(0) &= x_0 &x(0) &= x_0
$

右侧是仿连续时间状态空间模型的一种简化写法（状态空间模型描述的只是相邻状态之间的变化方式），只是符号，不重要。

连续时间情况下，状态空间模型就是一个微分方程，我们可以通过初值 $x_0$ 和已知的输入信号 $u(t)$ 求解出状态函数 $x(t)$ 的表达式（即求解一个 initial-value problem）。离散时间状态下模型是一个差分方程，可以通过初值 $x_0$ 和已知输入序列 $bold(u)_k$ 求解 $k$ 时刻的状态 $x(k)$，称为*状态解*（state solution），记为 $phi (k; x_0, bold(u)_k)$。这种记法方便显式地表达求解它所需的初值、输入信息，后面也会很常用。

对于线性离散时间动力学系统，

TODO
