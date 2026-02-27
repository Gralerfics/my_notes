#import "../generic.typ": *

#import "@preview/cetz:0.4.2"
#import "@preview/fletcher:0.5.8" as fletcher: diagram, node, edge

= Introduction & Mathematical Background

== Discrete-time Dynamical Systems

考虑*离散时间状态空间模型*如下：

$
x(k+1) &= f(x(k), u(k)) #h(5em) &x^+ &= f(x, u) \
y(k) &= h(x(k), u(k)) &y &= h(x, u) \
x(0) &= x_0 &x(0) &= x_0
$

右侧是仿连续时间状态空间模型的一种简化写法（状态空间模型描述的只是相邻状态之间的变化方式），只是符号，不重要。

连续时间情况下，状态空间模型就是一个微分方程，我们可以通过初值 $x_0$ 和已知的输入信号 $u(t)$ 求解出状态函数 $x(t)$ 的表达式（即求解一个 initial-value problem）。离散时间状态下模型是一个差分方程，可以通过初值 $x_0$ 和已知输入序列 $bold(u)_k$ 求解 $k$ 时刻的状态 $x(k)$，称为*状态解*（state solution），记为 $phi (k; x_0, bold(u)_k)$。这种记法方便显式地表达求解它所需的初值、输入信息。

对于线性离散时间动力学系统，

TODO
