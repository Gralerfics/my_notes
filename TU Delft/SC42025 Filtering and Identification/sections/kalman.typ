#import "../generic.typ": *

#import "@preview/cetz:0.4.2"
#import "@preview/fletcher:0.5.8" as fletcher: diagram, node, edge

= Kalman Filtering

// 状态空间模型的状态估计问题。

// 怎么解决？观测器（根据输入输出估计状态的都叫观测器）；卢恩伯格观测器，A - KC 那个。

// 如果带过程和测量噪声？发现卢恩伯格观测器无法确保渐近稳定 -> 卡尔曼滤波（也算观测器）。

// 简化一下模型（去掉 Du），阐述卡尔曼滤波问题（用 ... 求 ... 的 MVUE）；

// TODO，类比 SLS（为什么那样写（这样写使得状态包含 x(k) 和 x(k+1)，对应卡尔曼滤波问题里估计两个东西？）？那样的协方差 LLT 为什么和前面一致？）。

// 接上，两个问题，不能直接套 SLS 的解：1、需要 [x(k), x(k+1)]T 的 prior，但我们只有 x(k) 的（见卡尔曼滤波问题定义，即每轮一开始只有 hat(x)(k|k) 和 P(k|k)）；2、SLS 是去最小化两个状态组合起来的协方差，但这里我们是要单独最小化它们各自（？）。

// 不能直接套 SLS，先直接给出解（Conventional Kalman filter 的），证明放后面（看课件和书一致否，多半是课件还多个贝叶斯角度）。

// 书上后面还有 SRCF，平方根容积卡尔曼？似乎是用来处理非线性情况的数值方法，课上似乎没讲；这部分也对应前面 SLS 后 Section 4.5.4 有关于 SLS 的 square-root solution 的内容。

// TODO ... 低维状态的卡尔曼示例。

// 作业 Ex2 展示了多个 measurement 一起更新和依次更新，结果是一致的。

== Kalman Filtering Problem

$
x_(k+1) &= A x_k + B u_k + w_k, &quad &w_k ~ cal(N)(0, Q) \
y_k &=  C x_k + D u_k + v_k, & &v_k ~ cal(N)(0, R)
$

== Innovation Form Representation



== Asymptotic Observer Perspective


