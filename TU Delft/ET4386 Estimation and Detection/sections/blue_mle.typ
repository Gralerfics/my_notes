#import "../generic.typ": *

#import "@preview/cetz:0.4.2"
#import "@preview/fletcher:0.5.8" as fletcher: diagram, node, edge

= Practical Estimators

由前面的讨论我们有了一个基本的认识：想要得到 MVU 估计量需要有关于真实模型（概率分布）的知识，这在实际问题中通常是无法实现的；而且即便知道了概率分布，也不一定能找到方法导出 MVUE。

所以在实际情况下我们会采用一些 sub-optimal 的估计量，例如后面要介绍的 MLE、BLUE 和 LS 等。在一些特定情况下，这些 sub-optimal 估计量也可能就是 MVU 估计量，或者方差收敛到 CRLB。

== Maximum Likelihood Estimation (MLE)



== Best Linear Unbiased Estimation (BLUE)
