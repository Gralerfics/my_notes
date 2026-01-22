#import "../generic.typ": *

#import "@preview/cetz:0.4.2"
#import "@preview/fletcher:0.5.8" as fletcher: diagram, node, edge

=== Least Squares Estimator (LSE)

// 最优标准基于统计意义上的 MSE（话说前面哪里基于这个了来着？），而 LSE 的标准直接基于模拟输出和测量输出差的平方和，无需统计假设，但可能导致并非统计最优。

// 两个例子：A + w[n]；线性模型。

==== Geometrical Interpretation of Least Squares

==== Statistical Properties of Least Squares

==== Variations of Least Squares

===== Weighted Least Squares

===== Sequential Least Squares

===== Constrained Least Squares

===== Non-linear Least Squares

=== Summary

// Practical Estimators 总的 summary：各自的 motivation，所需要知道的已知量
