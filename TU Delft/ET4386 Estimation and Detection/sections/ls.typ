#import "../generic.typ": *

#import "@preview/cetz:0.4.2"
#import "@preview/fletcher:0.5.8" as fletcher: diagram, node, edge

=== Least Squares Estimator (LSE)

#Cre("TODO") 最优标准基于统计意义上的 MSE（话说前面哪里基于这个了来着？），而 LSE 的标准直接基于模拟输出和测量输出差的平方和，无需统计假设，但可能导致并非统计最优。

#Cre("TODO") 两个例子：A + w[n]；线性模型。

==== Geometrical Interpretation of Least Squares

#Cre("TODO") 经典图，关注一下符号。

==== Statistical Properties of Least Squares

#Cre("TODO") 偏不偏？MSE 是多少？什么时候统计最优？slides 没写，问一下确认一下。

==== Variations of Least Squares

===== Weighted Least Squares

#Cre("TODO") 加权最小二乘，slides 没写但 summary 里提了。

===== Sequential Least Squares

#Cre("TODO") 推导加入一个新测量数据后如何更新估计。

===== Constrained Least Squares

#Cre("TODO") 带约束的最小二乘优化问题，有点神秘，貌似大概是 QR 分解然后只管有值的部分，没明白怎么解耦了。

===== Non-linear Least Squares

#Cre("TODO") 局部线性化，迭代求解，F&I 里都有，这里关注一下符号的不同。

=== Summary

// Practical Estimators 总的 summary：各自的 motivation，所需要知道的已知量
