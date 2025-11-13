#import "../generic.typ": *

#import "@preview/cetz:0.4.2"
#import "@preview/fletcher:0.5.8" as fletcher: diagram, node, edge

= Fundamentals

（TODO）Book Verhaegen: Study Chapters 1, 2(2.1-2.6), 3, 4(4.1-4.3)

（TODO）Verhaegen: Sections 4.5.1 – 4.5.3 (SLS, unbiased minimum variance derivation)；Lecture notes: Sections 1–3；
    Additional background: Manipulating the multivariate Gaussian density, Thomas B. Sch ̈on and Fredrik Lindsten, Division of Automatic Control, Link ̈oping University, Sweden, January 2011.

== Linear Algebra

信号模型、最小二乘问题和解

矩阵运算

奇异值分解（代入前面最小二乘得到一些结论）

秩、基本子空间，和解的关系

（半）正定

QR分解

一般估计问题（信号模型、加权最小二乘和解）

== Probabilistic Theory and Random Process

随机变量，均值方差相关，及其样本估计，以及性质

（多维）高斯分布，记号，σ，PDF

（……………………………………）

随机过程，（多过程）WSS

白噪声性质，记号

数值计算均值方差，ergodic，noiseProperties.m（？）

（……………………………………）下面这部分应该放下个section了。

统计学估计参数theta，频率学派和贝叶斯学派，最小二乘的例子

非线性的情况，非线性最小二乘求解算法，nls.m（？）

P24，TODO，RLS、SLS、……

// == Symbols and Notations
