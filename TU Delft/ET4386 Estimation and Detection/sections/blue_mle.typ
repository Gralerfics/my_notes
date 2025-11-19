#import "../generic.typ": *

#import "@preview/cetz:0.4.2"
#import "@preview/fletcher:0.5.8" as fletcher: diagram, node, edge

= Practical Estimators

== Maximum Likelihood Estimation (MLE)

== Best Linear Unbiased Estimation (BLUE)

== Blackboard bkup

$
x = A + omega, quad omega ~ cal(N)(0, y_2 I)
$

$
ln 1(?): ln 1/(pi A)^(1/2) - (x-A)^2/A
$

$
ln l(?): ln 1/(pi A)^(N/2) - (sum (x[n] - A)^2)/A
$

$
s(x, A)(?): -N/(2A) + [(sum x^2[n])/A^2-N] attach(=, t: ?) T(A) [g(x) - A], quad g(x) ~ hat(A)_"MVUB"(?)
$
