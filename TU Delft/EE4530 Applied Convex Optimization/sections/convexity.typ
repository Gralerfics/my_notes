#import "../generic.typ": *

#import "@preview/cetz:0.4.2"
#import "@preview/fletcher:0.5.8" as fletcher: diagram, node, edge

= Convex Sets

#Cre("TODO 大多图都先不放，参考 Slides。")

== Affine and Convex Sets

=== Lines and Line Segments

*直线（line）*：穿过 $x_1, x_2$ 的直线是所有满足

$
x = theta x_1 + (1 - theta) x_2, theta in RR
$

的点的集合。

*线段（line segment）*：对应直线上 $0 <= theta <= 1$ 的部分。

=== Affine Sets

*仿射集（affine set）*：穿过集合中任意两点的直线都仍在该集合中，例如线性方程组的解（解的线性组合仍然是解）；反过来，任何仿射集也可以用一组线性方程来表达。

#blockquote([
    即集合中任意元素的线性组合都仍然在该集合中。
])

=== Convex Sets

*凸集（convex set）*：集合中任意两点连成的线段都仍在该集合中，即有

$
x_1, x_2 in C, quad 0 <= theta <= 1 quad => quad theta x_1 + (1 - theta) x_2 in C
$

*凸组合（convex combination）*：$x_1, dots, x_k$ 的凸组合是所有满足

$
x = theta_1 x_1 + theta_2 x_2 + dots + theta_k x_k, quad "with" theta_1 + dots + theta_k = 1, theta_i >= 0
$

的点的集合。

*凸包（convex hull）*：集合 $S$ 中所有点构成的所有凸组合的并集，记为 $bold("conv") S$。

=== Cones

*锥组合（conic combination）*：$x_1$ 和 $x_2$ 构成的锥组合是所有满足

$
x = theta_1 x_1 + theta_2 x_2, quad "with" theta_1 >= 0, theta_2 >= 0
$

的点的集合。

#blockquote([
    原点 $O$ 以及 $x_1$ 和 $x_2$ 构成的三角区域。
])

*凸锥（convex cone）*：集合中所有点构成的锥组合的并集。

*对偶锥（dual cone）*：锥 $K$ 的对偶锥是 $K^* = {y | y^T x >= 0, forall x in K}$。

例如（前三个都是自对偶即 self-dual 的）：
+ $(RR_+^n)^* = RR_+^n$；
+ $(SS_+^n)^* = SS_+^n$；
+ $({(x, t) | norm(x)_2 <= t})^* = {(x, t) | norm(x)_2 <= t}$；
+ $({(x, t) | norm(x)_1 <= t})^* = {(x, t) | norm(x)_infinity <= t}$。

== Some Important Examples

=== Hyperplanes and Halfspaces

*超平面（hyperplane）*：形如 ${x | a^T x = b}" "(a != 0)$ 的集合。

*半空间（halfspace）*：形如 ${x | a^T x <= b}" "(a != 0)$ 的集合。

#blockquote([
    其中 $a$ 为法向量（normal vector）。
])

超平面是 affine 并且 convex 的；半空间是 convex 的。

=== Euclidean Balls and Ellipsoids

*欧几里得球（Euclidean ball）*：中心在 $x_c$，半径为 $r$，满足

$
B(x_c, r) = {x | norm(x - x_c)_2 <= r} = {x_c + r u | norm(u)_2 <= r}
$

*椭球（ellipsoid）*：形如

$
{x | (x - x_c)^T P^(-1) (x - x_c) <= 1}, quad "with" P in SS_(++)^n
$

的集合，或 ${x_c + A u | norm(u)_2 <= 1}$，$A$ 为非奇异方阵。

=== Norm Balls and Norm Cones

*范数（norm）*：函数 $norm(dot)$ 满足
- $norm(x) >= 0 "and" norm(x) = 0 "iff" x = 0$；
- $norm(t x) = abs(t) norm(x), t in RR$；
- $norm(x + y) <= norm(x) + norm(y)$。

加下标表示具体的范数函数，例如 $norm(dot)_2$ 等。

*范数球（norm ball）*：形如 ${x | norm(x - x_c) <= r}$ 的集合。

*范数锥（norm cone）*：形如 ${(x, t) | norm(x) <= t}$ 的集合。

欧几里得范数锥也称二阶锥（second-order cone）。

范数球和范数锥都是 convex 的。

=== Polyhedra

*多面体（Polyhedra）*：有限个半空间和超平面的交集，也对应了有限多个线性不等式和等式组的解，即

$
A x prec.eq b, quad C x = d, quad "where" A in RR^(m times n), C in RR^(p times n)
$

=== Positive Semidefinite Cone

半正定 $n times n$ 对称矩阵集 $S_(+)^n$ 是一个凸锥。

== Operations That Preserve Convexity

要证明一个集合是否是凸集，可从定义证明，也可证明该集合是由一些简单凸集（如超平面）通过一些*不改变凸性的操作*变换得来的。

常见的不改变集合凸性的操作如下。

=== Intersection

任何凸集的交集都是凸集。

=== Affine Functions

凸集通过仿射变换 $f: RR^n -> RR^m$（$f(x) = A x + b "with" A in RR^(m times n), b in RR^m$）后的集合仍然是凸集。

例如：缩放、平移、投影。顺便，仿射变换的反函数也是 affine 的。

=== Linear-fractional and Perspective Functions 

#Cre("TODO 书上有但课件没展开")

= Convex Functions

== Basic Properties and Examples

=== Definition

*凸函数（convex function）*：若 $bold("dom") f$ 是 convex 的并且有：

$
f(theta x + (1 - theta) y) <= theta f(x) + (1 - theta) f(y), quad forall x, y in bold("dom") f, 0 <= theta <= 1
$

若满足：

$
f(theta x + (1 - theta) y) < theta f(x) + (1 - theta) f(y), quad forall x, y in bold("dom") f, x != y, 0 < theta < 1
$

则为严格凸函数，即 strictly convex 的。

*凹函数（concave function）*：若 $-f$ 是凸函数则 $f$ 是凹函数。

则函数 $f: RR^n -> RR$ 是 convex 的。

=== Extended-value Extensions

Extended-value extension $tilde(f)$ 是：

$
tilde(f)(x) = f(x), quad x in bold("dom") f; quad tilde(f)(x) = infinity, quad x in.not bold("dom") f
$

#blockquote([
    就是延拓函数的定义域，把原定义域之外的值都设为无穷。作用是延拓后可以省略定义里关于定义域为凸集的条件。
])

=== First-order Conditions

若 $f$ 的定义域是开集，且其梯度

$
nabla f(x) = ((partial f(x))/(partial x_1), (partial f(x))/(partial x_2), dots, (partial f(x))/(partial x_n))
$

在定义域中每点都存在，则 $f$ 是*可微（differentiable）*的。

对于定义域为凸集的可微函数 $f$，其为 convex 的一阶（充要）条件为：

$
f(y) >= f(x) + nabla f(x)^T (y - x), quad forall x, y in bold("dom") f
$

#blockquote([
    就是函数的一阶近似超平面永远在函数的下面。
])

=== Second-order Conditions

=== Examples

== Operations That Preserve Convexity

== The Conjugate Function

== Quasiconvex Functions

== Log-concave and Log-convex Functions
