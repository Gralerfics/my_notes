#import "../generic.typ": *

#import "@preview/cetz:0.4.2"
#import "@preview/fletcher:0.5.8" as fletcher: diagram, node, edge

= Convex Sets

#Cre("TODO 大多图都先不放，参考 Slides。")

== Affine and Convex Sets

=== Lines and Line Segments

#underline[直线（line）]：穿过 $x_1, x_2$ 的直线是所有满足

$
x = theta x_1 + (1 - theta) x_2, theta in RR
$

的点的集合。

#underline[线段（line segment）]：对应直线上 $0 <= theta <= 1$ 的部分。

=== Affine Sets

#underline[仿射集（affine set）]：穿过集合中任意两点的直线都仍在该集合中，例如线性方程组的解（解的线性组合仍然是解）；反过来，任何仿射集也可以用一组线性方程来表达。

#blockquote([
    即集合中任意元素的线性组合都仍然在该集合中。
])

=== Convex Sets

#underline[凸集（convex set）]：集合中任意两点连成的线段都仍在该集合中，即有

$
x_1, x_2 in C, quad 0 <= theta <= 1 quad => quad theta x_1 + (1 - theta) x_2 in C
$

#underline[凸组合（convex combination）]：$x_1, dots, x_k$ 的凸组合是所有满足

$
x = theta_1 x_1 + theta_2 x_2 + dots + theta_k x_k, quad "with" theta_1 + dots + theta_k = 1, theta_i >= 0
$

的点的集合。

#underline[凸包（convex hull）]：集合 $S$ 中所有点构成的所有凸组合的并集，记为 $bold("conv") S$。

=== Cones

#underline[锥组合（conic combination）]：$x_1$ 和 $x_2$ 构成的锥组合是所有满足

$
x = theta_1 x_1 + theta_2 x_2, quad "with" theta_1 >= 0, theta_2 >= 0
$

的点的集合。

#blockquote([
    原点 $O$ 以及 $x_1$ 和 $x_2$ 构成的三角区域。
])

#underline[凸锥（convex cone）]：集合中所有点构成的锥组合的并集。

#underline[对偶锥（dual cone）]：锥 $K$ 的对偶锥是 $K^* = {y | y^T x >= 0, forall x in K}$。

例如（前三个都是自对偶即 self-dual 的）：
+ $(RR_+^n)^* = RR_+^n$；
+ $(SS_+^n)^* = SS_+^n$；
+ $({(x, t) | norm(x)_2 <= t})^* = {(x, t) | norm(x)_2 <= t}$；
+ $({(x, t) | norm(x)_1 <= t})^* = {(x, t) | norm(x)_infinity <= t}$。

#blockquote([
    *对第四个例子的证明*：

    记 $K = {(x, t) | norm(x)_1 <= t}$，其对偶锥定义为：
    
    $
    K^* = {(y, s) | y^T x + s t >= 0, forall (x, t) in K}
    $
    
    注意这里具体化了关于内积的定义。
    
    要证明 $L = {(x, t) | norm(x)_infinity <= t}$ 是 $K$ 的对偶锥即 $K^*$，即证明集合相等，需要分别证明 $L subset.eq K^*$ 和 $K^* subset.eq L$。

    #underline[先证明] $L subset.eq K^*$。要证明一个集合是另一个的子集，即要证明该集合中任意元素都同时存在于另一个集合中。我们取 $forall (y, s) in L$，有 $norm(y)_infinity <= s$（这里用 $(y, s)$ 是考虑到对偶锥定义中还要逐个验证 $(x, t)$，这样会顺一点，实际上只是符号）。
    
    然后我们检查 $(y, s)$ 是否满足 $K^*$ 的条件，满足即可说明该元素也属于 $K^*$。由 Hölder 不等式可知 $abs(y^T x) <= norm(y)_infinity norm(x)_1$，由前我们还有 $norm(y)_infinity <= s$，于是对于 $forall (x, t) in K$，即 $norm(x)_1 <= t$ 有：

    $
    y^T x + s t >= -norm(y)_infinity norm(x)_1 + s t >= -s t + s t = 0
    $

    即得 $L subset.eq K^*$。
    
    #underline[接下来证明] $K^* subset.eq L$。那么取 $forall (y, s) in K^*$，满足 $y^T x + s t >= 0, forall (x, t) in K$，看看 $(y, s)$ 是否也在 $L$ 当中，即要证明 $norm(y)_infinity <= s$。还是利用 Hölder 不等式，直接：

    $
    y^T x + s t >= -norm(y)_infinity norm(x)_1 + s t >= (s - norm(y)_infinity) t >= 0
    $

    又因为 $t >= norm(x)_1 >= 0$ 的，于是有 $s - norm(y)_infinity >= 0$，即 $norm(y)_infinity <= s$，得证 $K^* subset.eq L$。
    
    二者结合#underline[最终证明] $L = K^*$。

    #underline[关于 Hölder 不等式]：若函数 $f(x)$ 和 $g(x)$ 在 $[a, b]$ 上连续，且 $p, q > 0$，$1/p + 1/q = 1$，则有：

    $
    integral_a^b abs(f(x)g(x)) dif x <= (integral_a^b abs(f(x))^p dif x)^(1/p) (integral_a^b abs(g(x))^q dif x)^(1/q)
    $
])

== Some Important Examples

=== Hyperplanes and Halfspaces

#underline[超平面（hyperplane）]：形如 ${x | a^T x = b}" "(a != 0)$ 的集合。

#underline[半空间（halfspace）]：形如 ${x | a^T x <= b}" "(a != 0)$ 的集合。

#blockquote([
    其中 $a$ 为法向量（normal vector）。
])

超平面是 affine 并且 convex 的；半空间是 convex 的。

=== Euclidean Balls and Ellipsoids

#underline[欧几里得球（Euclidean ball）]：中心在 $x_c$，半径为 $r$，满足

$
B(x_c, r) = {x | norm(x - x_c)_2 <= r} = {x_c + r u | norm(u)_2 <= r}
$

#underline[椭球（ellipsoid）]：形如

$
{x | (x - x_c)^T P^(-1) (x - x_c) <= 1}, quad "with" P in SS_(++)^n
$

的集合，或 ${x_c + A u | norm(u)_2 <= 1}$，$A$ 为非奇异方阵。

=== Norm Balls and Norm Cones

#underline[范数（norm）]：函数 $norm(dot)$ 满足
- $norm(x) >= 0 "and" norm(x) = 0 "iff" x = 0$；
- $norm(t x) = abs(t) norm(x), t in RR$；
- $norm(x + y) <= norm(x) + norm(y)$。

加下标表示具体的范数函数，例如 $norm(dot)_2$ 等。

#underline[范数球（norm ball）]：形如 ${x | norm(x - x_c) <= r}$ 的集合。

#underline[范数锥（norm cone）]：形如 ${(x, t) | norm(x) <= t}$ 的集合。

欧几里得范数锥也称二阶锥（second-order cone）。

范数球和范数锥都是 convex 的。

=== Polyhedra

#underline[多面体（Polyhedra）]：有限个半空间和超平面的交集，也对应了有限多个线性不等式和等式组的解，即

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

#underline[凸函数（convex function）]：若 $bold("dom") f$ 是 convex 的并且有：

$
f(theta x + (1 - theta) y) <= theta f(x) + (1 - theta) f(y), quad forall x, y in bold("dom") f, 0 <= theta <= 1
$

若满足：

$
f(theta x + (1 - theta) y) < theta f(x) + (1 - theta) f(y), quad forall x, y in bold("dom") f, x != y, 0 < theta < 1
$

则为严格凸函数，即 strictly convex 的。

#underline[凹函数（concave function）]：若 $-f$ 是凸函数则 $f$ 是凹函数。

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
