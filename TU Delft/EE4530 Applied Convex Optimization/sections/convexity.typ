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

#underline[仿射集（affine set）]：穿过集合中任意两点的直线都仍在该集合中。

即集合中任意元素的线性组合都仍然在该集合中，例如线性方程组的解（解的线性组合仍然是解）；反过来，任何仿射集也可以用一组线性方程来表达。

=== Convex Sets

#underline[凸集（convex set）]：集合中任意两点连成的线段都仍在该集合中，即有

$
forall x_1, x_2 in C, quad 0 <= theta <= 1 quad => quad theta x_1 + (1 - theta) x_2 in C
$

#underline[凸组合（convex combination）]：$x_1, dots, x_k$ 的凸组合是所有满足

$
x = theta_1 x_1 + theta_2 x_2 + dots + theta_k x_k, quad "with" theta_1 + dots + theta_k = 1, theta_i >= 0
$

的点的集合。

#underline[凸包（convex hull）]：集合 $S$ 中所有点构成的所有凸组合的并集，记为 $bold("conv") S$。

=== Cones

#underline[锥（cone）]：集合中任意元素的缩放都仍在集合中，即有

$
forall x in C, alpha >= 0 quad => quad alpha x in C
$

#underline[锥组合（conic combination）]：$x_1$ 和 $x_2$ 构成的锥组合是所有满足

$
x = theta_1 x_1 + theta_2 x_2, quad "with" theta_1 >= 0, theta_2 >= 0
$

的点的集合。

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
    *法向量与图形理解*：

    我们常说 $a$ 为法向量（normal vector），那么如何直观理解？具体地，显然超平面的约束 $a^T x = b$ 可化为如下形式：

    $
    a^T (x - b/(a^T a) a) = 0, quad "or" quad a^T (x - b/norm(a)_2^2 a)
    $

    其中，$a/norm(a)_2^2$ 是 $a$ 方向上的单位向量，所以 $b/norm(a)_2^2 a$ 就是 $a$ 方向上长度为 $b$ 的向量。
    
    于是超平面的约束可以看成所有 $x$ 在去掉这样一段偏移后，与 $a$ 的方向均垂直。这就会构成一个平面，其上所有线都与 $a$ 垂直，其自然就是法向量。同时，去除的那段偏移，代表平面是从原点朝着 $a$ 方向移动了距离 $b$ 得到的（由于 $b$ 的符号可能有正负，所以不能说明法向量与原点存在固定的方向关系）。

    此外，半空间的约束也可写作：

    $
    a^T (x - b/(a^T a) a) <= 0
    $

    说明所包含的半片空间是不同于法向量 $a$ 朝向的另一侧的那半片空间（如扯住 $a$ 拉开的半片帷幕）。
])

超平面是 affine 并且 convex 的；半空间是 convex 的。

=== Euclidean Balls and Ellipsoids

#underline[欧几里得球（Euclidean ball）]：中心在 $x_c$，半径为 $r$，满足

$
B(x_c, r) = {x | norm(x - x_c)_2 <= r} = {x_c + r u | norm(u)_2 <= 1}
$

#underline[椭球（ellipsoid）]：形如

$
{x | (x - x_c)^T P^(-1) (x - x_c) <= 1}, quad "with" P in SS_(++)^n
$

的集合，或 ${x_c + A u | norm(u)_2 <= 1}$，$A$ 为非奇异方阵，两种表示的关系是 $P = A A^T$，椭球的第 $i$ 条半主轴长度为 $sqrt(lambda_i (P))$。

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

// 证明，用半正定对应二次型恒非负证明：
// #image("/assets/image-35.png")
// 凸性也可以用一系列半空间交集证明：
// #image("/assets/image-36.png")
// #image("/assets/image-37.png") // 顺便，为什么是半空间

== Operations That Preserve Convexity

要证明一个集合是否是凸集，可从定义证明，也可证明该集合是由一些简单凸集（如超平面）通过一些*不改变凸性的操作*变换得来的。

常见的不改变集合凸性的操作如下。

=== Intersection

任何凸集的交集都是凸集。

=== Affine Functions

凸集通过仿射变换 $f: RR^n -> RR^m$（$f(x) = A x + b "with" A in RR^(m times n), b in RR^m$）后的集合仍然是凸集。

#Cre("TODO") 例如：缩放、平移、投影。顺便，仿射变换的反函数也是 affine 的。

#Cre("TODO") 线性不等式组的解集 ${x | x_1 A_1 + dots + x_m A_m prec.curly.eq B, A_i in SS^p, B in SS^p}$ 也是凸的。#underline[因为]约束里的不等式表示 $x_1 A_1 + dots + x_m A_m$ 这一串属于半负定对称矩阵，是凸集，而解集 $x$ 是它的一个逆仿射变换，也是凸的。

#Cre("TODO") 双曲锥（hyperbolic cone）即 ${x | x^T P x <= (c^T x)^2, c^T x >= 0, P in SS_+^n}$ 也是凸的。#underline[因为]首先我们可以把它写成类似二阶范数锥的形式：

$
{x | norm(P^(1\/2) x)_2 <= c^T x, c^T x >= 0, P in SS_+^n}
$

然后定义一个仿射函数：

$
f(x) = vec(P^(1\/2) x, c^T x)
$

将其与正儿八经的二阶锥建立联系：

$
{(x, t) | norm(x)_2 <= t, t >= 0}
$

已知二阶锥是凸集，所以用 $f^(-1)$ 就可以得到双曲锥是凸集。

=== Linear-fractional and Perspective Functions 

#Cre("TODO 书上有但课件没展开")

= Convex Functions

== Basic Properties and Examples

=== Definition

#underline[凸函数（convex function）]：若 $bold("dom") f$ 是 convex 的并且有：

$
f(theta x + (1 - theta) y) <= theta f(x) + (1 - theta) f(y), quad forall x, y in bold("dom") f, 0 <= theta <= 1
$

则函数 $f: RR^n -> RR$ 是 convex 的，注意加上定义域为凸集总共是两条约束。若它满足的是更严格的：

$
f(theta x + (1 - theta) y) < theta f(x) + (1 - theta) f(y), quad forall x, y in bold("dom") f, x != y, 0 < theta < 1
$

则为严格凸函数，即 strictly convex 的。

#underline[凹函数（concave function）]：若 $-f$ 是凸函数，则 $f$ 是凹函数。

=== Epigraph and Sublevel Set

#underline[上图（epigraph）]：函数 $f: RR^n -> RR$ 的上图为

$
bold("epi") f = {(x, t) in RR^(n+1) | x in bold("dom") f, f(x) <= t}
$

*函数上图为凸集和函数为凸函数互为充分必要条件*，可以理解为函数上方的整个空间是凸的，那么函数就是凸函数。

#underline[下图（hypograph）]：函数 $f: RR^n -> RR$ 的下图为

$
bold("epi") f = {(x, t) in RR^(n+1) | x in bold("dom") f, f(x) >= t}
$

#underline[下水平集（sublevel set）]：函数 $f: RR^n -> RR$ 的 $alpha$-下水平集为

$
C_alpha = {x in bold("dom") f | f(x) <= alpha}
$

即对应的函数值不高于某个平面的取值的集合。*凸函数的下水平集都是凸的，但反之不一定成立*。

=== Extended-value Extensions

Extended-value extension $tilde(f)$ 是：

$
tilde(f)(x) = f(x), quad x in bold("dom") f; quad tilde(f)(x) = infinity, quad x in.not bold("dom") f
$

#blockquote([
    就是延拓函数的定义域，把原定义域之外的值都设为无穷。作用是延拓后可以省略定义里关于定义域为凸集的条件。
])

=== First-order Conditions

若 $f$ 的定义域是开集（open），且其梯度

$
nabla f(x) = ((partial f(x))/(partial x_1), (partial f(x))/(partial x_2), dots, (partial f(x))/(partial x_n))
$

在定义域中每点都存在，则 $f$ 是*可微*（differentiable）的。

对于定义域为凸集的可微函数 $f$，其为 convex 的一阶（充要）条件为：

$
f(y) >= f(x) + nabla f(x)^T (y - x), quad forall x, y in bold("dom") f
$

#blockquote([
    就是函数的一阶近似超平面永远在函数的下面。
])

=== Second-order Conditions

若 $f$ 的定义域是开集，且其 Hessian 矩阵 $nabla^2 f(x) in SS^n$，

$
nabla^2 f(x)_(i j) = (partial^2 f(x))/(partial x_i partial x_j), quad i, j = 1, dots, n,
$

对任意 $x in bold("dom") f$ 存在，则 $f$ 是*二阶可微*（twice differentiable）的。

对于定义域为凸集的二阶可微函数 $f$，其为 convex 的二阶（充要）条件为：

$
nabla^2 f(x) succ.curly.eq 0, quad forall x in bold("dom") f
$

若正定则为 strictly convex 的。

=== Examples

#Cre("TODO") 

// R 上

R 上的仿射既凸又凹；

指数 e^(a x) 对于任意 R 上的 a 都凸；

正底数幂，指数在 [0, 1] 是凹，[-inf, 0] U [1, inf] 是凸；abs(x)^p 则可以在 R 上，p >=1 就凸；

负熵 x log x 在 R++（正）上为凸；

对数 log x 在 R++ 上为凹；

// R 向量、矩阵

仿射函数都是既凸又凹，$f(x) = a^T x + b$ 和 $f(X) = tr(A^T B) + b$；

范数都凸，$norm(x)_p$ 对于 $p >= 1$，以及矩阵的 spectral (maximum singular value) norm $f(X) = norm(X)_2 = sigma_max (X) = sqrt(lambda_max (X^T X))$；

还有二次函数（quadratic function）$f(x) = 1/2 x^T P x + q^T x + r$，其中 $P in SS^n$，有 $nabla f(x) = P x + q$ 和 $nabla^2 f(x) = P$，如果 $P succ.curly.eq 0$ *才*凸；

最小二乘目标函数 $f(x) = norm(A x - b)_2^2$，有 $nabla f(x) = 2 A^T (A x - b)$ 和 $nabla^2 f(x) = 2 A^T A$，对于任意 $A$ 都凸；

对于二次除以一次，例如 $f(x, y) = x^2 \/ y$ 有 $nabla^2 f(x) = 2/y^3 vec(y, -x) vec(y, -x)^T succ.curly.eq 0$，$y > 0$ *时*凸；

对于指数和的对数（log-sum-exp）例如 $f(x) = log sum_(k=1)^n exp(x_k)$ 是凸的，用二阶条件和 Cauchy-Schwarz 不等式证明；

对于几何均值（geometric mean）例如 $f(x) = (product_(k=1)^n x_k)^(1/n)$ 在 $RR_(++)^n$ 上是凹的，证明类似 log-sum-exp；

// #image("/assets/image-38.png")
// #image("/assets/image-39.png")
// #image("/assets/image-40.png")
// #image("/assets/image-41.png")

== Operations That Preserve Convexity

一些不改变凸性的操作：
+ 非负加权和；
+ 复合仿射函数，若 $f$ 凸则 $f(A x + b)$ 也凸；
+ 最大值，$f(x) = max{f_1 (x), dots, f_m (x)}$ 凸如果每个 $f_i (x)$ 都凸；
+ 上确界（supremum），若 $f(x, y)$ 对每个 $y in cal(A)$ 都在 $x$ 上凸，则 $g(x) = sup_(y in cal(A)) f(x, y)$ 凸（？？？）；
+ 复合（composition）标量函数，$f(x) = h(g(x))$，首先要求套外面的 $h$ 凸，然后 $tilde(h)$ 不减且 $g$ 凸，或者 $tilde(h)$ 不增且 $g$ 凹，可得 $f$ 凸；
+ 复合向量函数，改成 $f(x) = h(g(x)) = h(g_1 (x), g_2 (x), dots, g_k (x))$，要求变成针对每一个 $g_i$ 以及 $tilde(h)$ 的每一项，例如 log-sum-exp；
+ 下确界（infimum），若 $f(x, y)$ 在 $(x, y)$ 上凸且 $C$ 是一个凸集，则 $g(x) = inf_(y in C) f(x, y)$ 凸；
+ 共轭函数（conjugate function）：$f^*(y) = sup_(x in bold("dom") f) (y^T x - f(x))$ 是凸的，即便 $f$ 不凸也是（可以用前面的 pointwise supremum 证明）；
+ ……

// TODO 好像从某处开始已经不属于这个 section 了

== Quasiconvex Functions

// 定义域和任意-下水平集都凸，则为拟凸的；加个负号就拟凹；都的话就拟线性。

// 一堆 examples

== Log-convex and Log-concave Functions

// 如果 log f 是 convex 的则 f 是 log-convex 的。
