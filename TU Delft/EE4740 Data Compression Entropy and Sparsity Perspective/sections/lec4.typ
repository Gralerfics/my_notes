#import "../generic.typ": *

#import "@preview/cetz:0.4.2"
#import "@preview/fletcher:0.5.8" as fletcher: diagram, node, edge

= Introduction to Compressed Sensing (CS)

我们感兴趣的原始信号是 $bold(s) in RR^N$。由线性代数它可以在不同的正交基下进行表示。设 $bold(Phi), bold(Psi) in RR^(N times N)$ 为*单位正交矩阵*（orthonormal martix），表示两组不同的正交基，根据性质有：

$
bold(Phi)^top bold(Phi) = bold(Phi) bold(Phi)^top = bold(I)
$

即 $bold(Phi)^(-1) = bold(Phi)^top$，$bold(Psi)$ 类似。将 $bold(s)$ 分别投影到这两组基上，记新坐标分别为 $bold(g)$ 和 $bold(x)$，即有：

$
bold(s) = bold(Phi) bold(g) = bold(Psi) bold(x)
$ <equ:s_with_diff_basis>

基变换不损失信息，所以 $bold(s)$、$bold(g)$ 和 $bold(x)$ 都具有同样的信息量。我们假设选取的 $bold(Psi)$ 足够巧妙，让 $bold(x) = bold(Psi)^top bold(s)$ 中存在的零尽可能多（$K$ 个）。若 $K < N$，那么信号就是 *$K$-稀疏*（$K$-sparse）的，表示实际上用 $K$ 个元素（维度）的向量就可以表达信号中的信息，而这样的基 $bold(Psi)$ 称为*稀疏基*（sparsifying basis）。

另一个 $bold(Phi)$ 可称为*测量基*（measurement basis），用来表达一个测量的过程（假设它是观测原信号 $bold(s)$ 所用设备对信号产生的变换），总之我们可以特地选取。假设这样一种情形：我们实际能得到的只有 $bold(g) = bold(Phi)^top bold(s)$ 经过一个降采样算子（subsampling operator）降采样之后的信号，即 $bold(y) = cal(P)_Omega (bold(g))$，其中 $Omega$ 是采样索引的集合；由于前面提到过对于稀疏信号有效信息的维度低于实际维度，那么我们能否从这个残缺的信号 $bold(y)$ 中还原出 $bold(g)$ 乃至 $bold(s)$？

从信息论的角度来说，如果信号足够稀疏、$Omega$ 足够合适，降采样过程中就可能没有丢失什么信息，从而允许通过合理选取 $bold(Phi)$ 还原出原信号，这种还原的算法就称为*压缩感知*（compressed sensing，CS）算法。

== Information Theoretic Perspective

我们先从信息论的角度更严谨地阐述信号中所谓 "信息" 的概念。对于数字信号 $bold(x) in RR^N$，其每个元素 $x_i$ 的可能取值来源于有限符号集 $abs(cal(X))$，并且都是从中独立、均匀采样所得。我们将其建模为一个独立同分布随机过程 $bold(X) = {X_i}_(i = 1)^N$，总共有 $abs(cal(X))^N$ 种等可能实现，现可计算其熵为：

$
H(bold(X)) &= H(X_1, X_2, dots, X_N) \
&= - sum_(bold(x) in cal(X)^N) PP(bold(x)) log_2 PP(bold(x)) \
&= - abs(cal(X))^N (1/abs(cal(X))^N log_2 1/abs(cal(X))^N) \
&= log_2 abs(cal(X))^N \
&= N log_2 abs(cal(X))
$

现在考虑 $bold(x)$ 确定是 $K$-稀疏且最简的情况（符合 @equ:s_with_diff_basis 处的说明），即 $bold(x)$ 正好有 $K$ 个非零元素和 $N-K$ 个零。首先是一些假设：（1）既然有零元素那么一定得有 $0 in cal(X)$，顺便记排除掉 $0$ 之后的符号集为 $cal(X)'$，有 $abs(cal(X)') = abs(cal(X)) - 1$；（2）信号 $bold(x)$ 确为 $K$-稀疏的；（3）$bold(X)$ 依旧是从 $cal(X)$ 中独立、均匀地进行采样；（4）$N >> K$ 且 $N >> cal(X)$。这种情况下等可能实现的种数变为 $C_N^K abs(cal(X)')^K$，其中 $C_N^K$ 为组合数 $vec(N, K)$。由此，该稀疏信号的熵为：

$
H(bold(X)) &= H(X_1, X_2, dots, X_N) \
&= - sum_(bold(x)) PP(bold(x)) log_2 PP(bold(x)) \
&= - C_N^K abs(cal(X)')^K (1/(C_N^K abs(cal(X)')^K) log_2 1/(C_N^K abs(cal(X)')^K)) \
&= log_2 C_N^K abs(cal(X)')^K \
&= log_2 C_N^K + K log_2 abs(cal(X)') \
$

对组合数放缩一下（证明略）：

$
(N/K)^K <= C_N^K <= ((N e)/K)^K quad => quad 
K log_2 N/K <= log_2 C_N^K <= K log_2 (N e)/K
$

于是有：

$
K log_2 N/K + K log_2 abs(cal(X)') <= H(bold(X)) <= K log_2 (N e)/K + K log_2 abs(cal(X)')
$

由该式，可以粗略地说*一个 $N$ 维 $K$-稀疏数字信号的熵大于 $K log_2 N/K$*。

// 因为比较接近，可以假设 $abs(cal(X)') approx abs(cal(X))$ 以统一符号。

== Measurement Model and Recovery

=== Choice of a Good Measurement Basis

记某个测量样本为 $y_m$，有：

$
y_m
&= chevron.l bold(phi.alt)_i, bold(s) chevron.r
= chevron.l bold(phi.alt)_i, bold(Psi) bold(x) chevron.r
= sum_(j=1)^N chevron.l bold(phi.alt)_i, bold(psi)_j x_j chevron.r
= sum_(j=1)^N chevron.l bold(phi.alt)_i, bold(psi)_j chevron.r x_j
, quad m = 1, 2, dots, M \
&= chevron.l bold(phi.alt)_i, bold(psi)_1 chevron.r x_1 + chevron.l bold(phi.alt)_i, bold(psi)_2 chevron.r x_2 + dots + chevron.l bold(phi.alt)_i, bold(psi)_N chevron.r x_N
$ <equ:meas_y_x_lc>

其中 $chevron.l bold(x), bold(y) chevron.r = bold(x)^top bold(y)$ 表示内积，顺便如果考虑复向量那么都得换成 Hermitian；$bold(phi.alt)_i$ 为 $bold(Phi)$ 的第 $i$ 列，即 $bold(Phi)^top$ 的第 $i$ 行，与 $bold(s)$ 内积后得到 $bold(g)$ 的某个元素，也就是 $bold(y)$ 的某个元素。总结一下就是#underline[选一个]测量基向量 $bold(phi.alt)_i$，其在稀疏基上的投影和稀疏信号 $bold(x)$ 的线性组合即测量样本 $y_m$。

我们要从 $bold(y)$ 中还原出 $bold(x)$，很自然的想法就是让 $bold(y)$ 中每个元素都 "分布式地" 同时包含 $bold(x)$ 每个元素的一部分信息，这样删去 $bold(y)$ 中的个别元素对整体的影响就相对小。@equ:meas_y_x_lc 已经将 $y_m$ 表示为了 $bold(x)$ 元素的线性组合，如果能让式中的系数 $chevron.l bold(phi.alt)_i, bold(psi)_j chevron.r$ 均匀一些，也就相当于让信息分散得均匀一些。

反过来说，如果某个测量基向量 $bold(phi.alt)_i$ 和某个稀疏基向量 $bold(psi)_j$ 太像，即它对应的测量精准地提取出了某个稀疏分量的全部信息，那么删去这个测量就将损失巨大，如同把鸡蛋装进了一个篮子里。

我们可以根据如上分析定义一种指标来衡量一个 $bold(Phi)$ 好不好，称为*一致性*（coherence），衡量的就是测量基向量和稀疏基向量两两之间的最大相似度：

$
mu (bold(Phi), bold(Psi)) = sqrt(N) max_(i,j) abs(chevron.l bold(phi.alt)_i", "bold(psi)_j chevron.r)
$

如前所述，这个相似度自然是越小越好。可以证明 $1 <= mu (bold(Phi), bold(Psi)) <= sqrt(N)$，下界为 $1$ 是由于二者都是单位正交矩阵，$bold(phi.alt)_i$ 模长为 $1$，即有 $sum_(j=1)^N abs(chevron.l bold(phi.alt)_i", "bold(psi)_j chevron.r)^2 = 1$，而满足条件的让一致性最小的情况就是让所有 $chevron.l bold(phi.alt)_i, bold(psi)_j chevron.r = plus.minus 1/sqrt(N)$，结果为 $sqrt(N) dot 1/sqrt(N) = 1$；与之相反，如果存在某个 $bold(phi.alt)_i = plus.minus bold(psi)_j$，那么就是最坏的情况，一致性为 $sqrt(N)$。

现有如下*定理*：一个在正交基 $bold(Psi)$ 上的 $K$-稀疏信号，可以以 $1 - delta$ 的成功概率用 $M >= C K mu^2 (bold(Phi), bold(Psi)) log N/delta$ 个在正交基 $bold(Phi)$ 上的测量还原出来。// 当然，实际问题中我们通常也没法知道 $bold(Psi)$，仅用于理论分析。

// 每个对未知系数信号的测量提供 $log_2 abs(cal(X))$ 比特的信息，记要还原出 $bold(x)$ 需要的最少测量个数为 $M_o$。

举例说明来看实际应用中的选择，假设 $bold(Psi) = bold(I)_4 = [bold(e)_1, dots, bold(e)_4]$，其中 $bold(e)_i$ 是仅第 $i$ 个元素为 $1$、其他都为 $0$ 的列向量。一种可行的测量基是 Hadamard 基：

$
bold(Phi)^"hd" = 1/2 mat(
    delim: "[",
    1, 1, 1, 1;
    1, -1, 1, -1;
    1, 1, -1, -1;
    1, -1, -1, 1
)
$

计算可知对其任意 $i, j$ 有 $abs(chevron.l bold(phi.alt)^"hd"_i", "bold(e)_j chevron.r) = 1/sqrt(4)$，即一致性为最小值 $1$。此外，还可以采用 DFT 基：

$
bold(Phi)^"DFT" = 1/2 mat(
    delim: "[",
    1, 1, 1, 1;
    1, -j, -1, j;
    1, -1, 1, -1;
    1, j, -1, -j
)
$

大同小异，一致性也为下界 $1$。继续以前述 Hadamard 基为例，假设信号 $bold(x)$ 是 $K = 1$ 稀疏的，又因 $bold(Psi) = bold(I)_4$，所以信号 $bold(s)$ 只可能是四种：

$
bold(s)_1 = vec(alpha, 0, 0, 0), quad
bold(s)_2 = vec(0, alpha, 0, 0), quad
bold(s)_3 = vec(0, 0, alpha, 0), quad
bold(s)_4 = vec(0, 0, 0, alpha)
$

现在假设我们取 Hadamard 基的前三列作为测量基，可得测量样本 $bold(y) = vec(y_1, y_2, y_3) = vec(chevron.l bold(phi.alt)^"hd"_1", "bold(s) chevron.r, chevron.l bold(phi.alt)^"hd"_2", "bold(s) chevron.r, chevron.l bold(phi.alt)^"hd"_3", "bold(s) chevron.r)$。前述 $bold(s)$ 的四种情况分布对应测量结果：

$
bold(y)_1 = 1/2 vec(alpha, alpha, alpha), quad
bold(y)_2 = 1/2 vec(alpha, -alpha, alpha), quad
bold(y)_3 = 1/2 vec(alpha, alpha, -alpha), quad
bold(y)_4 = 1/2 vec(alpha, -alpha, -alpha)
$

现在比如说我们用这三个基测得 $bold(y) = [-0.45, 0.45, 0.45]^top$，可以利用上述分析匹配上 $bold(y)_4$，求得 $alpha = -0.9$，即原信号 $bold(s) = [0, 0, 0, -0.9]^top$，这就是一个手动还原信号的例子。

需要注意的是上述分析都建立在 $K=1$ 的假设上，如果稀疏假设不成立，我们是无法用上述方法还原任意 $bold(s) in RR^4$ 的。

=== Measurement Model with Random Projections

前面为了联系实际含义，符号较多且杂乱，之后使用更抽象的方式定义问题：

$
bold(y) = bold(A) bold(x) + bold(v)" " stretch(arrow)^" CS algorithm " " "hat(bold(x))
$

其中 $bold(y)$ 为测量向量，$bold(x)$ 为稀疏向量，$bold(A)$ 为 CS 矩阵（包含了 $bold(Psi)$ 到 $bold(Phi)$ 再到降采样的整个过程），$bold(v)$ 为噪声。还是定义 $Omega$ 为测量基索引，比如所取 $Omega = {4, 1, 3}$，那么就有：

$
bold(y)
= vec(y_1, y_2, y_3)
= vec(delim: "[", bold(phi.alt)_4^top, bold(phi.alt)_1^top, bold(phi.alt)_3^top) bold(Psi) bold(x)
= underbrace(vec(delim: "[", bold(phi.alt)_4^top bold(Psi), bold(phi.alt)_1^top bold(Psi), bold(phi.alt)_3^top bold(Psi)), bold(A)) bold(x)
= underbrace(mat(
    delim: "[",
    0, 0, 0, 1;
    1, 0, 0, 0;
    0, 0, 1, 0
) bold(Phi)^top bold(Psi), bold(A)) bold(x)
$

在压缩感知问题中 $bold(A)$ 的列数 $N$ 是要比行数 $M$ 更多的，所以这是一组欠定（underdetermined）方程。

将这些东西全都整合进 $bold(A)^(M times N)$ 矩阵之后就可以考虑另一种构造测量的方式，即完全不管什么 $bold(Phi)$、$bold(Psi)$ 的，直接随机产生整个 $bold(A)$：各元素都独立同分布地从 ${+1, -1}$ 上的 $1/2$ 伯努利分布采样，或者都独立同分布地从 $cal(N) (0, 1)$ 上采样，当满足 $M >= K log_2 N/K$ 时，这样随机生成的 $bold(A)$ 很大概率就是一个合适的 CS 矩阵。
