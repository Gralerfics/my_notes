#import "../generic.typ": *

#import "@preview/cetz:0.4.2"
#import "@preview/fletcher:0.5.8" as fletcher: diagram, node, edge
#import "@preview/algo:0.3.6": algo

= Least-Squares Estimates

== Notations and Symbols

按课件的偏好，我们给出统一的估计问题的模型。首先是*信号模型*：

$
y = F theta + L epsilon, quad epsilon ~ cal(N)(0, I_m)
$

其中 $theta$ 是待估计的参数，$L epsilon$ 是噪声项。规范一下维度，$F in RR^(m times n)" "(m>=n)$，$y in RR^m$， $theta in RR^n$，$epsilon in RR^m$，而 $L in RR^(m times m)$ 是非奇异方阵。*更清晰一些*：

$
[y]_(m times 1) = [F]_(m times n) [theta]_(n times 1) + [L]_(m times m) [epsilon]_(m times 1), quad epsilon ~ cal(N)([0]_(m times 1), I_m)
$

我们令 $epsilon$ 是服从*单位高斯分布*的（均值为 $0$，协方差为 $I$），再用一个分开的 $L$ 来允许其拓展到噪声协方差不为单位矩阵的其他情况。

我们也可以将该模型等效一下，用一个具有非单位协方差的新的 $epsilon$ 来替代 $L epsilon$ 整体。具体地，计算 $L epsilon$ 的均值：

$
E[L epsilon] = L dot E[epsilon] = L dot 0 = 0
$

和协方差：

$
E[(L epsilon - E[L epsilon])(L epsilon - E[L epsilon])^T] &= E[(L epsilon)(L epsilon)^T] \
&=E[L epsilon epsilon^T L^T] \
&=L E[epsilon epsilon^T] L^T \
&=L L^T
$

顺便，根据高斯分布的性质，$L epsilon$ 仍然服从高斯分布（$L$ 是线性变换）。于是*前述模型等效于*把 $L epsilon$ 合并成一个均值为零、协方差为 $L L^T$ 的新 $epsilon$：

$
y = F theta + epsilon, quad epsilon ~ cal(N)(0,L L^T)
$ <equ:lse_note_equv_model>

== The Linear Least-Squares Problem

我们先只考虑单位高斯分布的 $epsilon$，即基于上述信号模型中 $L = I_m$ 的情况（可省去），提出最基础的线性*最小二乘问题*：

$
min_theta epsilon^T epsilon = min_theta norm(epsilon)^2_2, quad s.t." "y = F theta + epsilon, quad epsilon ~ cal(N)(0, I_m)
$

顺便，"s.t." 指 "subject to"。

// #blockquote([
//     #Cre("TODO")1. 这里我们考虑 deterministic 的 $theta$，...
//     #Cre("TODO")2. $epsilon$ 含义的变化，...
// ])

对该问题的*最小方差无偏估计量*（minimum variance unbiased estimator）为：

#emphasis_equbox([
$
hat(theta) = (F^T F)^(-1) F^T y
$
])

其协方差为（计算证明详见书第 109 页，Lemma 4.4）：

#emphasis_equbox([
$
E[(theta - hat(theta))(theta - hat(theta))^T] = (F^T F)^(-1)
$
])

#Cre("TODO")证明这个解是 minimum variance unbiased estimator，书 P111。

== The Weighted Linear Least-Squares Problem

#underline[接下来考虑当 $L$ 不为 $I_m$ 时的情况]。新的问题为：

$
min_theta epsilon^T epsilon = min_theta norm(epsilon)^2_2, quad s.t." "#Cpu($y = F theta + L epsilon$), quad epsilon ~ cal(N)(0, I_m)
$ <equ:lse_wlls_basic>

将条件两侧同乘 $L^(-1)$ 可得到：

$
min_theta epsilon^T epsilon = min_theta norm(epsilon)^2_2, quad s.t." "#Cpu($L^(-1) y = L^(-1) F theta + epsilon$), quad epsilon ~ cal(N)(0, I_m)
$

该形式和上节中相同，所以将上节解中的 $y$ 和 $F$ 换为这里的 $L^(-1) y$ 和 $L^(-1) F$ 即可得到该问题的最小方差无偏估计量：

$
hat(theta) &= (F^T (L^T)^(-1) L^(-1) F)^(-1) F^T (L^T)^(-1) L^(-1) y \
&= (F^T (L L^T)^(-1) F)^(-1) F^T (L L^T)^(-1) y
$

若记 $W = (L L^T)^(-1)$，即为：

#emphasis_equbox([
$
hat(theta) = (F^T W F)^(-1) F^T W y
$
])

类似地，其协方差为：

#emphasis_equbox([
$
E[(theta - hat(theta))(theta - hat(theta))^T] = (F^T W F)^(-1)
$
])

#blockquote([
    *关于 "加权" 的注释*：

    由前 @equ:lse_note_equv_model 的思路将 $L epsilon$ 合起来。令 $epsilon' = L epsilon$，即有 $epsilon = L^(-1) epsilon'$，则 @equ:lse_wlls_basic 可化为：

    $
    min_theta (epsilon')^T W epsilon', quad s.t." "y = F theta + epsilon', quad epsilon' ~ cal(N)(0, L L^T)
    $

    不想看撇号，把变量换回来，问题等价于：

    $
    min_theta norm(epsilon)_W^2 := min_theta #Cbl($epsilon^T W epsilon$), quad s.t." "y = F theta + epsilon, quad epsilon ~ cal(N)(0, #Cbl($L L^T$))
    $

    在范数右下角标记 $W$ 来表示这种形式的二次型。从这个角度来看，这里的 $W$ 表达了某种 "权重"（#underline("W")eight）的意味，其决定了 $epsilon$ 不同元素（或交叉项）在优化问题中的相对影响力大小，故该问题类型称为加权线性最小二乘问题。
    
    具体一些，$W = (L L^T)^(-1)$ 其实就是等价后 $epsilon$ 的协方差的逆，表达的含义即#underline[ “哪个误差项的方差大，哪个就越不可靠，给它的权重就越小”]。
])

== Nonlinear Least-Squares Problem

顺便考虑一下非线性的情况：

$
min_theta norm(epsilon)_2^2, quad y = f(theta) + L epsilon, quad epsilon ~ cal(N)(0, I_m)
$

方法自然就是在 $hat(theta)$ 处将 $f(theta)$ 线性化。由泰勒展开，在 $theta = hat(theta)$ 附近，$f(theta)$ 可一阶近似为：

$
f(theta) approx f(hat(theta)) + lr((dif f(theta))/(dif theta)|)_(theta = hat(theta)) + dots
$

代入 $y = f(theta) + L epsilon$ 即有：

$
underbrace(y - f(hat(theta)), e) approx underbrace(lr((dif f(theta))/(dif theta)|)_(theta = hat(theta)), F) underbrace((theta - hat(theta)), Delta theta) + L epsilon
$

求解算法自然又是*迭代*，即选取一个初始的对结果的猜测 $hat(theta)^((0))$，在每一点都求解上式的最小二乘问题，根据 $Delta theta$ 的方向不断更新猜测 $hat(theta)^((i))$：

#algo(
    title: "Nonlinear least-squares algorithm",
    parameters: ($y$, $f(theta)$, $L$, $hat(theta)^((0))$, text("stopping_criterion"))
)[
    #import "@preview/algo:0.3.6": i, d, comment, code
    $i <- 0$ \
    while not stopping_criterion: #i \
        $F_i <- lr((dif f(theta))/(dif theta)|)_(theta = hat(theta)^((i)))$ \
        $e_i <- y - f(hat(theta)^((i)))$ \
        $Delta theta <-$ 最小二乘问题 $e_i = F_i Delta theta + L epsilon$ 的解 \
        $hat(theta)^((i+1)) <- hat(theta)^((i)) + Delta theta$ \
        $i <- i + 1$ #comment("并检查收敛性")
]

// #Cre("TODO") nls.m 在哪呢？

#pagebreak()
== The Stochastic Linear Least-Squares (SLS) Problem

#underline[按惯例先阐明接下来要讨论的问题和前面有什么不一样]。

*首先*，前面（加权）线性最小二乘问题中我们都认为 $theta$ 是 deterministic 的，即是确定的数值，而我们做的事情是去估计这个数值；*现在*，我们把 $theta$ 视为一个#underline[随机变量]，我们去估计的是最终这个随机变量的统计特征（最好就是直接能估计出分布函数，包含一切统计特征）。

那么#underline[到这里其实只是思维上的区别]，因为如果我们的目的还是得到一个最好的、单一的数值，那么就算我们把 $theta$ 看成随机变量，最后也是会取它的均值作为结果。

#underline[那么将 $theta$ 视为随机变量的好处在哪]？是在估计的过程之中，我们可以保持更多信息，例如猜测、倾向、可能性等。

#blockquote([
    *以卡尔曼滤波的图像解释为例说明*：
    
    Kalman 滤波器假设要估计的参数是一个随机变量，服从正态分布，在估计的过程中保留了其均值和协方差矩阵，每一步更新就像一个椭球在空间里转移、变形。
    
    其中，协方差矩阵就可以用来保存更多关于 "对当前估计各方面的自信程度" 的信息，在每一步的更新中提供更多参考。

    而如果我们不保留这些统计信息，只把参数当作 deterministic 的数值去估计，那么除了最小二乘也只能最小二乘了，参考前面非线性最小二乘算法。
])

#underline[回到问题上来]，在 Stochastic Least-Squares（SLS）这里，我们也将 $theta$ 看作随机变量，并且关注估计过程的最开始：我们对 $theta$ 这个参数有一些最初的假设和期待，即先验的（priori）知识。

例如，我们相信参数 $theta$ 的均值为 $mu_theta$、协方差矩阵为 $P_theta succ.curly.eq 0$，*这就是该问题相较前面基础的最小二乘问题多出来的已知条件*，最终导致 @equ:lse_sls_solution 中的解也大变样了。后面我们会提到，#underline[当先验知识无效时，这里的 SLS 问题将退化为前面的加权最小二乘形式]，具体见 @sec:lse_sls_alter_form_sol 末尾。

=== Bayesian Probability, Prior and Posterior <sec:lse_sls_bayesian_exp>

#Cgy([
这里稍微写一点关于贝叶斯概率、先验、后验的例子等，进一步说明一下这里加入先验条件对问题的影响，以及先验信息如何在许多问题中被忽略或默认假设，*或可跳过*。

#Cre("TODO") 连续抛两次硬币为正，第三次为正的概率；不同先验分布（均匀/Beta(1,1)、确信公平/delta、……）对结果的影响……

#Cre("TODO") 先验信息如何在许多问题中被忽略

#Cre("TODO") 前面加权最小二乘关于误差的协方差假设也是一种先验（？）
])

#pagebreak()
=== Problem Statement and The Solution

#underline[接下来正式定义 SLS 问题]：给定关于参数的先验期望 $mu_theta$ 和先验协方差 $P_theta succ.curly.eq 0$，及观察模型：

$
y = F theta + L epsilon, quad epsilon ~ cal(N)(0, I_m)
$

其中 $F$ 和 $L$ 是 deterministic 的，$L$ 为可逆方阵。且有正交法则（#Cre("TODO")？）：

$
E[(theta - mu_theta) epsilon^T] = 0
$

求综合这些信息之后的最优估计量 $hat(theta)$ 及其协方差 $E[(theta - hat(theta))(theta - hat(theta))^T]$。

我们*先直接给出该问题的解*为：

#emphasis_equbox([
$
hat(theta) = K y + (I_n - K F) mu_theta, quad "where" K = P_theta F^T (F P_theta F^T + W^(-1))^(-1)
$ <equ:lse_sls_solution>
])

顺便，把 $K$ 代进去写就是这样：

$
hat(theta) = P_theta F^T (F P_theta F^T + W^(-1))^(-1) y + (I_n - P_theta F^T (F P_theta F^T + W^(-1))^(-1) F) mu_theta
$

其协方差为：

#emphasis_equbox([
$
E[(theta - hat(theta))(theta - hat(theta))^T] &= P_theta - P_theta F^T (F P_theta F^T + W^(-1))^(-1) F P_theta \
&= #Cpu($(I_n - K F) P_theta$)
$ <equ:lse_sls_sol_corv_common>
])

该解的证明见 @sec:lse_sls_deri，接下来先推导一个它的等价形式，并借助其介绍一些理解内容。

#pagebreak()
=== An Alternative Form of the Solution and the Minimal-Covariance <sec:lse_sls_alter_form_sol>

若 $P_theta succ 0$（即正定），我们可以将 @equ:lse_sls_sol_corv_common 中的结果协方差*改写*为：

#emphasis_equbox([
$
E[(theta - hat(theta))(theta - hat(theta))^T] = #Cpu($(P_theta^(-1) + F^T W F)^(-1)$) =: P_(theta,"post")
$ <equ:lse_sls_sol_corv_post>
])

可以叫它*后验协方差* $P_"post"$，即结合了先验知识以及数据后得到的结果协方差。

#blockquote([
    *证明*：

    使用矩阵逆引理（Woodbury 公式）：

    $
    #Cbl($($)A + B C D#Cbl($)^(-1)$) = A^(-1) - A^(-1) B (C^(-1) + D A^(-1) B)^(-1) D A^(-1)
    $

    将后验协方差的表达式展开：

    $
    P_(theta,"post") &= #Cbl($($)P_theta^(-1) + F^T W F#Cbl($)^(-1)$) \
    &= P_theta - P_theta F^T (F P_theta F^T + W^(-1))^(-1) F P_theta \
    &= (I_n - K F) P_theta
    $

    即得到 @equ:lse_sls_sol_corv_common 的形式。
])

用这个后验协方差，我们可以#underline[将 @equ:lse_sls_solution 改写成一种更 "贝叶斯" 的形式]。先从 $K$ 的定义出发证明 $K = P_(theta,"post") F^T W$：

$
K &= P_theta F^T (F P_theta F^T + W^(-1))^(-1) \
K (F P_theta F^T + W^(-1)) &= P_theta F^T \
K F P_theta F^T + K W^(-1) &= P_theta F^T \
K F P_theta F^T W + K &= P_theta F^T W \
K &= P_theta F^T W - K F P_theta F^T W \
K &= (I_n - K F) P_theta F^T W \
#Cpu($K &= P_(theta,"post") F^T W$)
$ <equ:lse_sls_post_k>

然后用上述结论对 @equ:lse_sls_solution 作#underline[等价变换]：

$
hat(theta) &= K y + (I_n - K F) mu_theta \
&= underbrace(K, P_(theta,"post") F^T W) y + underbrace((I_n - K F) #Cbl($P_theta$), P_(theta,"post")) #Cbl($P_theta^(-1)$) mu_theta \
&= #Cpu($P_(theta,"post")$) F^T W y + #Cpu($P_(theta,"post")$) P_theta^(-1) mu_theta \
&= P_(theta,"post") (P_theta^(-1) mu_theta + F^T W y)
$

即最终解的等价形式如下：

#emphasis_equbox([
$
hat(theta) = underbrace((P_theta^(-1) + F^T W F)^(-1), P_(theta,"post")) (P_theta^(-1) mu_theta + F^T W y)
$ <equ:lse_sls_sol_bayes_like_form>
])

#blockquote([
    *后验协方差的构成*：

    观察这个后验协方差的形式，它实际上是先验估计 $theta ~ (mu_theta, P_theta)$ 和加权最小二乘最优估计量 $theta ~ ((F^T W F)^(-1) F^T W y, (F^T W F)^(-1))$ 这两种估计量的协方差结合后的结果：各自逆一下加起来，再逆回去。
    
    顺便，Slides 中都用 $~ cal(N)(dot)$ 来表示分布，但实际上这里从头到尾都#underline[没有作正态分布假设]，只是在使用均值和协方差两个统计特征，所以为了避免混淆，我这里就不写 $cal(N)$，而是只写二元组来表示具有相应均值和协方差的分布。

    #underline[该观点即 Lec2 Slides 第 30 页中的第一点。]
])

#blockquote([
    *随机最小二乘问题和加权最小二乘问题*：

    根据节首的分析，我们认定先验信息的加入是导致 SLS 问题的解发生变化的主要原因。那么显然，#underline[如果提供的先验信息没用]，或者说知不知道没区别，#underline[那么问题就应该退化为加权最小二乘问题]。我们用刚得到的解的等价形式 @equ:lse_sls_sol_bayes_like_form 可以很容易地验证这一点。

    首先，"无用的" 先验信息指非常 "平坦" 的先验分布，即对于参数 $theta$ 的初始假设中，取各种值的概率都相同，也就是说*先验协方差极大*，趋于无穷，即 $P_theta -> infinity$。

    于是，将 $P_theta -> infinity$ 即 $P_theta^(-1) -> 0$ 代入 @equ:lse_sls_sol_bayes_like_form：

    $
    hat(theta)_* &= (#Cbl($P_theta^(-1)$) + F^T W F)^(-1) (#Cbl($P_theta^(-1)$) #Cpu($mu_theta$) + F^T W y) \
    &= (F^T W F)^(-1) F^T W y
    $

    显然，#underline[这个就是前面加权最小二乘问题的解]。顺便我们可以发现，这种情况下先验均值 $mu_theta$ 也失去意义了（全都均匀），数学上则体现为被一起消去了。

    #underline[这部分即 Lec2 Slides 第 30 页中的第二点。]
])

#pagebreak()
=== The Recursive Least-Squares Algorithm

基于上述 SLS 问题，我们可以提出一种优化的迭代算法，即 RLS：

#algo(
    title: "Recursive least-squares (RLS) algorithm",
    parameters: ($y_(1:N)$, $F_(1:N)$, $W_(1:N)$, $hat(theta)_0$, $P_0$)
)[
    #import "@preview/algo:0.3.6": i, d, comment, code
    Initially supposing $theta ~ (hat(theta)_0, P_0)$ \
    for $k <- 1 "to" N$: #i \
        load new data $(y_k, F_k, W_k)$ \
        $K_k <- (P_(k-1)^(-1) + F_k^T W_k F_k)^(-1) F_k^T W_k$ #comment([@equ:lse_sls_post_k]) \
        $hat(theta)_k <- (I - K_k F_k) hat(theta)_(k-1) + K_k y_k$ #comment([@equ:lse_sls_solution]) \
        $P_k^(-1) <- P_(k-1)^(-1) + F_k^T W_k F_k$ #comment([@equ:lse_sls_sol_corv_post])
]

给定初始先验后，用每一个新数据对估计进行更新；#underline[对于每一条新数据，都将之前的阶段性估计结果视作当前的 “先验”]，用求解 SLS 问题的方式更新均值和协方差估计。

// #Cre("TODO") RLS_demo.m 在哪呢？

=== Derivation of Stochastic Least-Squares <sec:lse_sls_deri>

接下来，我们从最小方差无偏估计量的定义，以及贝叶斯估计两个角度分别*推导*以得到 @equ:lse_sls_solution 中的结果。

==== Unbiased and Minimum Variance <sec:lse_sls_deri_mvu>

#Cre("TODO")

1. Lec2 Slides 第 33 页起；
2. 书第 114 页起；
3. video_1 和 https://brightspace.tudelft.nl/d2l/le/content/775429/viewContent/4432122/View。

==== Bayesian <sec:lse_sls_deri_bayes>

#Cre("TODO")

1. Lec2 Slides 第 39 页起；
2. video_2 和 https://brightspace.tudelft.nl/d2l/le/content/775429/viewContent/4432124/View。

#Cre("TODO")

Lec3 Slides 第 3 页（From Bayesian derivations to least squares problems）？
