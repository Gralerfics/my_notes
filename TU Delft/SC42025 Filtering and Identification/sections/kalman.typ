#import "../generic.typ": *

#import "@preview/cetz:0.4.2"
#import "@preview/fletcher:0.5.8" as fletcher: diagram, node, edge

= Kalman Filtering

// 状态空间模型的状态估计问题。

// 怎么解决？观测器（根据输入输出估计状态的都叫观测器）；卢恩伯格观测器，A - KC 那个。

// 如果带过程和测量噪声？发现卢恩伯格观测器无法确保渐近稳定 -> 卡尔曼滤波（也算观测器）。

// 简化一下模型（去掉 Du），阐述卡尔曼滤波问题（用 ... 求 ... 的 MVUE）；

// TODO，类比 SLS（为什么那样写（这样写使得状态包含 x(k) 和 x(k+1)，对应卡尔曼滤波问题里估计两个东西？）？那样的协方差 LLT 为什么和前面一致？）。

// 接上，两个问题，不能直接套 SLS 的解：1、需要 [x(k), x(k+1)]T 的 prior，但我们只有 x(k) 的（见卡尔曼滤波问题定义，即每轮一开始只有 hat(x)(k|k) 和 P(k|k)）；2、SLS 是去最小化两个状态组合起来的协方差，但这里我们是要单独最小化它们各自（？）。

// 不能直接套 SLS，先直接给出解（Conventional Kalman filter 的），证明放后面（看课件和书一致否，多半是课件还多个贝叶斯角度）。

// 书上后面还有 SRCF，平方根容积卡尔曼？似乎是用来处理非线性情况的数值方法，课上似乎没讲；这部分也对应前面 SLS 后 Section 4.5.4 有关于 SLS 的 square-root solution 的内容。

// TODO ... 低维状态的卡尔曼示例。

// 作业 Ex2 展示了多个 measurement 一起更新和依次更新，结果是一致的。

== State Estimation and Asymptotic Observers

继 @sec:lse_sls_rls 一段注释中的讨论之后，我们正式考虑#underline[状态估计问题]；如果将待估计状态的物理意义设定为 “真实波形”，多余部分视为 “噪声”，则也可以称其为一类#underline[滤波问题]。如那段讨论中所述，我们可以大体将#underline[静态参数估计问题]通过添加一条表达状态不变的状态转移方程并入#underline[状态估计问题]这一大类中去，虽然不一定有必要。

本节我们先讨论控制理论中提到过的#underline[观测器]（observer），分析其作用和局限，再引出 Kalman 滤波问题和 Kalman 滤波器。在状态估计问题中，我们需要根据输入 ${u_k}$ 和输出 ${y_k}$ 序列估计系统的状态序列 ${x_k}$，在控制理论中完成这一任务的系统即被称为观测器。以最基本的状态空间模型为例：

$
x_(k+1) &= A x_k + B u_k \
y_k &=  C x_k + D u_k
$

我们可以#underline[先建立一个与原系统具有相同输入的复制版本的系统]：

$
hat(x)_(k+1) &= A hat(x)_k + B u_k \
hat(y)_k &=  C hat(x)_k + D u_k
$

若初始状态 $x_0$ 一致，则对于模型相同、输入相同的两个系统，在理想的情况下行为也将相同；若初始状态不一致，但系统本身渐近稳定，最终二者也将收敛到相同状态。具体地，两个系统的状态误差用 $e_k = hat(x)_k - x_k$ 表示，用 $hat(x)_(k+1) &= A hat(x)_k + B u_k$ 减去 $x_(k+1) &= A x_k + B u_k$ 即得到误差动力学：

$
e_(k+1) = A e_k
$

若此时 $A$ 是渐近稳定（asymptotically stable）的，则误差将随时间趋于零。接下来的问题是，#underline[我们无法控制这个 “观测器” 估计的状态收敛到真实状态的速度]。若初始状态不同，两个系统状态可能要等快稳定了才会收敛到一起。如果我们用这种跟不上瞬态响应的状态估计来设计控制器，效果想必不会太好。

别忘了我们还知道输出序列，可以考虑对两个系统输出的误差施加增益：

$
hat(x)_(k+1) &= A hat(x)_k + B u_k + K (y_k - hat(y)_k) \
&= A hat(x)_k + B u_k + K (y_k - C hat(x)_k - D u_k) \
&= (A - K C) hat(x)_k + B u_k + K y_k - K D u_k \
// &= (A - K C) hat(x)_k + B u_k + K y_k - K D u_k
$

时刻注意观测器作为一个子系统，$u_k$ 和 $y_k$ 对它来说是已知的输入。仍然两式相减，得到新的误差动力学为：

$
hat(x)_(k+1) - x_(k+1) &= [(A - K C) hat(x)_k + B u_k + K y_k - K D u_k] - [A x_k + B u_k] \
hat(x)_(k+1) - x_(k+1) &= (A - K C) hat(x)_k + B u_k + K (C x_k + D u_k) - K D u_k - A x_k - B u_k \
hat(x)_(k+1) - x_(k+1) &= (A - K C) hat(x)_k - (A - K C) x_k \
e_(k+1) &= (A - K C) e_k
$

类似地，若 $A - K C$ 是渐近稳定的，则观测器的状态估计可以随时间收敛到真实状态。与前面不同的是，#underline[只要原系统 $(A, C)$ 是能观的（observable），我们就能通过设计 $K$ 任意配置 $A - K C$ 的极点位置，从而设计观测器的收敛速度、稳定性等]，由此我们得到了经典的 #underline[Luenberger 观测器]。

接下来我们考虑带噪声的情况：

$
x_(k+1) &= A x_k + B u_k + w_k, &quad &w_k &~ &cal(N)(0, Q) \
y_k &=  C x_k + D u_k + v_k, & &v_k &~ &cal(N)(0, R)
$

此时我们继续应用 Luenberger 观测器，得到的误差动力学是这样的：

$
e_(k+1) = (A - K C) e_k - w_k + K v_k
$

#underline[此时即便 $A - K C$ 被设计为渐近稳定的，由于这两个噪声项的存在，误差最终也无法收敛到零]。

当然，由于噪声的随机性时时刻刻存在，我们本就不可能精确地消除噪声。之后要提出的 Kalman 滤波器从概念上来说也可以视为一款观测器，它也不是为了完全消除噪声，我们#underline[所能做的只有让这个误差 “尽可能地小”]。怎么样才算 “小”？类比最小方差无偏估计，我们依然希望得到的状态估计#underline[在统计意义上期望收敛到零，而方差尽可能地小]，显然这样的最优性是和噪声模型挂钩的。在下面的 Kalman 滤波问题中，对噪声的建模是必要的，且会影响最优解即 Kalman 滤波器的形式。

此外，*Kalman 滤波器*的引入不止是为了*更好地利用噪声模型*，它还可以更系统性地通过计算时变的增益，*处理时变系统的状态估计问题*，而 Luenberger 观测器在这方面比较受限。

== Kalman Filtering Problem <sec:km_kfp>

我们先声明一下之后会用到的符号，状态下标的竖线左侧是状态所处的时刻，竖线右侧是已知信息所达到的时刻。

例如，$hat(x)_(k|k)$ 指的是已知到 $k$ 为止所有的输入和输出序列时对状态 $x_k$ 的估计，所以这其实#underline[就是我们综合到时刻 $k$ 为止所有信息后对 $x_k$ 的最终估计 $hat(x)_k$]。再例如，$hat(x)_(k|k-1)$ 则是已知到 $k-1$ 为止所有输入和输出序列时对状态 $x_k$ 的估计，可以由 $hat(x)_(k-1|k-1)$ 再加一步基于模型的预测得到它，因为条件都只到 $k-1$ 时刻，没有新的时刻 $k$ 的输出信息，只能利用模型预测。

从贝叶斯估计的角度我们可以写得更形式化、更具体一点，例如：

$
hat(x)_(k|k) &= arg max_(x_k) p(x_k|u_(1:k), y_(1:k)) \
hat(x)_(k|k-1) &= arg max_(x_k) p(x_k|u_(1:k-1), y_(1:k-1))
$

类似地，还有表示估计量协方差的 $P_(k|k)$ 等。

这里还有一个问题，#underline[为什么要引入 $hat(x)_(k|k-1)$ 这些中间估计量]而不直接使用 $hat(x)_k$ 代表对 $x_k$ 的估计？这将在 @sec:km_kf 的一处注释中解释。

接下来我们#underline[正式规定 Kalman 滤波问题]。给定信号模型如下：

$
x_(k+1) &= A_k x_k + B_k u_k + w_k, &quad &w_k &~ &cal(N)(0, Q_k) \
y_k &= C_k x_k + v_k, & &v_k &~ &cal(N)(0, R_k)
$

其中，过程噪声（process noise）$w_k$ 和测量噪声（measurement noise）$v_k$ 建模为具有如下联合协方差矩阵的零均值白噪声序列：

$
E[mat(delim: "[", v_k; w_k) mat(delim: "[", v_l^T, w_l^T)] = mat(delim: "[", R_k, S_k^T; S_k, Q_k) Delta(k - l) succ.curly.eq 0
$

其中，$R_k succ 0$，$Delta(k - l)$ 是单位脉冲函数，这里表达了不同时刻之间的噪声不相关。

记在时刻 $k-1$ 对 $x_k$ 的最优估计为 $hat(x)_(k|k-1)$，与噪声不相关，满足如下性质：

$
&E[hat(x)_(k|k-1)] = E[x_k] \
&E[(x_k - hat(x)_(k|k-1))(x_k - hat(x)_(k|k-1))^T] = P_(k|k-1) succ.curly.eq 0
$

接下来在每个时刻 $k$，我们需要利用已知的 $u_k$、$y_k$ 等，以及前一轮递推得到的 $hat(x)_(k|k-1)$（对当前时刻来说是已知的常数，且满足上述性质），寻找针对 $x_k$ 和 $x_(k+1)$ 的线性估计量：

$
mat(delim: "[", hat(x)_(k|k); hat(x)_(k+1|k)) = M mat(delim: "[", y_k; -B_k u_k; hat(x)_(k|k-1))
$

其中 $M in RR^(2 n times (l + 2 n))$，使得该估计量（两项）是最小方差无偏的。具体地，即满足：

$
E[hat(x)_(k|k)] &= E[x_k] \
E[hat(x)_(k+1|k)] &= E[x_(k+1)]
$

同时最小化 $E[(x_k - hat(x)_(k|k))(x_k - hat(x)_(k|k))^T]$ 和 $E[(x_(k+1) - hat(x)_(k+1|k))(x_(k+1) - hat(x)_(k+1|k))^T]$，即 $P_(k|k)$ 和 $P_(k+1|k)$。

#blockquote([
    *关于模型中被简化的 $D_k u_k$ 项*：

    完整的测量模型应当是 $y_k = C_k x_k + D_k u_k + v_k, & &v_k &~ &cal(N)(0, R)$，书中包括上述问题陈述中为简化，省略了输入相关项。不过由于 $u_k$、$y_k$ 和 $D_k$ 是已知的信息，我们可以直接移项得到：

    $
    tilde(y)_k := y_k - D_k u_k = C_k x_k + v_k, & &v_k &~ &cal(N)(0, R_k)
    $

    所以若考虑这一项，提前计算 $tilde(y)_k = y_k - D_k u_k$，代入替换问题中和推导中的所有 $y_k$ 即可。
])

== Kalman Filter <sec:km_kf>

=== Solution to the Kalman Filtering Problem

我们先直接给出问题的解，即 Kalman 滤波器的形式。

记 $K_k = #Cpu($P_(k|k-1) C_k^T (C_k P_(k|k-1) C_k^T + R_k)^(-1)$)$，称为 Kalman 增益。

时刻 $k$ 的最优状态估计为：

$
hat(x)_(k|k) &= #Cpu($P_(k|k-1) C_k^T (C_k P_(k|k-1) C_k^T + R_k)^(-1)$) #Cbl($y_k$) \
&quad + [I_n - #Cpu($P_(k|k-1) C_k^T (C_k P_(k|k-1) C_k^T + R_k)^(-1)$) C_k] #Cbl($hat(x)_(k|k-1)$) \
&= K_k #Cbl($y_k$) + (I_n - K_k C_k) #Cbl($hat(x)_(k|k-1)$) \
&= hat(x)_(k|k-1) + K_k (y_k - C_k hat(x)_(k|k-1))
$ <eq:km_kf_xkk>

第一行将这个线性估计量写成了最朴素的关于 $y_k$ 和 $hat(x)_(k|k-1)$ 的#underline[线性组合]的形式，和问题定义中一致；第二行用 #underline[Kalman 增益]化简了式子，是常见的形式；最后一行将其写成了类似 #underline[“更新”] 的形式，其中有关 $y_k - C_k hat(x)_(k|k-1)$ 将在 @sec:km_innovation_form 中进一步阐述。

该估计量的协方差为：

$
P_(k|k) &= P_(k|k-1) - #Cpu($P_(k|k-1) C_k^T (C_k P_(k|k-1) C_k^T + R_k)^(-1)$) C_k P_(k|k-1) \
&= (I_n - K_k C_k) P_(k|k-1)
$ <eq:km_kf_Pkk>

接下来，在时刻 $k$ 对 $k+1$ 的最优估计为：

$
hat(x)_(k+1|k) &= (A_k P_(k|k-1) C_k^T + S_k) (C_k P_(k|k-1) C_k^T + R_k)^(-1) #Cbl($y_k$) \
&quad + B_k #Cbl($u_k$) \
&quad + [A_k - (A_k P_(k|k-1) C_k^T + S_k) (C_k P_(k|k-1) C_k^T + R_k)^(-1) C_k] #Cbl($hat(x)_(k|k-1)$) \

&= A_k underbrace([K_k #Cbl($y_k$) + (I_n - K_k C_k) #Cbl($hat(x)_(k|k-1)$)], hat(x)_(k|k)) + B_k #Cbl($u_k$) \
&quad + S_k (C_k P_(k|k-1) C_k^T + R_k)^(-1) (y_k - C_k hat(x)_(k|k-1)) \

&= #Cgr($A_k hat(x)_(k|k) + B_k u_k$) + S_k (C_k P_(k|k-1) C_k^T + R_k)^(-1) (y_k - C_k hat(x)_(k|k-1))
$ <eq:km_kf_xk1k>

类似地，第一行是朴素的线性组合形式；第二行和第三行将其整理了一下（尽量把 $S_k$ 单独提出去），提取出了和前面 $hat(x)_(k|k)$ 一致的部分，这样 $A_k hat(x)_(k|k) + B_k u_k$ 这部分就是状态转移模型，即表达模型随时间演化的部分，后面再加上与 $S_k$ 相关的项。

它的协方差为：

$
P_(k+1|k) &= A_k P_(k|k-1) A_k^T + Q_k \
&quad - (A_k P_(k|k-1) C_k^T + S_k) (C_k P_(k|k-1) C_k^T + R_k)^(-1) (C_k P_(k|k-1) A_k^T + S_k^T)
$ <eq:km_kf_Pk1k>

这大致可以看成是几项协方差以及变换后协方差的和的形式。我们也可以尽量把 $S_k$ 单独提出去，将剩下和 $S_k$ 无关的 $A_k P_(k|k) A_k^T + Q_k$ 部分，但其余项也很乱，这里就不写了。

#blockquote([
    *为什么要引入 $hat(x)_(k|k-1)$ 等中间估计量而不直接使用 $hat(x)_k$ 和 $hat(x)_(k-1)$*：

    首先，我们当然可以直接跳过 $hat(x)_(k|k)$，只要求 $hat(x)_(k+1|k)$ 最优，直接推导出从 $hat(x)_(k|k-1)$ 到 $hat(x)_(k+1|k)$ 的递推关系，即 one step ahead prediction。但 $hat(x)_(k|k)$ 才是我们需要的、结合了到 $k$ 为止所有已知信息的最优状态估计，我们终归要计算它。

    那么很自然地，我们想到能否直接利用从 $hat(x)_(k-1|k-1)$ 到 $hat(x)_(k|k)$ 的递推关系？如此直接使用 $hat(x)_(k-1)$ 和 $hat(x)_k$ 这些符号就可以了。其实也能做到，我们当然可以从 @eq:km_kf_xkk 的：
    
    // $
    // hat(x)_(k|k) = K_k y_k + (I_n - K_k C_k) hat(x)_(k|k-1)
    // $
    
    // 反过来得到：
    
    // $
    // hat(x)_(k|k-1) = (I_n - K_k C_k)^(-1) (hat(x)_(k|k) - K_k y_k)
    // $

    #Cre("TODO") 可以吗

    然而，通过 $hat(x)_(k|k)$ 倒解 $hat(x)_(k|k-1)$ 的过程让计算#underline[平白无故多出了一系列计算成本]，还不如正着推导时将其存下来，仅为简化符号多出这部分代价未免有些得不偿失。同时，引入中间估计量还提高了可解释性，将估计显式分为#underline[时间传播和测量融合]两个子阶段，更好地体现每一个输入数据和输出数据对估计及不确定性的具体影响。此外，#underline[有一系列场景需要用到这些中间估计量]，例如扩展 Kalman 滤波器中测量函数的线性化发生在#underline[测量前]状态估计的中心，即需要 $hat(x)_(k|k-1)$。
    
    综上，作为有实际物理意义的估计量，将它们留下比较好。
])

#blockquote([
    *几个多次出现的项的大致含义*：

    #underline[首先是] $(C_k P_(k|k-1) C_k^T + R_k)^(-1)$，其中 $P_(k|k-1)$ 是融合测量前对 $k$ 时刻估计的协方差，左右乘 $C_k$ 后的 $C_k P_(k|k-1) C_k^T$ 是通过测量模型之后的协方差，而 $R_k$ 是测量噪声的协方差。于是整个项的含义大致就是融合了测量模型和测量噪声之后的新协方差的逆，代表 “这次测量的可信度有多高”（取逆前则是 “有多不可靠”）。
    
    或者更严谨一点，测量前我们对状态的估计是 $hat(x)_(k|k-1)$，根据测量模型，我们对测量输出的预测就是 $hat(y)_k = C_k hat(x)_(k|k-1)$，而真实测量是 $y_k = C_k x_k + v_k$，二者相减得到残差：

    $
    y_k - hat(y_k) = C_k (x_k - hat(x)_(k|k-1)) + v_k
    $

    经过计算，它的协方差 $"Cov"(y_k - hat(y)_k) = C_k P_(k|k-1) C_k^T + R_k$ 就是那一项求逆里面的东西。

    类似地，关于 @eq:km_kf_Pk1k 中如的 $A_k P_(k|k-1) C_k^T + S_k$ 等项则是在描述，根据模型测量和我们关心的下一个估计天然有多相关，以及二者噪声之间的相关性。

    可以认识到，整个估计过程实际就是在一层层地剥离排除状态转移模型、测量模型以及二者噪声的影响，均反映在了估计量和协方差的形式上。
])

=== Similarity with Stochastic Least Squares

#Cre("TODO")

=== Case Without the Noise Cross-term

#Cre("TODO")

很多情况下，过程噪声和测量噪声之间并不相关，即可以假设 $S_k = 0$。在这样的前提下，Kalman 滤波器的形式将会更加简洁，我们直接将前文中涉及 $S_k$ 的项都去掉即可得到：

#Cre("TODO")

$

$

可以发现，时间演化和测量融合两个步骤之间的耦合项也消去了，于是我们可以将每一步估计明确地分为两个步骤，即测量更新（measurement update）：

#Cre("TODO")

$

$

以及时间更新（time update）：

#Cre("TODO")

$

$

=== Stationary Kalman Filter

#Cre("TODO") 只是时变时不变的问题吗？slides 里的定义似乎是 P 渐进稳定还是什么？

注意到前面解决的 Kalman 滤波问题涵盖了模型参数时变的情况。如果模型参数（包括噪声模型）不随时间变化，即符合：

$
x_(k+1) &= A x_k + B u_k + w_k, &quad &w_k &~ &cal(N)(0, Q) \
y_k &= C x_k + v_k, & &v_k &~ &cal(N)(0, R)
$

我们就可以直接设计好一个固定的 Kalman 增益，不用再随时间动态计算。

== Derivation of Kalman Filter

=== Unbiased and Minimum Variance

#Cre("TODO") 从 $hat(x)_(k|k-1)$ 到 $hat(x)_(k+1|k)$ 关系的证明参考：

// #image("/assets/image.png")
// #image("/assets/image-1.png")
// #image("/assets/image-2.png")
// #image("/assets/image-3.png")

=== Bayesian

#Cre("TODO")

== Extended Kalman Filter (EKF)

#Cre("TODO") 非线性就线性化。

== Innovation Form Representation <sec:km_innovation_form>

#Cre("TODO") 会在后面系统辨识的状态空间模型中被广泛采用。
