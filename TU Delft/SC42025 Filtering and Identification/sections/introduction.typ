#import "../generic.typ": *

#import "@preview/cetz:0.4.2"
#import "@preview/fletcher:0.5.8" as fletcher: diagram, node, edge

= Introduction

== Key Problems

课程涉及解决三类主要问题。*第一*是从数据中估计静态模型的参数（estimate static models from data）：

$
min_#Cre($theta$) norm(y - F #Cre($theta$))^2_W
$

将模型输出 $hat(y) = F theta$，与实际输出 $y$ 做差后求最小二乘，找到使误差最小的一组模型参数 $theta$。

#blockquote([
    *对一般线性模型的实例注释*：

    若定要从模型有输入和输出的角度去分辨，这里有输出 $y$ 和参数 $theta$，而无显式的输入 $x$ 会显得有一些奇怪。但实际上，#underline[$F$ 就是包含了所有数据信息的矩阵]，是一个精神上的输入 $x$，所以可以叫它 “data matrix”；相应地，$y$ 是与输入对应的测量值，就叫 "measurement vector"。
    
    // 或可反过来将输入数据 $x$ 视为一个滤波器 $F$，应用到模型参数 $theta$ 上。

    对于线性运算，我们*怎么写都是可以的*，例如这里的 $F theta$，或者将参数写成线性算子（矩阵）$Theta$ 应用到 $x$ 上写成 $Theta x$，区别只是把 $x$ 还是把 $theta$ 拆开拼成矩阵。举一个线性卷积（或者就 FIR 滤波器）的例子，同样是表达 $y_n = sum_(i = 1)^p theta_i x_(n-i+1)$，我们可以把它写成 $y = F theta$：

    $
    mat(
        delim: "[",
        #Cpu($y_n$) ;
        y_(n-1);
        y_(n-2);
        dots.v;
        y_p;
    )
    =
    mat(
        delim: "[",
        #Cre($x_n$), x_(n-1), x_(n-2), dots, #Cbl($x_(n-p+1)$) ;
        x_(n-1), x_(n-2), x_(n-3), dots, x_(n-p);
        x_(n-2), x_(n-3), x_(n-4), dots, x_(n-p-1);
        dots.v, dots.v, dots.v, dots.down, dots.v;
        x_p, x_(p-1), x_(p-2), dots, x_1;
    )
    mat(
        delim: "[",
        #Cre($theta_1$) ;
        theta_2;
        theta_3;
        dots.v;
        #Cbl($theta_p$) ;
    )
    := F theta
    $
    
    也可以写成 $y = Theta x$：

    $
    mat(
        delim: "[",
        #Cpu($y_n$) ;
        y_(n-1);
        y_(n-2);
        dots.v;
        y_p;
    )
    =
    mat(
        delim: "[",
        #Cre($theta_1$), theta_2, theta_3, dots, #Cbl($theta_p$), 0, 0, dots, 0;
        0, theta_1, theta_2, dots, theta_(p-1), theta_p, 0, dots, 0;
        0, 0, theta_1, dots, theta_(p-2), theta_(p-1), theta_p, dots, 0;
        dots.v, dots.v, dots.v, dots.down, dots.v, dots.v, dots.v, dots.down, dots.v;
        0, 0, 0, dots, 0, 0, 0, dots, theta_1;
    )
    mat(
        delim: "[",
        #Cre($x_n$) ;
        x_(n-1);
        x_(n-2);
        dots.v;
        #Cbl($x_(n-p+1)$) ;
        x_(n-p);
        x_(n-p-1);
        dots.v;
        x_1;
    )
    := Theta x
    $

    #underline[但由于我们在该问题中想要求解的是参数 $theta$]，基于上述认识，将其作为向量形式写入问题中会更方便一些，求解一个矩阵比求解一个向量不直观得多。
    
    // 我们预先设定了这是一个静态模型，即输入输出数据集给定，都是 deterministic 的，所以可以任意地对 $x$ 进行操作，将其组合成一个 $F$。$x$ 是包含所有输入信息的常量，$F$ 也是包含所有输入信息的常量，只是形式不一样。
])

*第二*是滤波（Filtering），已知模型和参数前提下，估计时变状态空间模型的状态（estimate the state of a time-varying state space model）：

$
#Cre($x_(k+1)$) &= A_k #Cre($x_k$) + B_k u_k + w_k quad quad &w_k ~ cal(N)(0, Q) \
y_k &= C_k #Cre($x_k$) + D_k u_k + v_k &v_k ~ cal(N)(0, R)
$

即模型参数 $A_k, B_k, C_k, D_k, Q, R$ 及输入 $u_(1:k)$ 和输出 $y_(1:k)$ *已知*，*要求*找到最符合当前观察结果（输出）的状态序列 $x_k$。模型参数 $A_k, B_k, C_k, D_k$ 的下标 $k$ 表明此处#underline[考虑时变的可能]。

相较输入和输出，#underline[状态往往隐藏在系统内部]，例如一个含噪声的测量系统，实际值作为隐藏的内部状态，测量值是经过噪声干扰后的输出，我们需要通过输出和对噪声的认识（模型和模型参数）估计出实际值，故称滤波。

*第三*是系统辨识（Identification），在给定输入和输出的情况下，同时估计状态空间模型参数和状态：

$
#Cre($x_(k+1)$) &= #Cre($A$) #Cre($x_k$) + #Cre($B$) u_k + w_k quad quad &w_k ~ cal(N)(0, #Cre($Q$)) \
y_k &= #Cre($C$) #Cre($x_k$) + #Cre($D$) u_k + v_k &v_k ~ cal(N)(0, #Cre($R$))
$

即仅输入 $u_(1:k)$ 和输出 $y_(1:k)$ *已知*，*要求*估计状态序列 $x_k$ 和模型参数 $A, B, C, D, Q, R$。估计是问题比较复杂，这里只考虑#underline[模型参数时不变]的情况。

对于一个未知系统我们只能进行输入并得到输出的测试，假定它满足状态空间模型的形式，我们要通过这些测试结果，将猜的模型的参数填出来，从现象中分辨出模型的行为模式，故称辨识。

视情况这个过程也可以叫校准（calibration），比如已经知道模型长什么样的传感器，在使用前需要通过类似的过程从测试数据中校准一下参数，本质也是系统辨识的过程。

#blockquote([
    *有关 “模型” 概念*：

    #Cre("TODO") well，Lec3 中还写了非线性的情况

    前面的问题分类中，对于笼统的 "模型" 概念，都将其具象化为了*状态空间模型*（state space models），而且还是*离散*的情况。
    
    但要注意的是，状态空间模型*只是一类模型描述*，一种通用的线性系统描述形式，通常用来描述动态系统与输入输出之间的联系（动态系统，dynamic systems，即当前状态与过去状态有关的系统）。这里的离散状态空间模型涵盖了时变/时不变的离散线性系统，*但它一般无法表示*：
    
    + #underline[非线性系统]：一般在某点处对系统进行线性化处理之后用状态空间模型表示，但无法直接描述。
    
    + #underline[非因果系统]：定义上就限制了它无法描述出 $x_(k+1)$ 同 $x_(k+2)$ 等未来状态相关的模型。#underline[顺便]，虽然状态方程中只有 $x_(k+1)$ 同 $x_k$ 的关系，看起来无法描述当前状态和之前多个状态的关系（即让系统的记忆长一些）；但我们可以把过去多个状态 ${x_k, x_(k-1), dots, x_(k-p+1)}$ 合在一起作为一个状态向量 $bold(x)_k$，#underline[类似于系统在带着一个记忆滑动窗口沿着时间演化]。

    + #underline[无限维系统]：状态向量长度是有限的，无法描述例如空间上无限维、时间上无限记忆等情况。#underline[例如]，状态是一个连续温度场等，常出现在PDE系统等；#underline[再例如]，前一条中提到的将过去多个状态综合成一个状态向量可以做到记忆的延长，但由于向量 $x_k$ 长度有限，这种技巧无法处理记忆无限长的情况。
    
    + ……
])

== About Filtering

课件举了一下传感器融合（sensor fusion）的例子，实际上也同时关于滤波和辨识问题；列了一下后面会讲的方法，涉及线性和非线性系统。

== About System Identification

#Cre("TODO") 见后具体章节吧。

=== Parametric System Identification

#Cre("TODO")

// ${(u(k), y(k))}_(k=1)^N$

// $
// y_k approx hat(y)(k|k-1, theta)
// $

=== Two Schools of Linear System Identification

两个方向的线性系统辨识方法，一是*预测误差法*（prediction error methods），二是*子空间法*（subspace methods）。如前面注释中提到的，可以有各种不同的模型假设，这里*前者*常采用传递函数模型：

$
hat(y)_k (theta) = hat(G)(q|theta) u_k + hat(H)(q|theta) y_k
$

然后通过求解非线性优化问题来辨识参数：

$
min_theta sum_(k=1)^N norm(y_k - hat(y)_k (theta))^2
$

*后者*则采用状态空间模型，避开非线性优化问题，通过线性代数的手段求解系统辨识问题。
