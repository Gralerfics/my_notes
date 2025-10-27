#import "../generic.typ": *

#import "@preview/cetz:0.4.2"
#import "@preview/cetz-plot:0.1.3"
#import "@preview/fletcher:0.5.8" as fletcher: diagram, node, edge

= Spectrum Estimation

功率谱的估计显然是个很有用的东西。例如，加性噪声、信号与噪声不相关前提下的非因果 Wiener 平滑滤波器有如下频率响应：

$
H(e^(j omega)) = (P_(d x) (e^(j omega))) / (P_x (e^(j omega))) = (P_d (e^(j omega))) / (P_d (e^(j omega)) + P_v (e^(j omega)))
$

如果目标信号与噪声的功率谱密度已知则可以直接求出其频率响应，若不知道则可以通过谱估计得到结果；再例如对窄频带信号的检测与追踪等。

总之本节考虑宽平稳随机过程的功率谱密度的估计。具体地，给定一个随机过程产生的随机信号，如何根据它有效地估计这个随机过程的谱？

最直接的想法是，由 Wiener–Khinchin Theorem，我们知道自相关函数就可以通过傅里叶变换得到功率谱：

$
P_x (e^(j omega)) = sum_(k=-infinity)^infinity r_x (k) e^(j k omega)
$

而对于一个遍历的随机过程产生的随机信号 $x[n]$，我们可以这样得到该过程的自相关函数：

$
r_x (k) = lim_(N->infinity) {1/(2N+1) sum_(n=-N)^N x[n+k] x^*[n]}
$

这是由于满足遍历性假设，我们用一条无限长的实现就可以无偏地估计出该过程的统计特征。

但显然这么做存在一些问题：第一，我们拥有的样本长度往往是有限的，例如地震波等本身就很短的信号，以及语音信号等只有很短时间内才近似满足平稳性假设的信号；第二，信号样本往往自己还包含噪声。

谱估计的方法可以分为两类：一类是无参数（Nonparametric）的方法，如从估计序列的自相关入手，变换得到功率谱；另一类是在对随机过程的模型有先验了解的情况下可以使用的有参（Parametric）估计的方法，从估计模型参数入手，再由模型计算功率谱。

== Nonparametric Spectrum Estimation

=== Periodogram <sec:se_nonparam_periodogram>

本章开头提到从样本中估计自相关函数，然后傅里叶变换得到功率谱的方法。即使现在样本长度有限（例如 $N$ 个），我们也就用这点样本来直接估计自相关函数：

$
hat(r)_x (k) = 1/N sum_(n=0)^(N-1-k) x[n+k] x^*[n], quad k=0, 1, dots, N-1
$

此处实际上是假定了 $x[n]$ 在 $[0, N-1]$ 以外的部分值都为 $0$，等效于将求和范围设置为只从 $0$ 到 $N-1-k$，即求和共有 $N-k$ 项。但由于前面除以的是 $N$，最终得到的是一个有偏（Biased）估计：

$
E{hat(r)_x (k)} &= 1/N sum_(n=0)^(N-1-k) E{x[n+k] x^*[n]} \
&= 1/N sum_(n=0)^(N-1-k) r_x (k) \
&= (N-k)/N r_x (k), quad k=0, 1, dots, N-1
$

当 $N$ 趋于无穷时期望就和实际值相等了，所以它是渐近无偏（Asympotically Unbiased）的。$k$ 为负的部分我们可以由 WSS 过程自相关函数的共轭对称性质翻转得到，于是可以将其重写为：

$
E{hat(r)_x (k)} = w_B (k) r_x (k)
$ <equ:se_E_hatrxk_bw_form>

其中：

$
w_B (k) = cases(
    (N-abs(k))/N\, quad &abs(k)<=N,
    0\, &abs(k)>N
)
$

这是一个 Bartlett (triangular) window。

#blockquote[
    #text(fill: red, "（TODO）")为何不直接换成除以 $N-k$ 得到无偏估计？

    也许和 Parseval 一致性有关：

    $
    1 / (2 pi) integral_(-pi)^pi hat(P) (omega) dif omega = hat(r)_x (0) = 1 / N sum_(n=0)^(N-1) abs(x[n])^2
    $

    此外，这里似乎有一些循环定义的问题，暂时不纠结。

    #text(fill: red, "（TODO）")还有降低方差的考虑。以及后面有参估计时会发现如果用除以 $N-k$ 的无偏估计会导致协方差矩阵不确保正定。所以这里的方差表现差不只是量变还会有质变是吗。
]

接下来再对它使用傅里叶变换估计功率谱：

$
hat(P) (e^(j omega)) = sum_(k=-N+1)^(N-1) hat(r)_x (k) e^(-j k omega)
$

以上这个两步走的过程主要是顺着定义来，了解原理后我们可以化简一下直接通过随机信号 $x[n]$ 求得上述频谱估计。首先我们对 $x[n]$ 在 $[0, N-1]$ 以外的部分值都为 $0$ 的假设可以用乘以一个 Dirichlet 核（@equ:fun_dsp_sa_dirichlet_kernel）的过程替代，即 $x_N [n] = x[n] w_R [n]$，其傅里叶变换为：

#emphasis_equbox([
$
X_N (e^(j omega)) = sum_(n=-infinity)^infinity x_N [n] e^(-j omega n) = sum_(n=0)^(N-1) x[n] e^(-j omega n)
$
])

由：

$
hat(r)_x (k) = 1/N sum_(n=-infinity)^infinity x_N [n+k] x_N^* [n] = 1/N x_N [k] * x_N^* [-k]
$

再由 DTFT 的性质，可得功率谱估计为：

#emphasis_equbox([
$
hat(P)_"per" (e^(j omega)) = 1/N X_N (e^(j omega)) X_N^* (e^(j omega)) = 1/N abs(X_N (e^(j omega)))^2
$ <equ:se_periodogram_estimator_from_X_N>
])

将满足该定义式的功率谱估计称为*周期图*（Periodogram）。

==== An Equivalent Perspective from Filter Banks

#text(fill: red, "（TODO）")书 396 页。先说结论。

我们定义一组长度为 $N$ 的 FIR 滤波器如下：

$
h_i [n] = 1/N e^(j n omega_i) w_R [n] = cases(
    1/N e^(j n omega_i)\, quad &0<=n<N,
    0\, &"otherwise"
)
$

==== Performance of the Periodogram

我们接下来评估用周期图作为功率谱估计的表现。首先我们一定希望从样本中计算出的周期图可以收敛到实际随机过程的功率谱。由于其随机性质，我们需要从统计意义上考虑收敛性，如均方收敛：

$
lim_(N->infinity) E{[hat(P)_"per" (e^(j omega)) - P_x (e^(j omega))]^2} = 0
$

要满足这一点，我们需要其均值渐近无偏、方差在样本量足够大时趋于零：

$
lim_(N->infinity) E{hat(P)_"per" (e^(j omega))} = P_x (e^(j omega)) \
lim_(N->infinity) "Var"{hat(P)_"per" (e^(j omega))} = 0
$

换句话说，我们希望周期图是对功率谱密度的一致估计（Consistent Estimate）。*先说结论*，以上第一个条件成立，但第二个条件不成立。

具体地，首先考虑*第一个条件*，即渐近无偏。从样本自相关开始，其期望如 @equ:se_E_hatrxk_bw_form 所示，再由周期图的推导可以得到：

$
E{hat(P)_"per" (e^(j omega))} &= E{sum_(k=-N+1)^(N-1) hat(r)_x (k) e^(-j k omega)}
= sum_(k=-N+1)^(N-1) E{hat(r)_x (k)} e^(-j k omega) \
&= sum_(k=-N+1)^(N-1) w_B (k) r_x (k) e^(-j k omega)
$

这就是自相关函数和窗函数乘积的傅里叶变换，于是由频域卷积性质得到：

#emphasis_equbox([
$
E{hat(P)_"per" (e^(j omega))} = 1/(2 pi) P_x (e^(j omega)) * W_B (e^(j omega))
$ <equ:se_periodogram_expectation>
])

其中 Bartlett 窗的频域表达式是：

$
W_B (e^(j omega)) = 1/N [sin(N omega\/2)/sin(omega\/2)]^2
$

#blockquote[
    注意这里的窗函数是加在自相关序列上的，所以称为滞后窗（Lag window），区别于之后直接加在数据上的数据窗（Data window）。无实际意义，仅作概念上的区分。
]

注意到 $N -> infinity$ 时，$W_B (e^(j omega))$ 收敛到一个脉冲函数（周期积分为 $2 pi$，面积集中到原点，即 $2 pi delta(omega)$，证明略），于是周期图满足渐近无偏的条件。

从频谱图像的角度，本来是用一堆理想的脉冲函数拼出整个谱，现在变成用 Bartlett 窗函数来拼出整个谱。样本越多，窗函数越接近理想脉冲信号，频谱估计越准确，如 @fig:se_bartlett_window_freq 所示。

#figure(
    caption: [The Fourier transforms of Bartlett windows]
)[
    #cetz.canvas({
        import cetz.draw: *
        import cetz-plot: plot

        let samples = 2048
        
        let func(N) = (w) => {
            let sden = calc.sin(w / 2)
            if (calc.abs(sden) < 1e-5) {
                N
            } else {
                (1.0 / N) * (calc.pow(calc.sin(N * w / 2), 2)) / (calc.pow(sden, 2))
            }
        }

        set-style(
            axes: (stroke: .5pt, tick: (stroke: .5pt)),
            legend: (stroke: none, orientation: ttb, item: (spacing: .1), scale: 80%)
        )

        plot.plot(
            size: (12, 6),
            x-tick-step: calc.pi / 6,
            x-format: plot.formats.multiple-of,
            x-label: $omega$,
            y-tick-step: 20,
            y-min: -0.2,
            y-max: 100.0,
            y-label: $W_B (e^(j omega))$,
            legend: "inner-north-east",
            {
                let domain = (-calc.pi / 3, calc.pi / 3)
                plot.add(
                    func(32),
                    domain: domain,
                    samples: samples,
                    label: $N = 32$,
                    style: (stroke: black)
                )
                plot.add(
                    func(64),
                    domain: domain,
                    samples: samples,
                    label: $N = 64$,
                    style: (stroke: blue)
                )
                plot.add(
                    func(256),
                    domain: domain,
                    samples: samples,
                    label: $N = 256$,
                    style: (stroke: purple)
                )
            }
        )
    })
] <fig:se_bartlett_window_freq>

该函数的主瓣带宽大约为 $(2 pi)/N$，太过靠近的两个脉冲信号在同窗函数卷积后可能导致两个峰叠加合并，无法清晰分辨。故定义分辨率为此处窗函数的 $6"dB"$ 带宽：

#emphasis_equbox([
$
"Res"[hat(P)_"per" (e^(j omega))] = (Delta omega)_(6"dB") = 0.89 (2 pi)/N
$
])

$-6"dB"$ 约为 $0.5$，即大约在这个位置两峰交叠处的值为峰值的一半，考虑 @equ:se_periodogram_expectation 中的卷积过程，这会导致两峰叠加后只剩一个峰而使结果无法分辨。

#blockquote[
    考虑在实际计算中的一个细节问题，即这里给出的 $Delta omega$ 的数值，是在傅里叶变换归一化到 $[-pi, pi]$ 区间后的测度意义下的。

    而通常我们的指标是在实际频谱单位下的，这时我们需要除以采样率。例如，一个采样率为 $10"kHz"$ 的信号，我们指出要求分辨率至少达到 $10"Hz"$，那么 $Delta omega$ 应该是 $2 pi times (10"Hz")/(10"kHz")$。
]

*第二个条件*，考虑方差是否趋于零。由于周期图同样本是二阶关系，现在又计算方差，就是对随机过程四阶矩的计算，这太过复杂；但我们可以考虑随机过程是方差为 $sigma_x^2$ 的高斯白噪声（Gaussian white noise）的特殊情况。经过一系列计算（详见参考书第 8.2.2 节，404 页）可得该情况下的方差为：

$
"Var"{hat(P)_"per" (e^(j omega))} = sigma_x^4
$

与 $N$ 无关，不会随其增长收敛到零。实际上，如果考虑普遍情况，我们有如下近似（对于高斯白噪声的情况就是 $(sigma_x^2)^2$）：

#emphasis_equbox([
$
"Var"{hat(P)_"per" (e^(j omega))} approx P_x^2 (e^(j omega))
$ <equ:se_periodogram_variance>
])

所以结论就是，第二个条件并不满足，即*周期图不是对功率谱密度的一致估计*。

=== Modified Periodogram

所以我们当然想着要改进一下，我们先不管前面的各种推导，直接考虑在定义上修一修，然后再验证效果。回过来观察周期图的定义式：

$
hat(P)_"per" (e^(j omega)) = 1/N abs(X_N (e^(j omega)))^2 = 1/N abs(sum_(n=-infinity)^(infinity) x[n] #text(fill: blue, $w_R [n]$) e^(-j omega n))^2
$

式中体现的是我们对原信号使用 Dirichlet 核（也就是方形窗，Rectangular window）后进行谱估计的过程。那么一个很直观的想法是，如果这里用别的窗函数会不会有什么效果？

#text(fill: red, "（TODO）")书 408、409 页在推期望和方差。注意书上好像有关于 $w_B (k)$ 的定义不一致性，我们这里采用归一化的 $w_B (k) = 1/N w_R (k) * w_R (-k) = sum_(n=-infinity)^infinity w_R (k) w_R (n-k)$，和前面保持一致，所以和书上略有系数的差异，要改。

我们定义修正周期图（Modified periodogram）为：

#emphasis_equbox([
$
hat(P)_M (e^(j omega)) = 1/(N U) abs(sum_(n=-infinity)^infinity x[n] w[n] e^(-j n omega))^2
$
])

其中 $N$ 为窗函数的长度，常数 $U$ 为窗函数的平均能量（后续会说明这是为了使修正周期图渐近无偏）：

#emphasis_equbox([
$
U = 1/N sum_(n=0)^(N-1) abs(w[n])^2
$
])

==== Performance of Modified Periodogram

类似地，我们评估修正周期图的表现。首先是偏差（Bias），由类似的推导有：

#emphasis_equbox([
$
E{hat(P)_M (e^(j omega))} = #text(fill: blue, $1/(2 pi N U)$) P_x (e^(j omega)) * #text(fill: blue, $abs(W(e^(j omega)))^2$)
$ <equ:se_modified_periodogram_expectation>
])

其中 $W(e^(j omega))$ 是 $w[n]$ 的傅里叶变换，故由前 $U$ 的设置有：

$
U = 1/N sum_(n=0)^(N-1) abs(w[n])^2 = 1/(2 pi N) integral_(-pi)^pi abs(W(e^(j omega)))^2 dif omega
$

即：

$
integral_(-pi)^pi #text(fill: blue, $1/(2 pi N U) abs(W(e^(j omega)))^2$) dif omega = 1
$

这使得 $E{hat(P)_M (e^(j omega))}$ 在 $N->infinity$ 时趋于功率谱密度，即渐近无偏，也是这样设置 $U$ 的目的。

接下来是方差，加数据窗并不有助于降低方差，故同 @equ:se_periodogram_variance 一样：

#emphasis_equbox([
$
"Var"{hat(P)_M (e^(j omega))} approx P_x^2 (e^(j omega))
$
])

即修正周期图同样不是对功率谱密度的一致估计。

==== Trade-off between Resolution and Confusion

不影响估计的偏差与方差，那么加数据窗的操作到底对什么有影响？

不同数据窗的傅里叶变换形态有差异，主要体现在主瓣（Main lobe）和旁瓣（Sidelobe）上。参考 @sec:fun_dsp_sa，前者将影响谱估计的分辨率（Resolution），后者则会引入旁瓣的干扰和混淆。

我们定义分辨率为数据窗主瓣的 $3"dB"$ 带宽，这个值越大说明越不清晰：

$
"Res"[hat(P)_"per" (e^(j omega))] = (Delta omega)_"3dB"
$

#blockquote[
    注意到，前文分析 Periodogram 时定义的分辨率是 Bartlett 窗的 $6"dB"$ 带宽而不是 $3"dB"$，但实际上这是一致的。

    这是由于前文分析的是加在自相关序列上的窗（即滞后窗），注意 @equ:se_periodogram_expectation 中同功率谱密度卷积的是 $W_B (e^(j omega))$；而此处分析的是数据窗，注意 @equ:se_modified_periodogram_expectation 中同功率谱密度卷积的是 $1/(N U) abs(W (e^(j omega)))^2$，二者之间存在平方关系。

    所以前者的 $-6"dB"$ 点同后者的 $-3"dB"$ 点是一致的，都是相对信号来说的半功率点。
]

常见窗函数旁瓣抑制和分辨率的大致值总结如 @tab:se_properties_of_data_windows 所示。

#figure(
    caption: "Properties of a few commonly used windows with length N"
)[
    #table(
        columns: 3,
        stroke: (x, y) => if y == 0 {
            (bottom: 0.7pt + black)
        },
        align: (x, y) => (
            if x > 0 { center }
            else { left }
        ),
        table.header(
            [], [Sidelobe (dB)], [Resolution]
        ),
        "Rectangular", $-13$, $0.89 (2 pi \/ N)$,
        "Bartlett", $-27$, $1.28 (2 pi \/ N)$,
        "Hanning", $-32$, $1.44 (2 pi \/ N)$,
        "Hamming", $-43$, $1.30 (2 pi \/ N)$,
        "Blackman", $-58$, $1.68 (2 pi \/ N)$,
    )
] <tab:se_properties_of_data_windows>

可以观察到往往旁瓣抑制效果越好，分辨率就越差（即$3"dB"$带宽越大），这是一个 Trade-off。

=== Periodogram Averaging

至此，由于方差不收敛，以上方法都无法得到关于功率谱密度的一致估计。而接下来我们通过几种对周期图做平均的方法来得到我们想要的一致估计。

考虑之前对于一个随机变量 $x$，我们通过收集其大量不相关的测量样本并计算样本均值的方式，得到了对该随机变量均值的一致估计 $E{x}$。类比之，理论上我们需要用随机过程 $x(n)$ 的多个不相关的实现（Uncorrelated realizations），分别求周期图后再平均，估计周期图的期望。

具体地，设我们有 $K$ 个不相关的实现 $x_i [n]$，每个长度为 $L$，有总样本点数为 $N = L K$。计算每个实现的周期图：

$
hat(P)_"per"^((i)) (e^(j omega)) = 1/L abs(sum_(n=0)^(L-1) x_i [n] e^(-j n omega))^2, quad i = 0, 1, dots, K-1
$

然后求平均得最终的谱估计：

$
hat(P)_x (e^(j omega)) = 1/K sum_(i=1)^K hat(P)_"per"^((i)) (e^(j omega)) = 1/N sum_(i=0)^(K-1) abs(sum_(n=0)^(L-1) x_i [n] e^(-j n omega))^2
$

下面照例评估其偏差和方差。首先，因为只是再做了一次平均，期望同前面是一样的：

#emphasis_equbox([
$
E{hat(P)_x (e^(j omega))} = E{hat(P)_"per"^((i)) (e^(j omega))} = 1/(2 pi) P_x (e^(j omega)) * W_B (e^(j omega))
$
])

故在 $L -> infinity$ 时渐近无偏。然后考虑方差，由于不同实现之间不相关，有：

#emphasis_equbox([
$
"Var"{hat(P)_x (e^(j omega))} &= 1/K^2"Var"{sum_(i=0)^(K-1) hat(P)_"per"^((i)) (e^(j omega))} \
&= 1/K"Var"{hat(P)_"per"^((i)) (e^(j omega))} approx 1/K P_x^2 (e^(j omega))
$
])

综上，使用这种基于平均的方法可以在 $L$ 和 $K$ 都趋于无穷时给出对功率谱密度的*一致估计*。

==== Bartlett's Method

一般实际情况中我们没有那么多独立的实现，但若我们有一条足够长的实现，并且其背后的随机过程满足遍历性假设，那么我们就可以把它切成小段当作不相关的实现来用，求得谱估计，称为 Bartlett 法。

设信号长度为 $N$，切成*不重叠*（尽量保证不相关）的 $K$ 段，每段长度 $L$。若我们再令：

$
x_i [n] = x[n + i L], quad n = 0, 1, dots, L-1; quad i = 0, 1, dots, K-1
$

符号就同前面的分析一致了，我们直接得到公式：

#emphasis_equbox([
$
hat(P)_B (e^(j omega)) &= 1/K sum_(i=0)^(K-1) (1/L abs(sum_(n=0)^(L-1) x[n + i L] e^(-j n omega))^2) \
&= 1/N sum_(i=0)^(K-1) abs(sum_(n=0)^(L-1) x[n + i L] e^(-j n omega))^2
$
])

*偏差和方差则同前文平均法的分析一致*。不过由于即便不重叠，切分出的序列间也一定存在相关性（遍历性只是保证在较长时间下相关性会逐渐消弭），所以方差可能比前面分析的还小一些。不过我们还是近似地认为片段之间不相关，取这个方差的近似值。

接下来我们再分析一下该法对分辨率的影响，由于原本长度为 $N$ 的序列被我们切成长度为 $L$ 的小段了，所以计算周期图时实际序列长度只有 $L$，故分辨率为：

$
"Res"[hat(P)_B (e^(j omega))] = 0.89 (2 pi)/L = 0.89 K (2 pi)/N
$

可以看到相比用整个长度为 $N$ 的序列计算周期图，分辨率变差了 $K$ 倍，这是代价。

==== Welch's Method

Barlett 法用周期图做平均，接下来我们沿用这个思想，用修正周期图做平均，称为 Welch 法。

修正周期图的想法是对数据加窗函数，而这一操作实际可以减弱相邻两个片段边缘处信号的相关性。所以我们可以考虑放宽要求，允许切分数据时有*重叠*（Overlapping）。每个片段依旧长度为 $L$，但每段的起点只间隔 $D$ 个样本（$D<L$ 时就有重叠了）：

$
x_i [n] = x[n + i D], quad n = 0, 1, dots, L-1; quad i = 0, 1, dots, K-1
$

于是总样本点数就为 $N = L + D(K - 1)$。通常，我们取 $50%$ 重叠，即 $D = L \/ 2$。数据窗加到每个片段上，窗长度同片段长度 $L$ 一致。于是我们得到 Welch 法的谱估计为：

#emphasis_equbox([
$
hat(P)_W (e^(j omega)) &= 1/K sum_(i=0)^(K-1) (1/(L U) abs(sum_(n=0)^(L-1) w[n] x[n + i D] e^(-j n omega))^2) \
&= 1/(K L U) sum_(i=0)^(K-1) abs(sum_(n=0)^(L-1) w[n] x[n + i D] e^(-j n omega))^2
$
])

照例分析指标，首先是偏差，与 Modified periodogram 一致，将 $N$ 换成 $L$：

#emphasis_equbox([
$
E{hat(P)_W (e^(j omega))} = E{hat(P)_M^((i))} = 1/(2 pi L U) P_x (e^(j omega)) * abs(W(e^(j omega)))^2
$
])

方差与重叠程度有关，难以计算，我们只考虑 $50%$ overlap 的情况：

#emphasis_equbox([
$
"Var"{hat(P)_W (e^(j omega))} approx 9/(8 K) P_x^2 (e^(j omega)) approx 9/16 L/N P_x^2 (e^(j omega))
$
])

系数 $9\/8$ 看似表明方差表现变差了一点，但由于片段数量 $K approx (2 N)/L$ 差不多翻了一倍，方差表现实际上是提升的。

=== Periodogram-based Methods Summary

参考书中还提到了一个 Blackman-Tukey 法，此处略去不表。前面除此之外的每种方法我们都分析了其均值、方差、分辨率等表现，在此做一个总结，如 @tab:se_periodogram_formula_summary 所示。

#figure(
    caption: "Properties of a few commonly used windows with length N"
)[
    #resize_box([
        #table(
            columns: 5,
            stroke: (x, y) => if y == 0 {
                (bottom: 0.7pt + black)
            },
            align: horizon,
            table.header(
                [], [Definition $hat(P)_x (e^(j omega))$], [Expectation], [Variance (approx.)], [Resolution $Delta omega$]
            ),
            "Periodogram",
            $display(1/N abs(sum_(n=-infinity)^(infinity) x[n] w_R [n] e^(-j omega n))^2)$,
            $display(1/(2 pi) P_x (e^(j omega)) * W_B (e^(j omega)))$,
            $display(P_x^2 (e^(j omega)))$,
            $display(0.89 (2 pi)/N)$,

            "Modified periodogram",
            $display(1/(N U) abs(sum_(n=-infinity)^infinity x[n] w[n] e^(-j omega n))^2)$,
            $display(1/(2 pi N U) P_x (e^(j omega)) * abs(W(e^(j omega)))^2)$,
            $display(P_x^2 (e^(j omega)))$,
            [See @tab:se_properties_of_data_windows],

            [Bartlett's \ $N = K L$],
            $display(1/N sum_(i=0)^(K-1) abs(sum_(n=0)^(L-1) x[n + i L] e^(-j omega n))^2)$,
            $display(1/(2 pi) P_x (e^(j omega)) * abs(W_B (e^(j omega)))^2)$,
            $display(1/K P_x^2 (e^(j omega)))$,
            $display(0.89 K (2 pi)/N)$,

            [Welch's ($50%$ overlap) \ $N = L + D(K - 1)$],
            $display(1/(K L U) sum_(i=0)^(K-1) abs(sum_(n=0)^(L-1) w[n] x[n + i D] e^(-j omega n))^2)$,
            $display(1/(2 pi L U) P_x (e^(j omega)) * abs(W(e^(j omega)))^2)$,
            $display(9/16 L/N P_x^2 (e^(j omega)))$,
            "Window dependent"
        )
    ])
] <tab:se_periodogram_formula_summary>

注意，其中修正周期图的 $U = 1/N sum_(n=0)^(N-1) abs(w[n])^2$，而 Welch 法中由于每个片段长度为 $L$，其 $U = 1/L sum_(n=0)^(L-1) abs(w[n])^2$。

参考书中定义了两个指标来衡量以上方法的表现，其一为变异性（Variablilty）：

$
cal(V) = ("Var"{hat(P)_x (e^(j omega))})/(E^2{hat(P)_x (e^(j omega))})
$

说白了就是归一化的方差。其二是品质因数（Figure of merit）：

$
cal(M) = cal(V) Delta omega
$

是变异性和分辨率的乘积，这个值越小越好。顺带一提，这里品质因数这种定义为两个量相乘的指标，一般都是把 Trade-off 的量乘起来，所以我们会发现前面这些无参估计方法的品质因数都差不太多。

=== Minimum Variance (MV) Spectrum Estimation



#text(fill: red, "（TODO）")书 8.3 节。

== Parametric Spectrum Estimation

=== For Autoregressive Moving Average Models

=== For Autoregressive Models

=== Frequency Estimation 

==== Problem Setup

==== Mulitple Signal Clasification (MUSIC)
