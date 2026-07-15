#import "../generic.typ": *

#import "@preview/cetz:0.4.2"
#import "@preview/fletcher:0.5.8" as fletcher: diagram, node, edge

= Papers

== Learning Fine-Grained Bimanual Manipulation with Low-Cost Hardware

https://blog.csdn.net/nenchoumi3119/article/details/147504821

https://zhuanlan.zhihu.com/p/1981673900472557572

=== ALOHA 低成本开源双臂遥操作系统

略，主（leader）和从（follower）机械臂。

=== Action Chunking with Transformers (ACT)

*数据采集*。在 ALOHA 上采集人类演示数据，记录 leader 机器人的关节坐标（joint positions），而非 follower 的（“人希望机器人去哪里” 而非 “此刻真实在哪里”）。

*观测*。当前 follower 的关节位置和 4 摄像头的图像。

*预测*。训练 ACT，依据当前观测预测一段动作序列。这里动作（action）是指下一个时间步双臂的目标关节位置。测试时采用验证集误差最低的策略生成动作序列。主要的挑战是累计误差，即前一个动作的误差可能导致偏离训练分布。

*动作分块（action chunking）*。以 $k$ 个时间步为单位从观测中生成动作，即策略分布为 $pi_theta (a_(t:t+k) mid(|) s_t)$ 而非 $pi_theta (a_t mid(|) s_t)$。一个动作块可以包含一段连续的小动作任务，并且可以处理演示中的非马尔可夫行为（比如演示中途停顿）等与时间相关的混杂因素（temporally correlated confounders），因为行为不止依赖于状态还依赖于时间步。与历史条件策略（history-conditioned policy）的区别在于，历史条件策略是根据历史信息进行预测，可能导致因果混淆（即学习到错误的模式，认为是历史某些因素决定了未来动作）；而动作分块主要还是依据当前观测进行预测，只是让输出可以囊括一段动作。

*时序集成（temporal ensemble）*。刻板地每 $k$ 步观测并预测一次显得不连续，可能导致运动不稳定。考虑每一步都查询策略生成动作块，对于指定时间步将其在不同动作块中的预测进行加权平均。权重为 $w_i = exp(-m i)$，$m$ 越小新观测合并得越快。

*建模人类数据*。人类行为在精度要求较低时更加随机，应该学习 “给定观测下合理动作序列的分布” 而非一个唯一的确定策略，考虑将其建模为一个生成模型，本文采用条件变分自编码器（CVAE）。普通自编码器（autoencoder）将输入压缩（编码）到隐变量 $z$ 再重建（解码）原输入；变分自编码器（VAE）则将输入编码为一个概率分布（具体例如高斯分布的参数），从中采样后重建输入；CVAE 则加入条件将输入编码为一个条件分布，以应对本文中需要依据观测（条件）生成动作的需求。

具体结构上，CVAE 包含一个编码器和一个解码器，前者只在训练后者时使用，测试时丢弃。CVAE 编码器在#underline[训练时]基于当前观测 $macron(o)_t$（不含图像的观测，主要是关节状态）和动作序列 $a_(t:t+k)$ 推断样本对应的潜变量 $z$，称为风格变量（style variable），最终记为 $q_(phi.alt) (z mid(|) a_(t:t+k), macron(o)_t)$。CVAE 解码器则依据当前观测 $o_t$ 和潜变量 $z$ 生成动作序列，即 $pi_theta (hat(a)_(t:t+k) mid(|) o_t, z)$。

*实现 ACT*。CVAE 编码器采用 transformer encoder 实现，输入为当前观测以及目标动作（加上位置编码），"[CLS]" 是类 BERT transformer encoder 中常用的特殊 token，在经过多次注意力机制后会得到包含其他所有输入 token 信息的特征，被用于预测 $z$ 的分布。CVAE 解码器包含数个 ResNet 图像编码器、一个 transformer encoder 和一个 transformer decoder，提取图像特征后带上位置编码和关节状态（具体一点，这里是两个 7 自由度机械臂的关节位置共计 14 个状态）以及 $z$ 送入 transformer encoder 编码；transformer decoder 接收一组固定的位置编码，解码 encoder 来的信息，再通过一个 MLP 映射到 $k times 14$ 的动作序列输出。

注意位置编码，CVAE 编码器的目标动作序列、CVAE 解码器中的特征图以及 decoder 的输入都需要位置编码用于保留顺序、空间信息。此外示意图中还省略了部分 MLP，例如 14 维关节状态和 512 维 token 间的互相映射等。

*训练*。初始化编码器 $q_(phi.alt) (z mid(|) a_(t:t+k), macron(o)_t)$ 和解码器 $pi_theta (hat(a)_(t:t+k) mid(|) o_t, z)$。反复迭代如下步骤：从数据集 $cal(D)$ 采样观测 $o_t$ 和长度为 $k$ 的动作序列 $a_(t:t+k)$，从编码器基于此输出的分布中采样得到 $z$，和 $o_t$ 一起经过解码器得到预测动作序列 $hat(a)_(t:t+k)$，计算两个动作序列的 MSE 作为重建误差 $cal(L)_"reconst"$，并计算编码器 $q_(phi.alt)$ 和标准正态分布 $cal(N) (0, I)$ 之间的 KL 散度作为正则项 $cal(L)_"reg"$，得到损失函数 $cal(L) = cal(L)_"reconst" + beta cal(L)_"reg"$，用 ADAM 更新参数 $theta$ 和 $phi.alt$。

*推理*。初始化一组 FIFO 缓存 $cal(B)[0:T]$，每个 $cal(B)[t]$ 是一个 FIFO，存储不同预测动作块中 $t$ 时刻的动作。反复迭代如下步骤：令 $z = 0$，基于 $o_t$ 解码预测 $hat(a)_(t:t+k)$，将其各元素分别塞进 $cal(B)[t:t+k]$，得到当前时刻的一组动作 $A_t = cal(B)[t]$，再应用 $a_t = sum_i exp(-m i) A_t [i] \/ sum_i w_i$ 即时序集成加权得到当前动作输出。

*实验和结果*。#Cre("TODO")

#pagebreak()
== $pi_0$: A Vision-Language-Action Flow Model for General Robot Control

https://arxiv.org/pdf/2410.24164

https://blog.csdn.net/nenchoumi3119/article/details/148689279

#figure(
    placement: none,
    caption: [$pi_0$ 架构示意图],
    image("../assets/pi0_overview.png")
)

=== 预训练和后训练

类比大型语言模型和视觉语言模型中预训练模型的概念。使用来自互联网的多样化文本和语料库可以训练出较广泛的指令遵循和问题解决能力，但其完全基于抽象描述，要用到机器人学习上需要引入机器人的数据集，$pi_0$ 中用了自己录制的 $pi$ dataset 和开源的 OXE。

预训练主要得到一个 base model，可以遵循语言指令并完成一些初级的任务，不需要其对某个具体任务表现很好。要处理复杂灵巧操作任务需要在后训练中用高质量精选数据集专门微调令其适配下游任务。仅有低质量预训练数据无法处理复杂精细的任务，而只有高质量数据无法训练模型从错误（偏离训练分布）中恢复，于是结合二者。

$pi_0$ 不是最早用预训练模型针对机器人微调的，此前的类似工作主要是将动作视作一种 token，用一贯的自回归的方式预测下一个动作 token。而且动作本身是连续数值，为了用上 LLM 的训练范式就需要将其离散化成 token，比如将相近的数值视为一个 token。而 $pi_0$ 让 VLM 通过 flow matching 直接生成连续动作轨迹。不过这些纠结都是为了能在前面塞进一个 VLM 语义底座，前面的 ACT 也是直接生成动作轨迹，但没有预训练 VLM。

$pi_0$ 的预训练模型 backbone 基于 PaliGemma VLM（视觉编码器 SigLIP + 大语言模型 Gemma，别的也可以），遵循标准 late fusion VLM recipe，也就是将图片单独先编码并投影到语言模型的 token embedding space 然后再和文字 token 一起送进 LLM。用特定机器人的输入和输出增强这个 backbone。

这里 token 表示输入/输出槽位的统称，可以是离散的（语言）或连续的（图像、动作），如果都有就参考 Transfusion 对两类 token 使用不同 loss 训练。由于处理图像/文字信息和生成动作信息所需学习的模式本质不同（比如前者学习常识、语义、空间关系，后者则是学习关节状态、动力学、接触等），最好是分开使用两套权重，所以才分为了 $pi_0$ 中的 VLM 和 action expert 两部分。预训练模型参数冻结不参与训练，只训练 action expert。

#Cre("TODO") 预训练 VLM 和 action expert 在每层自注意力中通过 KV 连接，后者看到的是 VLM prefix tokens 的 K 和 V。

=== Flow Matching

动作块 $A_t = [a_t, a_(t+1), dots, a_(t+H-1)] =: a_(t:t+H)$，此处取 $H = 50$。$o_t$ 是观测，包含多个图像输入，一道语言指令，以及机器人的本体感知（proprioceptive）状态，例如 $o_t = [I_t^1, dots, I_t^n, l_t, q_t]$，$I_t^i$ 是第 $i$ 个图像，$l_t$ 是一个语言 token 序列，$q_t$ 是关节角状态向量，图像和关节状态通过编码器（官方实现中图像用 SigLIP ViT，状态直接通过一个 Linear）和线性投影（直接 Linear）映射到和语言 token 一致的 embedding 空间中。

动作块有多大，action expert 里就设有几个动作 token slot，训练时最优化如下条件 flow matching loss：

$
L^tau (theta) = EE_(p(A_t mid(|) o_t), q(A_t^tau mid(|) A_t)) norm(v_theta (A_t^tau, o_t) - u (A_t^tau mid(|) A_t))^2
$

其中 $tau in [0, 1]$ 是 flow 的时间戳，$tau -> 0$ 表示接近纯噪声，$tau -> 1$ 表示接近真实目标分布（动作）。$q (A_t^tau mid(|) A_t) = cal(N) (tau A_t, (1-tau) I)$，这是由于文献指出 flow matching 用线性高斯路径插值就可以得到很好的效果，具体地，flow 的过程是要从初始分布 $cal(N) (0, I)$ 流动到真实目标分布 $cal(N) (A_t, 0)$。实际中，先采样 $epsilon.alt ~ cal(N) (0, I)$，计算 noisy actions $A_t^tau = tau A_t + (1 - tau) epsilon.alt$，然后训练让网络输出的速度场 $v_theta (A_t^tau, o_t)$ 匹配 denoising 速度场 $u (A_t^tau mid(|) A_t) = A_t - epsilon.alt$（因为线性插值路径上每点速度就是终点 $A_t$ 减去起点 $epsilon.alt$ 除以归一化总时间 $1$，和 $tau$ 无关）。

action expert 使用完整双向注意力掩码，所有动作 token 都能互相看见。在训练时，$tau$ 采样自一个 Beta 分布，从而更关注较早（更接近噪声）的时间步。

推理时，通过将向量场 $v_theta$ 从 $tau = 0$ 积分到 $tau = 1$ 生成动作。具体地，初始化 $A_t^0 ~ cal(N) (0, I)$，和观测一起送入 action expert 得到向量场 $v_theta (A_t^0, o_t)$，然后使用 forward Euler 积分计算下一个 $A_t^delta$，依此类推直到 $A_t^1$：

$
A_t^(tau+delta) = A_t^tau + delta v_theta (A_t^tau, o_t)
$

其中 $delta$ 为积分步长，$pi_0$ 中是分十步，即 $delta = 0.1$。

注意在每步积分中观测 $o_t$ 不变化，所以可以缓存 $K_"prefix"$ 和 $V_"prefix"$ 每次用。

=== Blockwise Causal Attention Mask

$pi_0$ 中将 tokens 划分为顺序的三部分（block）：image/language（prefix）、state 和 noisy actions，每个 block 内 full bidirectional attention（互相可见），但未来 block 不可见。

第一个 block 是送进预训练 VLM 的 token，这样的设计可以确保后面的 state 和 noisy actions 不污染 VLM 给出的语义，同时这也是前面提到的在 flow integration 中可以缓存 prefix KV 的前提。

第二个 block 是 state，令其不可见 noisy actions 表示其不会受 flow integration 的影响，也可以缓存

=== 官方实现

==== 位置编码 `posemb_sincos`

位置编码将一个位置信息 `pos`（函数中是一批）转换为一个 `embedding_dim` 维度的向量，常采用正余弦位置编码，向量的每一个元素（维度）分别从不同频率的正/余弦函数中采样得到。总之位置编码要求能区分不同位置，避免不同位置相似，同时维持一些位置的特征，例如连续性等，出于各种优势，正余弦位置编码中一半来自 sin 一半来自 cos。

具体地，$pi_0$ 中先按给定的周期范围生成待采样的周期序列，一半 sin 一半 cos，指数间隔：

```
fraction = jnp.linspace(0.0, 1.0, embedding_dim // 2)
period = min_period * (max_period / min_period) ** fraction
```

然后 `jnp.einsum("i,j->ij", pos, freq)` 等价于外积 `pos[:, None] * freq[None, :]`，外面套上 sin 和 cos，最后直接拼接在一起。

==== #Cre("TODO")

详见代码注释

=== 实验和结果

#Cre("TODO")

#pagebreak()
== Universal Manipulation Interface: In-The-Wild Robot Teaching Without In-The-Wild Robots

https://zhuanlan.zhihu.com/p/2003882878099014705

数据采集常用手持夹爪等传感器设备，虽然实验者可以用它们做出任意动作，但很多并不能迁移成有效的机器人策略。具体地，人类太复杂精细的动作机器人不见得学得会，这导致此前很多实验为了保证 transferability，所用的数据大多还是刻意简单化的动作，比如抓取、放下等。可以说数据中 visual diversity 很高（覆盖了很多视觉场景），但 action diversity 还是较低（能学到的技能有限）。

文章讨论了一些以往制约动作迁移的因素：

1、视觉上下文不充分。比如腕部摄像头虽然方便，但是缺少全局信息，还有部分视野被夹爪挡住，都会影响动作学习。

2、动作不精确。许多手持设备用单目视觉估计动作，但也受限于其模糊的尺度、运动模糊、纹理缺乏等，难以精确估计。

3、延迟差异。训练阶段所用的采集的数据中，观测和动作可以对齐同步，由此学到的策略就难以处理推理时来自传感器、推理、执行器等的延迟，造成 out-of-distribution inputs。

4、不充分的策略表示。policy representations 也就是参数化策略函数的方式，比如 MLP 等。简单的模型无法表达人类复杂的多模态动作，比如同场景下人类可能有很多种合理动作。特别是在更大规模收集的数据上，这样的模型无法准确学习行为模式，可能学出一种多半不好使的 “平均” 情况，参见 action chunking。通常会考虑学习随机策略的概率分布。

文章也针对这些问题提出一些修改设计：

1、改进物理接口，使其更适合人类示范采集（人类操作自然，同时也能采集足够信息供策略学习），文章用鱼眼摄像头提供更大的 field of view 及视觉上下文，在两侧放了两个反射镜提供隐式的双目信息。

2、设计更好的 policy interface，即如何表示输入的观测和输出的动作。策略接口最好是硬件无关的（hardware-agnostic），即不和机器人状态空间强绑定。此外加入 inference-time latency matching（推理时延匹配）。还有使用 relative trajectory，不预测绝对位置而是预测相对当前位置的运动，更容易迁移。

最后这套系统称为 universal manipulation interface（UMI）。

#Cre("TODO")

#pagebreak()
== FAST: Efficient Action Tokenization for Vision-Language-Action Models

https://zhuanlan.zhihu.com/p/1910755399646287695

$pi_0$ 是基于 flow matching 的 diffusion VLA，常见的还有一类是基于自回归（AR）的 VLA，即和原本 transformer 类似的预测 next token 的形式。但由于动作信号是连续的，需要将动作化为离散的 token。最简单的方法是逐维度、逐时间步分 bin，举个例子，长度 $H$ 的 action chunk，每个 action 有 $D$ 维，可以把每个维度的连续数值归一化到 $[-1, 1]$，分成 $256$ 个离散区间，每个用一个离散 token 表达；要预测一个 chunk，就要预测 $H times D$ 次 next token 然后拼起来，这样做问题不少：速度慢和量化误差是必然有的，此外同一个 action 不同维度被分开预测破坏了其间的耦合，而且这样硬拆散的操作让生成过程看不到完整轨迹的结构，丢失了很多性质，总之在处理高频控制灵巧操作时表现很差。

对此，可以寻找更好的 tokenization 方法，将冗余的动作信号*压缩*成较少数量的高信息密度的 token。FAST 是频率空间动作序列 token 化（frequency-space action sequence tokenization），能大幅加快训练速度。

=== Action Representations for VLA Training

#Cre("TODO")

=== FAST

*问题陈述*。模型需要训练的是策略 $pi (a_(1:H) mid(|) o)$，即将观测 $o$ 映射到一组未来动作序列 $a_(1:H)$ 上，这里假设采用 action chunking，输出的是一组长度 $H$ 的动作序列。action tokenization 的目标是定义一组映射 $cal(T)_a : a_(1:H) -> [T_1, dots, T_n]$，将一串 $abs(cal(A))$ 维的连续动作映射到一组 $n$ 个离散 tokens $T in abs(cal(V))$ 上。

*分词算法如何影响 VLA 训练*。文章举了个例子，创建了一个数据集，目标是学习预测由四个随机点插值来的三次样条曲线（cubic spline）。实验主要对比 naive tokenization（binning）和 DCT tokenization（基于离散余弦变换将序列转换到频域，然后对其系数进行离散化），调整的参数是采样率（也就是一下输出的动作块的长度），动作序列越长越能适应高频控制/数据集，理论上来说动作可以更平滑。实验结果来看，naive 方法随着采样率提高表现迅速下滑，主要问题在于随着采样率提高相邻 token 之间变得非常相似，marginal information 和与之关联的学习信号强度趋于零，极大降低训练收敛速度，直至完全无法在合理时间内收敛。

*基于 DCT 的时间序列压缩*。文章称任何足够有效的压缩方法应用在动作目标上都可以提高 VLA 的训练速度。DCT 将连续信号分解为一组不同频率余弦函数的和，相比基于学习的压缩算法，DCT 这类分析算法更简单且高效。

*FAST*。首先把动作数值归一化，方式是对训练集中每个动作维度，将其第 1 和第 99 个分位数映射到 $[-1, dots, 1]$。然后对每个动作维度做 DCT 得到一系列系数（低频到高频），在这步可以通过 scale-and-round 操作直接忽略部分不显著的系数，程度通过超参数控制。一个动作维度算一行，把系数都拼起来得到频率矩阵，并且 rounding 过后的矩阵通常是稀疏的，接下来就是要将其转为更稠密的 action tokens，现在还太长。首先将其 flatten 为一个序列，方式是一列一列数，即低频成分排在前面（先低频控制形状效果好一些）。接下来用一个训练的 BPE（byte pair encoding）tokenizer 将其无损压缩成一串 action token，合理的 BPE 可以将零元素以及经常出现 tokens 组合压缩掉。输出不一定固定长度，需要用常用方法例如 padding 应对一下。

#pagebreak()
== GR00T N1: An Open Foundation Model for Generalist Humanoid Robots

#Cre("TODO")

#pagebreak()
== OpenVLA: An Open-Source Vision-Language-Action Model

#Cre("TODO")

#pagebreak()
== $pi_(0.5)$: A Vision-Language-Action Model with Open-World Generalization

#Cre("TODO")
