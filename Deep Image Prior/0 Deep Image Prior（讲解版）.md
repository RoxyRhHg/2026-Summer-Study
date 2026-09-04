---
tags:
  - 暑假学习
  - Deep-Image-Prior
  - DIP
  - 逆问题
  - 图像恢复
created: 2026-09-03
updated: 2026-09-03
---

# 0 Deep Image Prior（讲解版）

> [!cite] Zotero 资料基线
> [[Ulyanov 等 - 2018 - Deep Image Prior|Ulyanov et al. 2018 — Deep Image Prior]] · [在 Zotero 中打开](zotero://select/library/items/DNFAYSCE)
> 
> 本章只使用该 Zotero 条目及其 MinerU Markdown。论文明确使用 PyTorch 与 Adam，但本章不额外引入外部 PyTorch 教程。

Deep Image Prior 最反直觉的一点是：`网络没有在数据集上训练过，网络结构本身仍然可以充当图像先验。`

它不是先训练一个网络，再拿网络去恢复新图像；而是对每一张待恢复图像，重新随机初始化一个生成网络，并只利用当前这张退化图像优化网络参数。

## 从显式正则化到网络重参数化

传统图像恢复常写成

$$x^*=\operatorname*{argmin}_x E(x;x_0)+R(x),$$

其中：

- $x_0$：观测到的退化图像；
- $x$：待恢复图像；
- $E(x;x_0)$：由任务决定的数据一致性项；
- $R(x)$：显式图像先验，例如 TV。

DIP 的关键动作不是把 $R(x)$换成另一个显式公式，而是直接把图像写成网络输出：

$$x=f_\theta(z).$$

其中 $z$通常是固定的随机张量，$\theta$是网络权重。于是优化变量从像素 $x$变成网络参数 $\theta$：

$$\theta^*=\operatorname*{argmin}_\theta E(f_\theta(z);x_0),\qquad x^*=f_{\theta^*}(z).$$

表面上看，显式正则项 $R(x)$消失了；实际上，允许被优化器快速访问到的图像集合被网络结构和优化路径重新限制。详细解释见 [[1 网络重参数化为什么能形成先验]]。

## 先验到底藏在哪里

DIP 不声称随机网络不能表示噪声。相反，高容量网络最终也可能把随机噪声拟合得很好。

论文观察到的关键差异是`拟合速度`：

- 自然图像和自然图像的大尺度结构被卷积生成网络较快拟合；
- 随机打乱像素、白噪声等结构被拟合得明显更慢。

作者把这种现象描述为网络参数化对自然信号具有较低阻抗、对噪声具有较高阻抗。

因此 DIP 的先验不仅存在于“哪些图像能表示”，更存在于：

> `梯度优化从随机初始化出发时，会按什么顺序经过图像空间。`

网络结构改变了优化轨迹，使它往往先恢复低频、大尺度和自相似结构，之后才逐渐拟合高频噪声。见 [[2 noise impedance与early stopping]]。

## 为什么 early stopping 是方法的一部分

对超分辨率这类欠定问题，优化可能停在一个满足数据一致性的合理图像附近；但对去噪，若目标只是

$$\min_\theta\lVert f_\theta(z)-x_0\rVert_2^2,$$

而 $x_0$本身包含噪声，那么长期训练的真正极限就是把噪声也拟合进去。

因此去噪中的良好结果往往出现在优化轨迹中间：

`先拟合自然图像结构 → 经过较好的重建 → 继续优化 → 开始拟合噪声。`

所以 early stopping 不是纯粹为了节省计算，而是 DIP 隐式正则化机制的重要组成部分。

## 网络结构为什么真的重要

论文主要采用 fully-convolutional 的 encoder-decoder / hourglass 结构，常配合 skip connections。卷积、上采样与多尺度结构会诱导空间平稳性、自相似性以及局部相关性。

但“skip 越多越好”并不成立。论文不同任务中观察到：

- 多尺度 skip 结构能形成丰富的多尺度纹理；
- 某些过宽的 skip connection 会让网络拟合过快，使先验过弱；
- 大孔洞 inpainting 中，论文甚至去掉 skip connections；
- 架构太强或太弱都会改变优化速度和恢复偏好。

因此 DIP 中网络架构本身就是正则化设计的一部分，而不只是实现容器。

## 三个核心任务如何统一

### 去噪

直接用

$$E(x;x_0)=\lVert x-x_0\rVert_2^2,$$

代入 $x=f_\theta(z)$。网络早期先拟合自然结构，之后逐渐拟合噪声，所以需要 early stopping。

### 超分辨率

令 $d(\cdot)$为已知降采样算子：

$$E(x;x_0)=\lVert d(x)-x_0\rVert_2^2.$$

DIP 通过网络参数化在无穷多个满足低分辨率观测的高分辨率候选中偏向自然图像结构。论文强调，方法只使用当前 LR 图像，没有外部训练集。

### Inpainting

令 $m$为已知像素的二值 mask：

$$E(x;x_0)=\lVert(x-x_0)\odot m\rVert_2^2.$$

数据项完全不约束缺失区域，缺失部分之所以能被填出来，主要依赖网络结构对空间纹理和上下文的隐式偏好。

三个任务只需更换数据一致性项，网络先验机制保持不变，详见 [[3 DIP在去噪超分辨率与修复中的统一写法]]。

## 它与 PnP 的区别

DIP 和 [[0 PnP-ADMM（讲解版）|PnP-ADMM]]、[[0 PnP-FISTA（讲解版）|PnP-FISTA]]都属于“弱化显式手工正则化”的思路，但方式完全不同：

| 方法 | 先验放在哪里 | 每轮优化什么 |
| --- | --- | --- |
| 经典正则化 | 显式 $R(x)$ 或 $g(x)$ | 图像变量 $x$ |
| PnP | denoiser 算子 $\mathcal D$ | 图像状态 + 算法状态 |
| DIP | 生成网络参数化 $x=f_\theta(z)$ | 网络参数 $\theta$ |

PnP 的 denoiser通常是预先设计或训练好的模块；DIP 原论文的核心恰恰是不需要外部训练数据，针对单张图像从随机权重开始拟合。

## 论文级实现信息

论文明确说明：

- 主要使用 fully-convolutional hourglass / encoder-decoder；
- 常用 LeakyReLU；
- downsampling 多用 stride convolution；
- upsampling 采用 bilinear 或 nearest neighbor，论文中 transposed convolution 表现较差；
- 大多数卷积使用 reflection padding；
- $z$常由 $U(0,0.1)$随机初始化并固定；
- 训练时常对 $z$叠加小的 Gaussian perturbation；
- 使用 Adam optimizer；
- 框架为 PyTorch；
- 单图需要反复 forward/backward，论文明确指出耗时可达数分钟。

可直接复现的论文参数整理见 [[4 DIP论文级复现清单与实验边界]]。

## 本章最重要的边界

- DIP 不是“随机网络天然输出好图”，而是`随机初始化网络 + 特定架构 + 针对单图优化轨迹`共同形成先验。
- 网络最终能够拟合噪声，所以不能把“训练 loss 更低”简单等同于“恢复更好”。
- early stopping 对去噪等任务尤其关键。
- 论文结果不能证明 DIP 在所有任务上优于训练式模型；作者自己指出它速度慢，而且专门训练的方法在实际应用中通常更强更快。
- DIP 的“prior”是隐式的，不对应论文中某个写出的概率密度或显式 $R(x)$。

整条理解链可以压缩为：

`显式正则化 x → 用生成网络重参数化 x=fθ(z) → 网络结构改变可达路径与拟合速度 → 自然结构先于噪声被拟合 → 任务数据项 + 架构偏置 + early stopping 共同完成单图恢复。`
