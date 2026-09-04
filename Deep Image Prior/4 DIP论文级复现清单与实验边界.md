---
tags:
  - 暑假学习
  - Deep-Image-Prior
  - 复现
  - PyTorch
  - 实验设计
created: 2026-09-03
updated: 2026-09-03
---

# 4 DIP 论文级复现清单与实验边界

> [!cite] Zotero 资料基线
> [[Ulyanov 等 - 2018 - Deep Image Prior|Ulyanov et al. 2018 — Deep Image Prior]] · [在 Zotero 中打开](zotero://select/library/items/DNFAYSCE)
> 
> 本笔记整理的是论文明确给出的实现信息，不补充论文之外的 PyTorch API 教程。

## 最小复现闭环

DIP 复现的核心流程是：

- 准备单张退化观测 $x_0$；
- 根据任务定义可微数据项 $E$；
- 建立随机初始化的 fully-convolutional generator $f_\theta$；
- 固定随机输入 $z$；
- 反复计算 $x=f_\theta(z)$；
- 通过任务观测模型计算 $E(x;x_0)$；
- 反向传播只更新 $\theta$；
- 保存中间输出，而不是只保存最终最低 loss 的输出；
- 对需要 early stopping 的任务，在过拟合噪声之前选择结果。

## 论文默认网络信息

论文主要采用 encoder-decoder / `hourglass` 网络：

- fully convolutional；
- 多层 downsampling / upsampling；
- 可带 skip connections；
- 激活使用 LeakyReLU；
- downsampling 主要通过 stride convolution；
- upsampling 在 bilinear 与 nearest neighbor 中选择；
- 作者尝试 transposed convolution，但论文报告结果更差；
- 除少数特定实验外，卷积采用 reflection padding。

这些设计不是普通工程细节，因为 [[1 网络重参数化为什么能形成先验]] 中的隐式 prior 正是由架构产生的。

## 输入 $z$

论文主要考察两种输入：

- `random`：随机张量；
- `meshgrid`：二维坐标网格，在大孔洞 inpainting 中作为额外平滑先验。

多数任务中使用随机输入，并保持 $z$固定；作者没有默认同时优化 $z$。

论文的 super-resolution 默认配置使用

$$z\in\mathbb R^{32\times W\times H},\qquad z\sim U(0,0.1).$$

训练过程中还常使用输入扰动：每轮在 $z$上加零均值 Gaussian noise。

## 论文 super-resolution 默认参数

MinerU 技术细节给出的默认 SR 配置为：

- 输入通道：32；
- 五级 down / up；
- $n_u=n_d=[128,128,128,128,128]$；
- $k_u=k_d=[3,3,3,3,3]$；
- skip 通道 $n_s=[4,4,4,4,4]$；
- skip kernel $k_s=[1,1,1,1,1]$；
- 输入扰动标准差 $\sigma_p=1/30$；
- `num_iter = 2000`；
- learning rate `0.01`；
- upsampling：bilinear。

对 8× SR，论文把

$$\sigma_p=\frac1{20}$$

并把迭代次数改为 4000。

SR 的 forward operator $d$由 Lanczos2 low-pass filtering 与 resampling 构成，并作为固定可微层使用。

这些数值属于论文实验配置，不应被写成所有 DIP 任务的普适最佳参数。

## Inpainting 的论文调整

- `text inpainting`：沿用 SR 的大部分超参数，迭代 6000 次；
- `large-hole inpainting`：采用 meshgrid 输入，去掉 skip connections，迭代 5000 次。

这再次说明架构和输入都是 prior 的组成部分，而不是只需固定一个网络模板。

## Denoising 的论文调整

论文中一组 denoising 实验使用与 SR 相近的超参数，但迭代数约为 1800；并使用 exponential sliding window 对最近输出进行平均，以改善恢复结果。

这里不能简单得出“1800 就是正确 early-stop 点”。这个数字来自论文特定数据和设置，真正需要观察的是拟合轨迹中的 signal-first / noise-later 现象。

## 优化器与框架

论文明确使用：

- `Adam` optimizer；
- `PyTorch` framework。

作者指出，每张图都需要重复深网 forward/backward，因此运行时间通常为数分钟。这也是 DIP 相比预训练 feed-forward 网络的主要实用缺点之一。

## loss destabilization 的保护

论文观察到：当 loss 下降到一定区域时，优化有时会突然 destabilize，表现为 loss 大幅增加和输出模糊，之后又可能重新下降。

作者采用的处理是：

> 若相邻迭代 loss 的变化超过某个阈值，就回滚到上一轮参数。

因此复现时至少应保存最近的稳定参数，而不是只维护当前权重。

## 建议记录的实验日志

下面属于基于论文机制设计的复现记录格式，不是论文强制格式：

- iteration；
- data loss；
- 当前输出图；
- 若有 ground truth，记录 PSNR / SSIM；
- 输入扰动强度 $\sigma_p$；
- learning rate；
- 网络结构配置；
- 当前是否发生 loss spike；
- 最佳视觉 / 指标迭代位置。

这样才能真正观察 [[2 noise impedance与early stopping]]，而不是只得到一个最终图片。

## 复现完成标准

第一轮不必追求完全复制论文所有表格，至少完成：

- 一个去噪任务：看到中期恢复后、长期过拟合噪声的轨迹；
- 一个超分辨率任务：严格把 downsampling operator 写入 loss，而不是直接拿 HR ground truth 监督；
- 一个 inpainting 任务：loss 只在 mask 已知区域计算；
- 保存并比较多个迭代时刻的输出；
- 能解释为什么网络从随机初始化也能形成 prior。

只有这几项都验证，才算真正复现了 DIP 的核心机制，而不是“用 PyTorch 跑了一个 encoder-decoder”。

## 论文边界

- 论文没有提供“所有任务统一最优”的网络与超参数。
- DIP 很慢，且作者明确承认专门训练的 feed-forward 方法在实际应用中通常更快、更强。
- 不同架构的 prior 强弱不同；skip connection 的作用依任务而异。
- 高容量网络最终能拟合噪声，不能把训练到完全收敛当作默认策略。
- 本笔记没有引入论文外的 PyTorch API、环境版本或第三方实现细节；后续真正在本机复现时，应把“代码工程问题”和“论文方法事实”分开记录。