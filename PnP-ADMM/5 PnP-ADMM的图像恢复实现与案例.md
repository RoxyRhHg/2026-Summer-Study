---
tags:
  - 暑假学习
  - PnP-ADMM
  - 图像恢复
  - 超分辨率
  - 单光子成像
created: 2026-09-03
updated: 2026-09-03
---

# 5 PnP-ADMM 的图像恢复实现与案例

> [!cite] Zotero 资料基线
> [[Chan 等 - 2017 - Plug-and-Play ADMM for Image Restoration - Fixed Point Convergence and Applications|Chan et al. 2017 — PnP-ADMM]] · [在 Zotero 中打开](zotero://select/library/items/2HLFYW77)

## PnP 的工程瓶颈其实在数据一致性步

[[3 PnP替换与denoiser]] 把先验步改造成了可替换的 denoiser，但这并不意味着整个逆问题都变成“调用一次去噪器”这么简单。

PnP-ADMM 仍然需要反复求解

$$x^{k+1}\in\operatorname*{argmin}_x\left[f(x)+\frac{\rho_k}{2}\lVert x-(v^k-u^k)\rVert_2^2\right].$$

其中 $f$ 完整保留了成像系统、观测噪声和数据一致性。换句话说：

- `v 步`因为 PnP 而获得通用去噪接口；
- `x 步`仍然必须针对具体 forward model 设计高效求解器。

Chan 等人的论文专门选择超分辨率与单光子成像，展示两种不同的加速思路：前者利用线性算子的特殊结构，后者利用目标函数的可分性。

## 图像超分辨率

### 前向模型

论文把超分辨率写成抗混叠模糊后再降采样：

$$f(x)=\lVert SHx-y\rVert_2^2,$$

其中：

- $x$：待恢复的高分辨率图像；
- $H$：循环卷积矩阵，对应抗混叠模糊核；
- $S$：降采样矩阵；
- $y$：低分辨率观测。

记 $G=SH$，则 PnP-ADMM 的 $x$ 子问题是一个带二次近端项的最小二乘问题：

$$\min_x\ \lVert Gx-y\rVert_2^2+\frac{\rho}{2}\lVert x-\widetilde x\rVert_2^2,$$

其中 $\widetilde x$ 由当前 $v,u$ 状态决定。该问题等价于解一个线性系统；真正的计算难点是怎样避免每轮都做高维矩阵求逆或内层迭代。

### 两个容易的特殊情况

- `非盲去模糊`：$S=I$。因为 $H$ 是循环卷积，可由离散 Fourier 变换对角化，所以 $x$ 步能够用 FFT 快速实现。
- `插值`：$H=I$。此时 $S^TS$ 是由 0/1 构成的对角结构，更新可以转成逐元素运算。

这两个例子说明：PnP 的数据步是否快，取决于 forward operator 有没有可利用的代数结构。

### 一般超分辨率为什么更难

当 $G=SH$ 同时包含卷积与降采样时，$H^TS^TSH$通常不能像纯卷积那样直接被 FFT 对角化。常见做法是：

- 再引入额外变量分裂，但会增加对偶变量与参数；
- 在每次外层 PnP-ADMM 中再运行 conjugate gradient，但会产生内层迭代成本。

论文选择另一条路线：先利用 Sherman–Morrison–Woodbury 恒等式，把原来的高维逆转化为观测空间中的较小逆；再研究

$$GG^T=SHH^TS^T.$$

在标准 $K$ 倍均匀降采样、$H$ 为循环卷积的条件下，论文使用 `polyphase decomposition` 分析“上采样 → 滤波 → 下采样”结构，最终把 $GG^T$等效成可预计算的有限冲激响应滤波器。这样就能重新利用 FFT，避免每轮运行 conjugate gradient。

这里需要保留三个适用边界：

- 结论依赖标准均匀降采样结构；
- $H$需要满足论文使用的循环卷积结构；
- 周期边界条件下论文给出的闭式实现是精确的；其他边界需要另外处理。

因此 `polyphase` 不是 PnP-ADMM 的通用组成部分，而是针对这个 forward model 的快速 $x$ 步实现。

## 单光子成像：另一种加速来自可分性

论文第二个例子是 Quanta Image Sensor（QIS）单光子成像。这里不是普通的“图像 + 加性 Gaussian 噪声”，而是：

`潜在强度 x → 光子到达率 → Poisson 光子计数 → 阈值量化 → 二值观测 y`。

论文考虑空间过采样的 QIS：每个传统像素区域由多个 jot 进行单光子探测，并由二值输出反推潜在图像。数据项 $f(x)$来自该量化 Poisson 观测模型的负对数似然。

代入 PnP-ADMM 后，关键结构不是 FFT，而是 `可分离`：$x$ 子问题可以按像素 $x_j$独立拆开。于是高维优化变成许多一维问题，每个像素只需求一个标量方程的根。

论文进一步指出，可以把与观测计数、$\rho$及当前近端中心有关的一维解预先组织成 lookup table，从而把迭代中的求解成本降下来。

这个案例说明了与超分辨率不同的另一条原则：

> PnP 的数据一致性步不一定靠矩阵对角化加速；如果 likelihood 本身可分，也可以把高维问题拆成大量独立低维问题。

## 实验结果应该怎样理解

论文实验统一采用 BM3D 作为 denoiser，重点不是比较不同 denoiser，而是验证 continuation PnP-ADMM 与快速数据步实现。

### 超分辨率实验

- 使用 10 张标准测试灰度图像；
- 包括不同降采样倍数、bicubic 或 Gaussian 抗混叠模型，以及有无噪声的配置；
- 在部分配置中，PnP 方法的 PSNR 优于比较方法；
- 但论文也明确展示了 `model mismatch` 的影响：当竞争方法只按 bicubic 模型实现，而真实退化采用 Gaussian 模糊时，知道真实 forward model 的方法天然占优。

所以不能把表格简单总结成“PnP 对所有超分辨率算法都更强”。更可靠的结论是：PnP 能把较强的通用 denoiser 与精确的数据模型结合，而 forward model 是否匹配会显著影响结果。

### QIS 实验

- 测试不同空间过采样规模；
- 与最大似然估计和 TV-ADMM 进行比较；
- 每个配置使用多次独立 photon-arrival Monte Carlo 试验；
- 论文报告的 PnP 结果在其 QIS 设置下取得更高的平均 PSNR，并在示例图中表现出更低的可见噪声。

这些结果只支撑论文中的特定 QIS 模型、参数和比较设置，不应泛化为“PnP 对所有 Poisson 逆问题都优于传统方法”。

## 从两个案例提炼出的实现原则

- `PnP 只把先验接口模块化`：数据一致性仍必须认真建模。
- `先看 forward model 的结构`：循环卷积适合 FFT；采样结构可能适合 polyphase；似然可分则适合逐变量求解。
- `避免无必要的 inner solver`：如果能找到闭式、FFT 或 lookup table，就不要把每个 PnP 外层迭代再套一个昂贵的迭代优化器。
- `数据模型匹配很重要`：强 denoiser 不能弥补错误的 forward model。
- `理论与实现要分开`：[[4 Bounded denoiser与fixed-point convergence]] 讨论为什么迭代稳定；本笔记讨论每一步怎样算得快。两者共同决定一个 PnP 方法是否真正可用。
