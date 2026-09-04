---
tags:
  - 暑假学习
  - PnP-FISTA
  - FISTA
  - Plug-and-Play
  - 逆问题
created: 2026-09-03
updated: 2026-09-03
---

# 0 PnP-FISTA（讲解版）

> [!cite] Zotero 资料基线
> - [[Gavaskar 等 - 2021 - On Plug-and-Play Regularization using Linear Denoisers|Gavaskar et al. 2021 — Linear Denoisers PnP]] · [Zotero](zotero://select/library/items/R3JAW5G3)
> - [[Sinha 等 - 2024 - FISTA Iterates Converge Linearly for Denoiser-Driven Regularization|Sinha & Chaudhury 2024 — PnP-FISTA]] · [Zotero](zotero://select/library/items/IJAI4UMU)
> - FISTA 前置回顾：[[1 FISTA（讲解版）]]；PnP-ADMM 对照：[[0 PnP-ADMM（讲解版）]]。
> - 本章数学事实只按上述 Zotero 条目及其 MinerU Markdown 核对。

PnP-FISTA 最容易被一句“把 FISTA 的 prox 换成 denoiser”讲得过于简单。算法形式确实只改一行，但真正需要理解的是：换掉 proximal operator 以后，原本 FISTA 的显式优化目标、$O(1/k^2)$ 目标值保证以及迭代点是否收敛，都必须重新说明。

## 从 FISTA 的结构看，PnP 替换为什么自然

经典 FISTA 解决

$$\min_x F(x)=f(x)+g(x),$$

其中 $f$ 光滑，$g$ 可以不可微但 proximal operator 易算。其核心更新可以理解为

$$z_k=y_k-\gamma\nabla f(y_k),$$

$$x_k=\operatorname{prox}_{\gamma g}(z_k),$$

再用 Nesterov 外推构造下一轮 $y_{k+1}$。

因此 FISTA 本来就有两个接口：

- `数据一致性接口`：梯度步只使用 $f$；
- `结构先验接口`：proximal step 只使用 $g$。

PnP-FISTA 直接把第二个接口改为 denoiser：

$$x_k=\mathcal D(z_k).$$

于是它和 [[3 PnP替换与denoiser|PnP-ADMM 的 PnP 替换]]有同一个思想：把显式先验计算替换成更强的图像去噪模块。区别在于父算法不同——PnP-FISTA 仍保留梯度步和 Nesterov 外推，而不是 ADMM 的变量分裂和对偶更新。完整迭代见 [[1 从FISTA到PnP-FISTA]]。

## 第一层理论问题：这个 denoiser 还是不是某个 prox

如果 denoiser 真的是某个凸函数 $g$ 的 proximal mapping，那么 PnP-FISTA 其实可以重新解释成普通 FISTA，经典优化理论就有机会接回来。

对对称线性 denoiser，这个问题早已有部分答案；但实际 kernel denoiser，例如标准 NLM 的归一化权重矩阵，往往不是对称矩阵。Gavaskar 等人 2021 的关键推进是：

> 非对称线性 denoiser 不一定是标准 Euclidean prox，但可能是某个加权几何中的 `scaled proximal map`。

令 $H\succ0$，定义

$$\lVert x\rVert_H=\sqrt{x^THx}.$$

相应的 scaled proximal operator 是

$$\operatorname{prox}_{g,\lVert\cdot\rVert_H}(y)=\operatorname*{argmin}_x\left\{\frac12\lVert x-y\rVert_H^2+g(x)\right\}.$$

论文证明，一类线性 denoiser $W$ 在适当条件下可写成这种 scaled prox；特别地，kernel denoiser $W=D^{-1}K$ 可以使用 $H=D$。这样，“非对称”不再意味着一定没有凸正则化解释，而是说明正确的几何可能不是标准 Euclidean 几何。详见 [[2 线性denoiser与scaled proximal map]]。

## 第二层：既然 prox 的几何变了，算法本身也必须变

这是最重要的一步。

不能只证明 $W$ 是 $H$-scaled prox，然后仍然运行原来的 Euclidean PnP-FISTA。Gavaskar 等人专门给出反例说明：对非对称 denoiser，标准 PnP 算法可能不能继承父算法的收敛性质。

因此 scaled PnP-FISTA 把梯度步也改到同一个 $H$ 几何中：

$$x_{k+1}=\mathcal D\left(y_{k+1}-\rho^{-1}H^{-1}\nabla f(y_{k+1})\right).$$

Nesterov 的 $t_k$ 与外推形式保持不变，但普通梯度 $\nabla f$ 前多了 $H^{-1}$。这不是随意的预条件，而是为了让“梯度步 + denoiser prox”在同一个内积空间中组成真正的 FISTA。

在论文条件下，scaled PnP-FISTA 重新对应一个显式凸目标 $f+g$，并得到

$$f(x_k)+g(x_k)\le p^\star+O(1/k^2).$$

这恢复的是`目标值收敛率`。若把 FISTA 的动量序列改为 Chambolle–Dossal 型 $t_{k+1}=1+k/a$、$a>2$，论文还可得到迭代点收敛到某个 minimizer，但这一结果本身没有宣称线性速率。详见 [[3 scaled PnP-FISTA与目标值收敛]]。

## 第三层：$O(1/k^2)$ 不是“迭代点线性收敛”

这两个结论很容易混在一起。

经典 FISTA 的 $O(1/k^2)$ 说的是

$$F(x_k)-F^\star,$$

也就是目标函数值误差。

Sinha 与 Chaudhury 2024 研究的是更直接的问题：

$$\lVert x_k-x^\star\rVert,$$

也就是图像迭代本身离极限点还有多远。

他们把范围收窄到线性逆问题

$$f(x)=\frac12\lVert Ax-b\rVert_2^2$$

和一类可验证的线性 denoiser $W$。在对称情形下，关键条件包括：

$$\sigma(W)\subset[0,1],$$

以及

$$\ker(A)\cap\operatorname{fix}(W)=\{0\}.$$

第二条很有直觉：不能存在某个非零方向既完全看不见于测量系统 $A$，又完全不被 denoiser 改变。否则算法在那个方向上缺乏收缩来源。

对 PnP-FISTA 的步长，论文要求

$$0<\gamma<\frac{1}{\lambda_{\max}(A^TA)}.$$

在这些条件以及动量系数 $\alpha_k\to1$ 下，作者把二阶带动量递推提升到 $\mathbb R^{2n}$ 的状态空间，证明极限转移矩阵的 spectral radius 严格小于 1，因此存在 $0\le\beta<1$ 使

$$\lVert x_k-x^\star\rVert\le C\beta^k$$

从某个迭代开始成立。这才是`线性迭代收敛`。详见 [[4 PnP-FISTA的线性迭代收敛与谱分析]]。

## 为什么谱半径会成为核心

当 $f$ 是二次损失、denoiser 是固定线性算子时，梯度步和 denoiser 组合本身就是一个 affine map。加上 FISTA 动量以后，只需把

$$z_k=\begin{bmatrix}x_k\\x_{k-1}\end{bmatrix}$$

当成新状态，就能写成

$$z_k=R_kz_{k-1}+s.$$

随着 $\alpha_k\to1$，$R_k$趋于某个固定矩阵 $R_\infty$。若

$$\rho(R_\infty)<1,$$

误差就会被反复乘上严格收缩的线性算子，从而产生几何级数式衰减。这就是 2024 论文能给出比一般 FISTA 更强结论的原因：问题被限制到足够线性的范围后，整个算法可以直接做谱分析。

## 非对称 kernel denoiser 仍然要回到 scaled 几何

Sinha 2024 也把结果扩展到非对称 NLM 型 kernel denoiser。做法与 2021 论文一致：在 $D$ 加权内积中运行 scaled PnP-FISTA，梯度变为

$$\nabla_D f(y)=D^{-1}A^T(Ay-b).$$

相应步长条件变成

$$0<\gamma<\frac{1}{\lambda_{\max}(D^{-1/2}A^TAD^{-1/2})}.$$

因此“scaled”不是某篇论文的技巧，而是一条贯穿线性非对称 denoiser 理论的主线：换了 prox 的几何，就要同步换梯度和步长的几何。

## 加速不等于迭代点一定更快

这里还有一个对 FISTA 很重要的修正认识。

Nesterov 动量保证经典 FISTA 有更快的`目标值最坏情形速率`，但这不意味着每个问题里 $x_k$ 到 $x^\star$ 的距离都比 ISTA 更快。Sinha 2024 的实验专门比较不同 $\alpha_k$，发现某些动量选择下 PnP-FISTA 的迭代点收敛甚至可以比无动量 PnP-ISTA 慢。

所以以后看到“FISTA 加速”必须问清楚：

- 加速的是目标值还是迭代点？
- 是最坏情形理论还是具体谱半径？
- denoiser 是否让父算法的目标函数解释仍成立？

## 本章收束

PnP-FISTA 的理论链可以压缩为：

`FISTA 的 prox 接口 → denoiser 替换 → 一般最优性解释丢失 → 线性 denoiser可在加权几何中恢复 scaled prox → 算法也改为 scaled PnP-FISTA → 恢复显式目标和 O(1/k²) 目标值收敛 → 在线性逆问题与更强谱条件下进一步证明迭代点全局线性收敛。`

这也解释了为什么阶段 4 不能简单记成“PnP-ADMM 换成 FISTA 更快”：PnP-FISTA 对 $f$ 的光滑性要求更强，理论也更依赖 denoiser 结构；它的优势是每轮结构简单并保留动量机制，但能声称什么收敛结论必须严格匹配所用 denoiser 与 forward model。