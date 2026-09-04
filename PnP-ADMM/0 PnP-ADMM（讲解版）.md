---
tags:
  - 暑假学习
  - PnP-ADMM
  - 讲解版
  - 图像恢复
created: 2026-09-03
updated: 2026-09-03
---

# 0 PnP-ADMM（讲解版）

> [!cite] Zotero 资料基线
> - [[Boyd 等 - 2011 - Distributed Optimization and Statistical Learning via the Alternating Direction Method of Multipliers|Boyd et al. 2011 — ADMM]] · [在 Zotero 中打开](zotero://select/library/items/B7AYLH4I)
> - [[Chan 等 - 2017 - Plug-and-Play ADMM for Image Restoration - Fixed Point Convergence and Applications|Chan et al. 2017 — PnP-ADMM]] · [在 Zotero 中打开](zotero://select/library/items/2HLFYW77)
> - 本章数学事实只按上述 Zotero 条目及其 MinerU Markdown 核对；已有学习笔记仅用于排版与双链组织。

PnP-ADMM 的主线可以压缩成一句话：保留 ADMM 的数据一致性与对偶协调机制，把显式先验的 proximal 更新换成现成 denoiser；因此模块更自由，但收敛终点的含义也随之改变。

## 从一个难拆的目标到两个可交换模块

图像恢复先从 MAP 形式出发：给定观测$y$，求潜在图像$x$，

$$\widehat x\in\operatorname*{argmin}_x f(x)+\lambda g(x).$$

$f$负责回答“候选图像经过前向模型后是否解释观测”，$g$负责回答“候选图像是否像合理的图像”。困难在于，这两种结构混在同一个变量上：数据项可能便于线性代数求解，先验项却可能适合另一种专门算法。

ADMM 引入副本$v$：

$$\min_{x,v}\ f(x)+\lambda g(v)\quad\text{s.t.}\quad x=v.$$

这样，$x$专注数据一致性，$v$专注图像先验，等式约束再要求二者最终达成一致。完整的增广拉格朗日配方与 scaled form 推导见 [[1 ADMM变量分裂与scaled form]]；这里直接使用结果：

$$x^{k+1}\in\operatorname*{argmin}_x\left[f(x)+\frac{\rho}{2}\lVert x-(v^k-u^k)\rVert_2^2\right],$$

$$v^{k+1}\in\operatorname*{argmin}_v\left[g(v)+\frac{1}{2\sigma^2}\lVert v-(x^{k+1}+u^k)\rVert_2^2\right],\qquad \sigma=\sqrt{\frac{\lambda}{\rho}},$$

$$u^{k+1}=u^k+x^{k+1}-v^{k+1}.$$

第一步在数据项与“别离当前共识太远”之间折中；第二步恰好是一个 Gaussian 形式的去噪问题；第三步累计$x$与$v$的不一致，推动后续迭代纠正偏差。这个结构解释了 ADMM 为什么适合成为 PnP 的骨架。

## 真正的 plug-and-play 发生在先验步

若$g$是总变差等显式先验，第二步求的是$g$的 proximal mapping。PnP 的动作很小，却改变了理论问题：

$$v^{k+1}=\mathcal D_{\sigma_k}(x^{k+1}+u^k).$$

我们不再先写出$g$，而是把已有图像 denoiser $\mathcal D_{\sigma_k}$直接插入。算法于是成为两个模块的循环：

- `反演模块`：知道成像系统与观测$y$，由$f$保证数据一致性；
- `去噪模块`：知道什么样的局部或非局部结构像自然图像，但不必给出闭式先验；
- `对偶协调`：记住两个模块尚未消除的分歧。

这种替换的收益是模块化：改变前向模型不必重造 denoiser，改变 denoiser也不必重写数据项。更细的变量角色与$\sigma_k$解释见 [[3 PnP替换与denoiser]]。

但不能倒推说“任何 denoiser 都等于某个先验的 proximal mapping”。对一般非线性 denoiser，对应的$g$可能未知，甚至未必存在。因此“仍然在最小化$f+\lambda g$”不是无条件事实。

## 为什么经典 ADMM 收敛定理不能直接照搬

经典 ADMM 的定理针对显式优化问题。若两个目标函数闭、正常、凸，且未增广拉格朗日函数存在鞍点，则原始残差趋于零、目标值趋于最优值、对偶变量趋于对偶最优点。这里的每个结论都能指回一个明确的目标与最优性条件，详见 [[Boyd 等 - 2011 - Distributed Optimization and Statistical Learning via the Alternating Direction Method of Multipliers|Boyd et al. 2011 — ADMM]] 和 [[2 primal-dual residual与停止准则]]。

PnP 替换后，一般 denoiser 未必提供这样的显式目标。若进一步假设 denoiser 非扩张且具有对称梯度，可以尝试恢复 proximal mapping 解释；然而这类条件对自适应非线性 denoiser很难证明，论文甚至给出非局部均值出现扩张的反例。

所以 [[Chan 等 - 2017 - Plug-and-Play ADMM for Image Restoration - Fixed Point Convergence and Applications|Chan et al. 2017 — PnP-ADMM]] 换了一个更弱也更现实的问题：不先证明算法找到某个未知目标的最优解，而先证明整个迭代会进入稳定状态。

## bounded denoiser 与 continuation 如何让状态稳定

论文要求 denoiser 满足

$$\frac{1}{n}\lVert\mathcal D_\sigma(x)-x\rVert_2^2\le\sigma^2C.$$

这称为 bounded denoiser。它没有比较两个输入之间的 Lipschitz 比例，只限制一次去噪不能把输入改变得超过$O(\sigma\sqrt n)$。当$\sigma$变小时，denoiser 渐近接近恒等映射。

接着引入 continuation：令$\rho_k$不减，并保持

$$\sigma_k=\sqrt{\frac{\lambda}{\rho_k}}.$$

单调方案每次把$\rho_k$乘以固定$\gamma>1$，于是$\sigma_k\to0$；自适应方案只在状态变化没有按比例下降时增大$\rho_k$。这给出了两种稳定机制：

- 若反复提高$\rho_k$，数据步与去噪步的改变量被越来越小的上界控制；
- 若保持$\rho_k$，更新规则已经检测到状态变化量按$\eta<1$几何收缩。

因此不能把自适应证明粗略写成“$\rho_k$一定趋于无穷”。真正关键的是，无论落在哪个分支，相邻状态的距离最终都受一个可求和的几何序列控制。

再加上图像域$[0,1]^n$上$f$的梯度一致有界，论文证明$\theta^k=(x^k,v^k,u^k)$是 Cauchy 序列，因而收敛到某个$\theta^\star$。完整条件与证明骨架见 [[4 Bounded denoiser与fixed-point convergence]]。

## 最关键的边界：稳定点不等于优化最优点

现在可以严格区分两句话：

- `经典 ADMM 收敛`：在经典假设下，残差、目标值和对偶变量与一个显式凸优化问题的最优性相联系。
- `本文 PnP fixed-point convergence`：在 bounded denoiser、$f$梯度有界和 continuation 条件下，$(x^k,v^k,u^k)$收敛到稳定状态，且$x^k-v^k\to0$。

第二句话没有自动包含第一句话。fixed point 可以告诉我们迭代停止变化，却不能单凭这一点告诉我们它是哪个显式目标的全局最优解。PnP 的理论价值在于对一般 denoiser 给出可验证的稳定性边界，而不是悄悄恢复一个尚未证明存在的先验。

## 两种停止问题也必须分开

经典 ADMM 监测

$$r^{k+1}=x^{k+1}-v^{k+1},\qquad s^{k+1}=-\rho(v^{k+1}-v^k),$$

分别对应原始可行性与对偶可行性。PnP 论文为 fixed-point 目标监测

$$\Delta_{k+1}=\frac{1}{\sqrt n}\left(\lVert x^{k+1}-x^k\rVert_2+\lVert v^{k+1}-v^k\rVert_2+\lVert u^{k+1}-u^k\rVert_2\right).$$

$\Delta_{k+1}$小表示整个迭代状态已基本不动。虽然$u^{k+1}-u^k=x^{k+1}-v^{k+1}$也包含一致性信息，它仍不是经典 primal/dual residual 判据的改名版本。

## 本章收束

PnP-ADMM 的逻辑链是：

`MAP 的数据项与先验耦合 → 一致性分裂形成反演步与 proximal 去噪步 → 用 denoiser 替换显式 proximal → 显式目标解释可能消失 → 用 bounded 条件与 continuation 控制迭代增量 → 得到 fixed-point convergence。`

掌握这一链条后，使用 PnP 时最重要的不是背下三行迭代，而是每次都问清三件事：数据项是否忠实于观测模型，denoiser 满足什么可验证条件，当前声称的是优化最优性还是仅仅 fixed-point 稳定性。Chan 论文如何把这一框架落到超分辨率和单光子成像，见 [[5 PnP-ADMM的图像恢复实现与案例]]；更现代的深度 denoiser 工程化扩展见 [[6 DPIR与深度denoiser prior]]，注意 DPIR 使用 HQS 而不是 ADMM。
