---
tags:
  - 暑假学习
  - 学习路径
  - FISTA
  - 近端梯度
  - 逆问题
  - 文献总结
created: 2026-09-04
updated: 2026-09-04
title: FISTA原论文文献总结
paper: A Fast Iterative Shrinkage-Thresholding Algorithm for Linear Inverse Problems
authors: Amir Beck; Marc Teboulle
year: 2009
citekey: beckFastIterativeShrinkageThresholding2009
---

# FISTA 原论文文献总结

> [!cite] 原始文献
> [[Beck 等 - 2009 - A Fast Iterative Shrinkage-Thresholding Algorithm for Linear Inverse Problems|Beck & Teboulle, 2009 — FISTA]] · [Zotero](zotero://select/library/items/EJXDP5IE)
>
> Amir Beck, Marc Teboulle. *A Fast Iterative Shrinkage-Thresholding Algorithm for Linear Inverse Problems*. SIAM Journal on Imaging Sciences, 2009, 2(1): 183–202.
>
> DOI: `10.1137/080716542`

## 一、这篇论文真正要解决什么问题

论文的起点不是“怎样设计一个带动量的算法”，而是一个非常具体的矛盾：

- 线性逆问题规模很大，矩阵可能是稠密的；
- 内点法等高精度优化方法在这种规模下代价太高；
- ISTA 每轮只需矩阵向量乘法和一次 shrinkage，足够便宜；
- 但 ISTA 的函数值收敛很慢。

所以作者真正提出的问题是：

> 能不能在几乎不增加 ISTA 单轮计算量的前提下，把它的全局收敛速度明显提高？

论文给出的答案就是 FISTA。

它最重要的贡献不是简单加入一个“惯性项”，而是证明：对一类凸复合问题，FISTA 保持 ISTA 的基本近端步计算量，却把最坏情形函数值收敛率从 $O(1/k)$ 提升到 $O(1/k^2)$。

---

## 二、从线性逆问题为什么会走到 $\ell_1$ 正则化

论文从线性观测模型出发：

$$Ax=b+w.$$

其中：

- $A$：已知前向算子，例如模糊算子；
- $x$：待恢复的真实信号或图像；
- $b$：观测数据；
- $w$：噪声或扰动。

最直接的方法是最小二乘：

$$\min_x\lVert Ax-b\rVert_2^2.$$

但逆问题中 $A$ 常常病态甚至近似奇异，单纯追求数据拟合会把噪声一起放大。因此需要正则化。

论文重点考虑的是

$$\min_x F(x)=\lVert Ax-b\rVert_2^2+\lambda\lVert x\rVert_1.$$

在图像复原里，$x$还可以理解成图像在小波域中的系数。$\ell_1$项鼓励这些系数稀疏，因此能够在满足观测的同时抑制噪声并保留重要结构。

但论文并没有把理论限制死在 $\ell_1$ 上。后面的真正分析对象被提升为一般复合问题：

$$\min_x F(x)=f(x)+g(x),$$

其中：

- $f$：凸、可微，并且梯度 Lipschitz 连续；
- $g$：凸，可以不可微；
- 最优解集合非空。

这一步很重要：FISTA 的理论本质上不是“软阈值算法的技巧”，而是一个`加速近端梯度方法`。

与已有学习笔记的联系见 [[近端梯度为什么自然出现]]。

---

## 三、论文最核心的数学对象不是 FISTA，而是局部二次模型 $Q_L$

论文定义

$$Q_L(x,y)=f(y)+\langle x-y,\nabla f(y)\rangle+\frac{L}{2}\lVert x-y\rVert_2^2+g(x).$$

这里 $y$ 是当前展开点，$x$ 是下一步候选变量。

这个模型可以拆成三部分理解：

- $f(y)+\langle x-y,\nabla f(y)\rangle$：用一阶 Taylor 模型近似光滑项；
- $\frac{L}{2}\lVert x-y\rVert^2$：限制线性化模型不能无限往下跑，同时利用 Lipschitz 条件控制模型误差；
- $g(x)$：不可微正则项不做线性化，直接保留。

这与 [[为什么在线性化模型中加入二次惩罚]] 的核心逻辑完全一致。

定义

$$p_L(y)=\operatorname*{argmin}_x Q_L(x,y).$$

整理平方项后有

$$p_L(y)=\operatorname*{argmin}_x\left\{g(x)+\frac{L}{2}\left\lVert x-\left(y-\frac1L\nabla f(y)\right)\right\rVert_2^2\right\}.$$

用现代近端算子记法就是

$$p_L(y)=\operatorname{prox}_{g/L}\left(y-\frac1L\nabla f(y)\right).$$

因此论文里的 $p_L$ 就是我们现在熟悉的`梯度步 + proximal step`。

---

## 四、为什么 $L$ 必须和 Lipschitz 条件联系起来

若 $\nabla f$ 的 Lipschitz 常数为 $L(f)$，则对任何 $L\ge L(f)$，有

$$f(x)\le f(y)+\langle x-y,\nabla f(y)\rangle+\frac{L}{2}\lVert x-y\rVert_2^2.$$

于是

$$F(x)\le Q_L(x,y).$$

这说明 $Q_L$ 不是任意构造出来的局部模型，而是 $F$ 的一个可控上界。

所以选择 $L$ 本质上是在保证：

> 我们虽然只使用 $f$ 的一阶信息，但加上的二次项足够大，使局部模型不会低估真实函数。

这也是论文中 backtracking 的依据：如果 $L(f)$ 不知道，就不断增大 $L$，直到满足

$$F(p_L(y))\le Q_L(p_L(y),y).$$

因此 backtracking 并不是为了“找一个经验上好用的步长”，而是在自动寻找一个足以支撑局部上界关系的曲率尺度。

对应前置笔记：[[Lipschitz 条件]]。

---

## 五、Lemma 2.3：整篇论文真正的证明发动机

论文后面 ISTA 和 FISTA 的收敛率证明都依赖 Lemma 2.3。

在满足

$$F(p_L(y))\le Q_L(p_L(y),y)$$

时，对任意 $x$ 有

$$F(x)-F(p_L(y))\ge\frac{L}{2}\lVert p_L(y)-y\rVert_2^2+L\langle y-x,p_L(y)-y\rangle.$$

它的意义不是公式本身，而是它把两类量联系起来了：

- 左边：目标函数值差；
- 右边：迭代点之间的欧氏距离和内积。

于是后面就可以把函数值下降问题转化成距离平方的望远镜求和。

这也是理解论文证明最重要的一点：

> 作者不是直接对 FISTA 的三条递推式做神秘代数，而是在不断把目标函数误差改写成能够前后抵消的“能量差”。

---

## 六、ISTA 在论文里扮演什么角色

一般 ISTA 写成

$$x_k=p_{L_k}(x_{k-1}).$$

如果是 $\ell_1$ 正则化，则 $p_L$ 进一步变成软阈值：

$$x_k=\mathcal T_{\lambda/L}\left(x_{k-1}-\frac1L\nabla f(x_{k-1})\right).$$

软阈值本身的推导见 [[软阈值分析]]。

论文首先证明 ISTA 的全局函数值上界：

$$F(x_k)-F(x^*)\le\frac{\alpha L(f)\lVert x_0-x^*\rVert_2^2}{2k}.$$

其中：

- 常步长时 $\alpha=1$；
- backtracking 时 $\alpha=\eta$。

因此 ISTA 是

$$F(x_k)-F^*=O(1/k).$$

这一步不是论文的附属内容，而是必须先建立的比较基线。

如果目标误差要求小于 $\varepsilon$，ISTA 的迭代复杂度大致是

$$k=O(1/\varepsilon).$$

问题于是变得非常明确：有没有办法保持每轮还是一次 $p_L$，却把这个复杂度进一步降低？

---

## 七、FISTA 到底改了 ISTA 的哪一步

ISTA 是直接在上一个点做近端梯度步：

$$x_k=p_L(x_{k-1}).$$

FISTA 不再从 $x_{k-1}$ 出发，而是从一个外推点 $y_k$ 出发：

$$x_k=p_L(y_k).$$

然后引入标量序列

$$t_{k+1}=\frac{1+\sqrt{1+4t_k^2}}{2},\qquad t_1=1,$$

再令

$$y_{k+1}=x_k+\frac{t_k-1}{t_{k+1}}(x_k-x_{k-1}).$$

所以标准 FISTA 的完整结构是：

- 在 $y_k$ 上做一次与 ISTA 完全相同的近端梯度步；
- 用 $x_k-x_{k-1}$ 形成一个外推方向；
- 用 $\frac{t_k-1}{t_{k+1}}$ 控制外推强度。

这意味着 FISTA 的主要单轮成本仍然与 ISTA 相同，额外增加的只是几个向量加减与标量运算。

这正是论文标题里 `Fast` 的关键：不是靠每一步算得更精，而是在几乎相同单步成本下减少达到同等精度所需的迭代次数。

---

## 八、$t_k$ 递推为什么是这个奇怪形式

最容易误解的是：$t_k$ 不是作者凭经验调出的“动量参数”。

它被设计成满足

$$t_k^2=t_{k+1}^2-t_{k+1}.$$

也就是

$$t_{k+1}^2-t_{k+1}-t_k^2=0.$$

解这个二次方程的正根，正好得到

$$t_{k+1}=\frac{1+\sqrt{1+4t_k^2}}2.$$

为什么非要满足这个关系？

因为在 Lemma 4.1 的证明中，作者要把 $v_k=F(x_k)-F(x^*)$ 前面的系数整理成相邻两项的

$$t_k^2v_k-t_{k+1}^2v_{k+1},$$

只有上述递推关系成立，目标函数误差才能和另一个距离型“能量”一起形成望远镜结构。

所以：

> $t_k$ 的递推首先是为证明中的能量递推服务的，动量系数只是它最终在算法上的表现。

---

## 九、FISTA 的证明为什么要构造一个新的“能量”

论文定义

$$v_k=F(x_k)-F(x^*),$$

以及

$$u_k=t_kx_k-(t_k-1)x_{k-1}-x^*.$$

然后证明

$$\frac{2}{L_k}t_k^2v_k-\frac{2}{L_{k+1}}t_{k+1}^2v_{k+1}\ge\lVert u_{k+1}\rVert_2^2-\lVert u_k\rVert_2^2.$$

把它移项，可以看成某种组合能量

$$\frac{2}{L_k}t_k^2v_k+\lVert u_k\rVert_2^2$$

不会增长。

这就是 FISTA 证明的核心结构。

真正被控制的并不是单独的 $F(x_k)$，而是

- 带 $t_k^2$ 权重的函数值误差；
- 一个由当前点、上一点和最优点共同构成的距离项。

二者组合起来形成了可以跨迭代传递的势函数。

因此从现代角度看，FISTA 的证明已经带有非常明显的 `Lyapunov / energy function` 思想。

---

## 十、为什么最终是 $O(1/k^2)$

论文进一步证明

$$t_k\ge\frac{k+1}{2}.$$

另一方面，由能量递推可得

$$\frac{2}{L_k}t_k^2\big(F(x_k)-F(x^*)\big)\le\lVert x_0-x^*\rVert_2^2.$$

于是

$$F(x_k)-F(x^*)\le\frac{L_k\lVert x_0-x^*\rVert_2^2}{2t_k^2}.$$

再利用 $t_k\ge(k+1)/2$，得到论文最核心的 Theorem 4.4：

$$F(x_k)-F(x^*)\le\frac{2\alpha L(f)\lVert x_0-x^*\rVert_2^2}{(k+1)^2}.$$

因此

$$F(x_k)-F^*=O(1/k^2).$$

如果要求目标误差不超过 $\varepsilon$，FISTA 所需迭代次数约为

$$k=O(1/\sqrt\varepsilon).$$

而 ISTA 是

$$k=O(1/\varepsilon).$$

这就是论文意义上的“加速”。

---

## 十一、必须特别区分：论文证明的是函数值收敛率，不是迭代点收敛率

论文的核心结论是

$$F(x_k)-F(x^*)=O(1/k^2).$$

这里控制的是`目标函数值误差`。

它并没有在本文中证明标准 FISTA 的迭代点满足类似

$$\lVert x_k-x^*\rVert=O(1/k^2)$$

这样的结论。

也不能从 $O(1/k^2)$ 自动推断：

- $x_k$ 每一轮都更接近 $x^*$；
- 迭代点单调收敛；
- FISTA 的迭代点在所有问题上一定比 ISTA 更快。

后来的 PnP-FISTA 文献之所以专门研究 `iterate convergence`，就是因为这与原论文的 objective convergence 是不同的问题。

这一点与 [[0 PnP-FISTA（讲解版）]] 中“目标值加速不等于迭代点加速”的区分直接相连。

---

## 十二、标准 FISTA 也不保证目标值逐轮单调下降

ISTA 中论文明确指出 $F(x_k)$ 是非增的。

但标准 FISTA 的 $y_k$ 是由外推产生的，可能越过当前点，因此不能把 ISTA 的单调性直接搬过来。

FISTA 的理论保证是全局上界

$$F(x_k)-F^*\le O(1/k^2),$$

而不是

$$F(x_{k+1})\le F(x_k)\quad\text{对每个 }k\text{都成立}.$$

所以实际画出 FISTA 目标值曲线时出现局部振荡，并不与原论文的 $O(1/k^2)$ 保证矛盾。

---

## 十三、常步长与 backtracking 两个版本有什么关系

### 常步长 FISTA

若已知 $L(f)$，可以直接固定

$$L=L(f),$$

每轮运行

$$x_k=p_L(y_k).$$

### Backtracking FISTA

若 $L(f)$ 不容易计算，则从 $L_{k-1}$ 开始，以比例 $\eta>1$ 增大，直到满足

$$F(p_{\bar L}(y_k))\le Q_{\bar L}(p_{\bar L}(y_k),y_k).$$

然后设置 $L_k=\bar L$。

论文证明 backtracking 版本同样保持 $O(1/k^2)$，只是理论常数中多出 $\alpha=\eta$。

所以 backtracking 没有改变 FISTA 的核心加速结构，只是在每一轮自动寻找一个合法的局部曲率上界。

---

## 十四、论文中的图像去模糊实验说明了什么

作者使用小波稀疏表示的图像去模糊问题比较：

- ISTA；
- FISTA；
- monotone TwIST。

典型实验包括 `cameraman` 和一个简单测试图像。

论文报告的现象非常明确：

- 在相同迭代次数下，FISTA 的目标函数值明显更低；
- FISTA 更快得到视觉质量较好的重建；
- 某些实验中，FISTA 几百次迭代达到的精度，ISTA 或 MTWIST 需要数千甚至更多迭代；
- 在一个最优值已知的最小二乘测试中，10000 次迭代后 FISTA 的误差约达到 $10^{-7}$，而 ISTA 和 MTWIST 仍约在 $10^{-3}$、$10^{-4}$ 量级。

这些实验支持论文的核心主张：

> 加速不仅存在于理论最坏情形界里，在实际图像逆问题上也能表现为显著的迭代次数优势。

但需要注意，这些实验是 2009 年论文中的初步数值验证，不能把“FISTA 总是比任何后续方法快”作为论文结论。

---

## 十五、这篇论文相对于 Nesterov 的真正创新在哪里

作者明确说明，FISTA 的加速思想受到 Nesterov 1983 年光滑凸优化最优一阶方法的启发。

真正的新意在于：

- 原始 Nesterov 加速主要针对光滑凸目标；
- 本文处理 $f+g$，其中 $g$ 可以不可微；
- 每轮仍然只需要一次主要的 projection-like / proximal 操作；
- 仍然能够得到 $O(1/k^2)$ 的函数值复杂度。

所以 FISTA 可以理解成：

`Nesterov 加速思想 + forward-backward / proximal gradient 结构`。

它把加速一阶方法从光滑优化自然推进到了非光滑复合凸优化。

---

## 十六、论文的适用条件

原论文的主要理论框架要求：

- $f$ 为凸函数；
- $f$ 可微；
- $\nabla f$ Lipschitz 连续；
- $g$ 为连续凸函数，可以不可微；
- 问题存在最优解；
- $p_L(y)$ 可以有效计算。

对经典 $\ell_1$ 问题，$p_L$ 就是软阈值，因此非常便宜。

如果 $g$ 的 proximal 子问题本身很难求，FISTA 每轮便宜这一优势就会消失。

此外，原文还指出其分析可扩展到实 Hilbert 空间和带凸约束的情形，但具体计算是否仍然便宜，取决于相应的 proximal / projection 操作是否简单。

---

## 十七、论文没有解决什么

这篇论文非常经典，但阅读时必须把后来的理论与原论文分开。

原论文没有解决：

- 标准 FISTA 迭代点在最一般条件下是否收敛；
- 迭代点的线性收敛率；
- 非凸 $f+g$ 的一般加速理论；
- 任意 learned denoiser 替代 prox 后的收敛；
- 自适应重启、monotone FISTA 等后来变体；
- 强凸情况下更强的速率；
- PnP-FISTA 中 denoiser 是否对应显式 regularizer。

因此以后看到这些结论时，应当把它们视为 FISTA 之后的扩展，而不是 Beck–Teboulle 2009 已经证明的内容。

---

## 十八、从论文角度重新理解 FISTA 的三行公式

如果只记算法，FISTA 很容易变成三条孤立公式：

$$x_k=p_L(y_k),$$

$$t_{k+1}=\frac{1+\sqrt{1+4t_k^2}}2,$$

$$y_{k+1}=x_k+\frac{t_k-1}{t_{k+1}}(x_k-x_{k-1}).$$

但按照原论文真正的逻辑，它们应该这样理解：

- `第一行`：仍然是 ISTA 的近端梯度步，负责处理 $f+g$；
- `第二行`：专门构造能产生 $t_k^2$ 权重的递推关系；
- `第三行`：把前两次迭代信息组合成下一轮展开点，使 Lemma 2.3 可以形成跨迭代的能量递推。

所以 FISTA 的加速并不是“因为用了历史方向所以更快”这么简单。

更准确地说：

> 作者设计了一套特殊的外推结构，使目标误差乘上一个约为 $k^2$ 增长的权重后仍能被统一能量上界控制，因此目标函数误差自然衰减为 $O(1/k^2)$。

---

## 十九、整篇论文的论证链

可以把 Beck–Teboulle 2009 压缩成下面一条链：

`病态线性逆问题 → 需要正则化 → ℓ1 带来非光滑复合目标 → ISTA 用局部二次上界做近端梯度更新 → Lemma 2.3 建立函数值与距离的关键不等式 → 证明 ISTA 只有 O(1/k) → 问能否保持相同近端步成本但加速 → 把近端步从 x_{k-1} 改在特殊外推点 y_k 上 → 构造 t_k 使加权函数误差形成望远镜能量递推 → t_k≈k/2 → 得到 FISTA 的 O(1/k²) 函数值收敛率 → 图像去模糊实验验证明显加速。`

---

## 二十、读完论文后最应该留下的几个结论

- FISTA 的研究对象本质上是凸复合优化 $f+g$，不是只针对 $\ell_1$。
- $p_L(y)$ 就是今天所说的 proximal-gradient step。
- 二次项来自 $\nabla f$ 的 Lipschitz 上界，不是随意添加的稳定项。
- Lemma 2.3 是 ISTA 与 FISTA 两个速率证明共同的基础。
- ISTA 的理论函数值速率是 $O(1/k)$。
- FISTA 每轮主要计算量基本与 ISTA 相同，但函数值速率提升到 $O(1/k^2)$。
- $t_k$ 的递推是为了形成证明所需的能量关系，而不只是经验动量。
- 原论文中的“加速”是`目标函数值最坏情形复杂度加速`。
- $O(1/k^2)$ 不等于标准 FISTA 的迭代点距离也按 $O(1/k^2)$ 收敛。
- 标准 FISTA 不要求目标函数值逐轮单调。

这篇论文真正重要的地方，是把“近端梯度每轮很便宜”和“Nesterov 式 $O(1/k^2)$ 加速”非常干净地结合在了一起，从而奠定了今天大量复合优化和逆问题一阶算法的基本模板。
