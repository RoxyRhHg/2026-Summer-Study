---
tags:
  - 暑假学习
  - 学习路径
  - PnP-FISTA
  - FISTA
  - Plug-and-Play
  - 逆问题
created: 2026-09-03
updated: 2026-09-04
---

# PnP-FISTA：从 prox 替换到 scaled 几何与迭代收敛

> [!cite] Zotero 资料基线
> - [[Gavaskar 等 - 2021 - On Plug-and-Play Regularization using Linear Denoisers|Gavaskar 等 2021：线性 denoiser 与 scaled proximal map]] · [Zotero](zotero://select/library/items/R3JAW5G3)
> - [[Sinha 等 - 2024 - FISTA Iterates Converge Linearly for Denoiser-Driven Regularization|Sinha 与 Chaudhury 2024：PnP-FISTA 迭代收敛]] · [Zotero](zotero://select/library/items/IJAI4UMU)
> - FISTA 前置：[[1 FISTA（讲解版）]]；PnP 前置：[[0 PnP-ADMM（讲解版）]]。
> - 本篇只按上述 Zotero MinerU 原文梳理“为什么从 PnP-ADMM 继续走到 PnP-FISTA、为什么 scaled geometry 必须出现、为什么 objective convergence 与 iterate convergence 必须分开”。

## 先抓住新的问题：PnP-ADMM 已经能插 denoiser，为什么还要换成 FISTA 骨架

在 [[0 PnP-ADMM（讲解版）]] 中，我们已经解决了最关键的 PnP 接口问题：ADMM 通过变量分裂把先验步骤单独暴露出来，然后把

$$\operatorname{prox}_g\quad\longrightarrow\quad\mathcal D$$

替换成现成 denoiser。

那为什么还要继续学 PnP-FISTA？

因为 PnP-ADMM 每轮的数据步通常是一个完整子问题：

$$x^{k+1}=\operatorname*{argmin}_x\left\{f(x)+\frac{\rho}{2}\lVert x-c_k\rVert^2\right\}.$$

如果$f$对应复杂 forward model，这一步可能需要解线性方程、做内层迭代，或者专门利用$A$的结构。

而如果$f$本来就是光滑的，例如

$$f(x)=\frac12\lVert Ax-b\rVert_2^2,$$

那么 [[1 FISTA（讲解版）|FISTA]] 已经给我们一个更轻的求解框架：

$$\text{一次梯度}\quad+\quad\text{一次 prox}\quad+\quad\text{一次历史外推}.$$

于是新的问题自然出现：

> 能不能像 PnP-ADMM 一样，把 FISTA 里的 prox 接口也交给 denoiser，同时保留它简单的一阶数据步和动量结构？

形式上答案非常简单。但真正的难点恰恰从这一步之后才开始：`算法能跑`不等于`它仍然拥有 FISTA 的理论含义`。

---

## 一、FISTA 里真正可替换的接口是哪一行

经典 FISTA 解决

$$\min_x F(x)=f(x)+g(x),$$

其中$f$光滑，$g$可以不可微但 prox 易算。

把一轮更新按职责拆开，可以写成三步。

先根据历史点构造外推点：

$$y_k=x_{k-1}+\alpha_{k-1}(x_{k-1}-x_{k-2}).$$

然后只根据数据项$f$做梯度步：

$$z_k=y_k-\gamma\nabla f(y_k).$$

最后把结构项$g$放回来：

$$x_k=\operatorname{prox}_{\gamma g}(z_k).$$

所以 FISTA 本来就有一个非常清晰的模块边界：

- $\nabla f$：知道 forward model 和观测；
- $\operatorname{prox}_{\gamma g}$：知道结构先验；
- $\alpha_k$：只利用历史迭代点决定下一轮从哪里出发。

PnP-FISTA 做的核心替换就是

$$\operatorname{prox}_{\gamma g}(z_k)\quad\longrightarrow\quad\mathcal D(z_k),$$

从而得到

$$x_k=\mathcal D\left(y_k-\gamma\nabla f(y_k)\right).$$

注意这里没有把整个 FISTA 都“神经网络化”。被换掉的只是`先验接口`：

$$\text{gradient：数据一致性}\longrightarrow\text{denoiser：图像先验}\longrightarrow\text{momentum：利用历史信息}.$$

这就是 PnP-FISTA 最直观的工程结构。

但如果只停在这里，就会产生一个危险错觉：

> 既然只是把 prox 换成一个效果更好的去噪器，FISTA 原来的$O(1/k^2)$应该也还在吧？

不能这样推。

---

## 二、为什么“能运行”不等于“仍然是 FISTA”：原来的证明到底依赖了什么

经典 FISTA 的$O(1/k^2)$不是由“有梯度 + 有动量 + 有一个图像修正模块”这三个外观共同保证的。

它真正依赖的是：

$$x_k=\operatorname{prox}_{\gamma g}\left(y_k-\gamma\nabla f(y_k)\right)$$

确实是同一个显式目标

$$f(x)+g(x)$$

上的 accelerated proximal-gradient step。

如果把 prox 换成任意 denoiser

$$x_k=\mathcal D\left(y_k-\gamma\nabla f(y_k)\right),$$

那么第一件必须重新问的事情是：

> 这个$\mathcal D$到底是不是某个$g$的 proximal map？

如果答案是“是”，我们就可能把 PnP-FISTA 重新解释成普通 FISTA，只不过这个$g$没有一开始显式写出来。

如果答案是“不知道”，那原来的目标函数解释就断了；既然连$F(x)=f(x)+g(x)$是什么都不一定知道，就不能直接搬用

$$F(x_k)-F^\star=O(1/k^2).$$

所以 PnP-FISTA 的第一条理论路线并不是马上研究动量，而是先回头追问 denoiser 本身：

> 它有没有可能其实就是一个 prox？

---

## 三、对称线性 denoiser 还能像普通 prox，真正麻烦的是非对称 kernel denoiser

对于一类对称、半正定、谱合适的线性 denoiser

$$\mathcal D(x)=Wx,$$

已有结果可以把$W$解释成某个凸函数的标准 Euclidean proximal map。

这意味着在这些特殊情形下，PnP 的“黑盒先验”并没有真的彻底脱离凸优化：它只是把一个隐式正则项藏进了矩阵$W$里。

问题出在实际很常见的一类 kernel denoiser。

例如 NLM 型线性化后常写成

$$W=D^{-1}K,$$

其中$K$是对称 kernel matrix，$D$是行和归一化矩阵。

因为左乘了$D^{-1}$，$W$通常不再关于标准欧氏内积对称。

如果我们把“prox”只理解成

$$\operatorname{prox}_g(y)=\operatorname*{argmin}_x\left\{\frac12\lVert x-y\rVert_2^2+g(x)\right\},$$

就很容易得出一个过早结论：

> $W$不对称，所以它肯定不是某个凸$g$的 prox。

Gavaskar 2021 的关键推进正是指出：问题可能不在$W$本身，而在我们默认使用了错误的`距离几何`。

---

## 四、为什么“非对称”会逼出 scaled geometry：prox 本来就依赖你用什么距离衡量“靠近”

回忆 prox 的本质：

$$\operatorname{prox}_g(y)=\operatorname*{argmin}_x\left\{g(x)+\frac12\lVert x-y\rVert_2^2\right\}.$$

它一半在看$g(x)$，另一半在看$x$离输入$y$有多远。

但“多远”并不一定非要由普通欧氏距离

$$\lVert x-y\rVert_2^2$$

衡量。

给定$H\succ0$，定义加权内积和加权范数

$$\langle x,y\rangle_H=x^THy,\qquad\lVert x\rVert_H=\sqrt{x^THx}.$$

于是可以定义

$$\operatorname{prox}_{g,H}(y)=\operatorname*{argmin}_x\left\{g(x)+\frac12\lVert x-y\rVert_H^2\right\}.$$

这就是 `scaled proximal map`。

直觉上，$H$改变的是空间里不同方向的长度尺度：在普通欧氏几何中看起来“不对称”的算子，换到合适的加权内积以后，可能恰好就是一个合法的 proximal map。

Gavaskar 2021 证明了一类线性 denoiser 的 scaled-proximal characterization。对 kernel denoiser

$$W=D^{-1}K,$$

可以取

$$H=D,$$

从而把$W$解释成某个闭、正常、凸函数的$D$-scaled proximal map。

这一步解决了第一个理论缺口：

> 非对称 kernel denoiser 并不一定没有隐式凸正则项；它可能只是属于另一个 Hilbert-space geometry。

但新的问题马上又出现：

> 既然 prox 已经换到$H$几何里了，能不能仍然保留原来的 Euclidean gradient step，只把最后一行换成 scaled prox？

答案是不行。因为这样会让同一轮里的两个步骤生活在不同的几何里。

---

## 五、为什么 prox 的几何一换，梯度也必须一起换：从局部模型推一次就知道了

这是 PnP-FISTA 最容易被一句“加一个$H^{-1}$”掩盖掉的地方。

在普通 Euclidean FISTA 里，近端梯度来自局部模型

$$f(y)+\langle\nabla f(y),x-y\rangle+\frac{\rho}{2}\lVert x-y\rVert_2^2+g(x).$$

最小化时，把一阶项和二次项配方，中心会变成

$$y-\rho^{-1}\nabla f(y).$$

现在既然 prox 使用的是$H$距离，那么一致的局部模型就应该改成

$$f(y)+\langle\nabla f(y),x-y\rangle+\frac{\rho}{2}\lVert x-y\rVert_H^2+g(x).$$

而

$$\lVert x-y\rVert_H^2=(x-y)^TH(x-y).$$

对$x$求一阶条件时，二次项贡献的是

$$\rho H(x-y).$$

因此局部二次部分的中心不再是

$$y-\rho^{-1}\nabla f(y),$$

而是

$$y-\rho^{-1}H^{-1}\nabla f(y).$$

所以 scaled gradient step 并不是“为了预条件而顺手加一个$H^{-1}$”，它是同一个$H$-几何局部模型自己推出来的。

这就得到 scaled PnP-FISTA：

$$x_{k+1}=\mathcal D\left(y_{k+1}-\rho^{-1}H^{-1}\nabla f(y_{k+1})\right).$$

动量更新仍然可以保持原来的 FISTA 形式：

$$t_{k+1}=\frac{1+\sqrt{1+4t_k^2}}2,$$

$$y_{k+1}=x_k+\frac{t_k-1}{t_{k+1}}(x_k-x_{k-1}).$$

现在整个步骤终于重新处在同一个$H$几何里：

$$\text{$H$-gradient step}\longrightarrow\text{$H$-scaled prox / denoiser}\longrightarrow\text{FISTA momentum}.$$

这就是为什么 Gavaskar 2021 不只是“证明 denoiser 是 scaled prox”，而是还必须修改 PnP 算法本身。

---

## 六、这样修完以后，FISTA 的哪部分理论真的能接回来

如果$W$满足论文中的 scaled-prox 条件，并且$f$在$H$几何下是凸且$\rho$-smooth，即

$$\left\lVert H^{-1}\nabla f(x)-H^{-1}\nabla f(y)\right\rVert_H\le\rho\lVert x-y\rVert_H,$$

那么 scaled PnP-FISTA 可以重新解释成一个真正的 accelerated proximal-gradient method。

此时存在相应凸正则项$g$，并能恢复经典类型的目标值结论：

$$f(x_k)+g(x_k)-p^\star=O(1/k^2).$$

这件事非常重要，因为它说明：

> PnP 并不必然意味着“从此没有目标函数”。对一类线性 denoiser，只要找对几何，显式凸优化解释可以被重新接回来。

但这里必须立刻限制结论边界。

Gavaskar 2021 这里恢复的首先是`objective convergence`：

$$F(x_k)-F^\star=O(1/k^2).$$

它不是在说

$$\lVert x_k-x^\star\rVert=O(1/k^2).$$

而且如果改用 Chambolle–Dossal 型动量

$$t_{k+1}=1+\frac{k}{a},\qquad a>2,$$

可以进一步得到迭代点$x_k$本身收敛到 minimizer；但这个结论依然没有自动给出“迭代点线性收敛”的速率。

所以理论问题并没有结束，只是从“有没有目标函数”推进到了下一层：

> 即使目标值已经有$O(1/k^2)$，图像序列$x_k$本身到底怎么收敛？

这正是 Sinha 2024 继续研究的问题。

---

## 七、为什么目标值$O(1/k^2)$还不够：我们实际看到的是$x_k$，不是一串抽象函数值

经典 FISTA 里最容易混淆的两句话是：

$$F(x_k)-F^\star=O(1/k^2)$$

和

$$x_k\to x^\star.$$

它们不是同一个命题。

目标值很接近最优值，不代表迭代点在每个方向上都已经离某个$x^\star$很近；尤其当最优解不唯一，或者问题缺少强凸性时，这种区别会更明显。

而在图像重建里，我们真正保存、显示和继续处理的是$x_k$。

所以新的问题应该写得更直接：

> PnP-FISTA 产生的图像序列本身会不会收敛？如果会，它离极限图像的距离能不能按几何级数下降？

Sinha 与 Chaudhury 2024 为了回答这个更强的问题，没有继续试图覆盖任意黑盒 CNN denoiser，而是主动缩小问题范围。

这不是退步，而是为了换取可验证的强结论。

---

## 八、为什么要把问题限制成“二次数据项 + 固定线性 denoiser”：因为这时整个 PnP 步会变成 affine map

考虑线性逆问题

$$b=A\xi+e,$$

并取二次数据项

$$f(x)=\frac12\lVert Ax-b\rVert_2^2.$$

于是

$$\nabla f(x)=A^T(Ax-b).$$

再假设收敛分析阶段 denoiser 是固定线性算子

$$\mathcal D(x)=Wx.$$

那么一次 PnP 数据步 + 去噪步就是

$$x_k=W\left(y_k-\gamma A^T(Ay_k-b)\right).$$

展开：

$$x_k=W(I-\gamma A^TA)y_k+\gamma WA^Tb.$$

记

$$P_\gamma=W(I-\gamma A^TA),\qquad q=\gamma WA^Tb,$$

就得到非常简单的 affine map：

$$x_k=P_\gamma y_k+q.$$

这一步是整个谱分析能够成立的关键。

普通非线性 denoiser 会让这个映射随输入复杂变化；固定线性$W$则把“forward model + denoiser”的联合行为压缩成一个矩阵$P_\gamma$。

而 FISTA 还有动量：

$$y_k=x_{k-1}+\alpha_{k-1}(x_{k-1}-x_{k-2}).$$

所以$x_k$不仅依赖$x_{k-1}$，还依赖$x_{k-2}$。

这意味着我们不能只把$x_k$当作一阶状态；自然的下一步是把连续两轮合成一个更大的状态。

---

## 九、为什么要升到$2n$维：动量把一阶迭代变成了二阶动力系统

定义

$$z_k=\begin{bmatrix}x_k\\x_{k-1}\end{bmatrix}\in\mathbb R^{2n}.$$

因为$x_k$由$x_{k-1}$和$x_{k-2}$共同决定，所以在这个升维状态下，原来的“二阶递推”重新变成一阶状态转移：

$$z_k=R_kz_{k-1}+s.$$

这里$R_k$包含

- 数据算子$A$；
- denoiser$W$；
- 步长$\gamma$；
- 当前动量系数$\alpha_k$。

对标准 FISTA，$\alpha_k\to1$，因此$R_k$也会趋向一个极限矩阵

$$R_k\to R_\infty.$$

于是“迭代点会不会收敛”的问题，被转化成了一个线性动力系统问题：

> 极限状态转移矩阵$R_\infty$会不会把误差持续缩小？

这时谱半径自然成为核心。

但在真正看$R_\infty$之前，还要先保证基础算子$P_\gamma$没有某个永远无法被压缩的方向。

这就是 Sinha 2024 的两个关键假设为什么会出现。

---

## 十、为什么$\ker(A)\cap\operatorname{fix}(W)=\{0\}$不是一个抽象条件，而是在排除“谁都管不了”的方向

论文对对称线性$W$采用的核心条件包括

$$\sigma(W)\subset[0,1]$$

以及

$$\ker(A)\cap\operatorname{fix}(W)=\{0\}.$$

第一条说$W$的特征值位于$[0,1]$：denoiser 不会在它自己的特征方向上制造大于1的线性放大。

第二条更值得理解。

如果某个非零方向$d$满足

$$d\in\ker(A),$$

那么

$$Ad=0.$$

也就是说，forward model 完全看不见这个方向。沿$d$移动不会改变观测，数据一致性步骤无法纠正它。

如果同一个$d$又满足

$$d\in\operatorname{fix}(W),$$

那么

$$Wd=d.$$

denoiser 也完全不改变这个方向。

于是这个方向就变成了：

> 数据模块看不见，denoiser 也不管。

如果存在这样的非零$d$，整个 PnP 系统就缺少让它收缩的机制。

因此

$$\ker(A)\cap\operatorname{fix}(W)=\{0\}$$

真正排除的是`共同盲区`：任何非零方向都不能同时逃过数据项和 denoiser 的约束。

这让$A$和$W$第一次不再被看成两个独立模块，而是被放到同一个“联合可收缩性”问题里。

---

## 十一、为什么谱半径$<1$就会直接带来迭代点的线性收敛

对于固定矩阵迭代

$$e_{k+1}=Re_k,$$

如果

$$\rho(R)<1,$$

那么所有特征模态最终都会被反复乘上绝对值小于1的因子，误差呈几何级数衰减。

PnP-FISTA 的状态矩阵$R_k$虽然随$k$变化，但当

$$\alpha_k\to1$$

时有

$$R_k\to R_\infty.$$

Sinha 2024 证明在其假设和步长条件下，极限矩阵满足

$$\rho(R_\infty)<1.$$

于是可以得到存在$0\le\beta<1$和常数$C$，使从某个迭代开始

$$\lVert x_k-x^\star\rVert\le C\beta^k.$$

这才叫`迭代点的线性收敛`。

这里的“linear”不要理解成误差按$1/k$下降，而是指对数尺度上呈直线，也就是几何级数

$$\beta^k.$$

对 PnP-FISTA，论文给出的步长范围是

$$0<\gamma<\frac1{\lambda_{\max}(A^TA)}.$$

现在终于可以看清这篇论文为什么选择谱分析：

$$\text{二次$f$ + 线性$W$}\rightarrow\text{一次 PnP 变成 affine map}\rightarrow\text{动量使其成为二阶递推}\rightarrow\text{升到$2n$维状态}\rightarrow\text{检查极限矩阵谱半径}\rightarrow\text{迭代点几何收缩}.$$

---

## 十二、非对称 kernel denoiser 为什么又必须回到$D$-scaled geometry

到这里会发现 Gavaskar 2021 的 scaled geometry 并不是只为“恢复一个显式正则项”服务。

Sinha 2024 在处理非对称 NLM 型 kernel denoiser 时，同样需要回到加权内积

$$\langle x,y\rangle_D=x^TDy.$$

对应的数据梯度变成

$$\nabla_D f(x)=D^{-1}A^T(Ax-b).$$

于是 scaled PnP-FISTA 使用

$$x_k=W\left(y_k-\gamma D^{-1}A^T(Ay_k-b)\right).$$

步长条件相应变成

$$0<\gamma<\frac1{\lambda_{\max}(D^{-1/2}A^TAD^{-1/2})}.$$

这和前面的因果链完全一致：

> 如果$W$只有在$D$几何中才是正确的 proximal / contraction 对象，那么数据梯度、平滑性和谱分析也都必须放到同一个$D$几何里。

所以 `scaled` 不是附加技巧，而是一条贯穿非对称线性 denoiser 理论的结构原则：

$$\text{denoiser 属于哪个几何}\Longrightarrow\text{整套 PnP 迭代就必须在哪个几何里自洽}.$$

---

## 十三、为什么“FISTA 加速”绝不能直接翻译成“$x_k$一定比 ISTA 更快”

现在可以把两种“快”严格分开。

经典 FISTA 的加速结论首先是

$$F(x_k)-F^\star=O(1/k^2),$$

它比 ISTA 的

$$O(1/k)$$

更快。

但 Sinha 2024 关心的是

$$\lVert x_k-x^\star\rVert,$$

也就是迭代点本身离极限图像多远。

这两个量没有一般的一一对应关系。

甚至在论文实验中，不同动量序列下，PnP-FISTA 的迭代点收敛可能比无动量的 PnP-ISTA 更慢。

这并不推翻 FISTA 的经典加速理论，因为经典理论根本没有承诺：

> 对任何具体问题，带 Nesterov 动量的$x_k$都一定以更小的几何因子靠近$x^\star$。

所以以后看到“accelerated”必须继续问三个问题：

- 加速的是`目标值`还是`迭代点`？
- 说的是`最坏情形复杂度`还是当前线性系统的`谱收缩率`？
- 当前 denoiser 是否真的满足使父算法理论成立的条件？

把这三个问题分开，PnP-FISTA 的很多表面矛盾就会消失。

---

## 十四、把 PnP-FISTA 和 PnP-ADMM 放在一起：它们不是“谁更高级”，而是在不同地方付成本

现在可以做一个真正有意义的对照。

| 问题 | PnP-ADMM | PnP-FISTA |
| --- | --- | --- |
| 数据模块 | 解一个带二次罚项的$x$子问题 | 对光滑$f$做一次梯度步 |
| 先验模块 | denoiser | denoiser |
| 额外状态 | $v,u$ | $x_{k-1}$与动量状态 |
| 对$f$的要求 | 可处理更一般、甚至非光滑的凸$f$ | 通常要求$f$光滑并有合适 Lipschitz / scaled smoothness |
| 单轮主要代价 | 数据子问题可能较重 | 通常是一梯度 + 一 denoiser |
| 典型理论入口 | Chan 2017：bounded + continuation → fixed point | Gavaskar 2021：scaled prox → objective；Sinha 2024：线性情形 → iterate linear convergence |
| 最容易误读 | fixed point 当成显式目标最优点 | $O(1/k^2)$目标值当成迭代点一定更快 |

所以阶段 1 学 PnP-ADMM 的目的不是为了得出“ADMM 更稳”，阶段 2 学 PnP-FISTA 也不是为了得出“FISTA 更快”。

真正需要建立的是选择意识：

- 如果$f$的 proximal / 子问题容易解，ADMM 的分裂结构很自然；
- 如果$f$光滑且梯度便宜，FISTA 型 PnP 每轮可能更轻；
- 如果要声称某种收敛性，必须重新检查当前 denoiser、几何、forward model 和所讨论的收敛对象。

---

## 十五、把整条问题链收回来：PnP-FISTA 为什么会自然出现

整个过程可以从 PnP-ADMM 之后继续往下走。

第一个问题是：`PnP-ADMM 已经能插 denoiser，但数据子问题可能较重；如果$f$光滑，能不能换成更轻的一阶骨架？`

于是从 FISTA 出发，把

$$\operatorname{prox}_{\gamma g}\longrightarrow\mathcal D$$

得到 PnP-FISTA。

第二个问题是：`能跑以后，经典 FISTA 的$O(1/k^2)$还能不能直接继承？`

不能。必须先知道 denoiser 是否仍然是某个 prox。

第三个问题是：`非对称 kernel denoiser 不是普通 Euclidean prox，是否就彻底没有优化解释？`

Gavaskar 2021 的回答是：不一定。它可能是$H$-scaled proximal map。

第四个问题是：`既然 prox 换了几何，为什么算法不能只换 denoiser？`

因为同一个局部模型要求梯度也同步变成

$$H^{-1}\nabla f,$$

否则 gradient 与 prox 不在同一个几何里，父算法理论无法重新接上。

第五个问题是：`scaled PnP-FISTA 恢复了$O(1/k^2)$以后，是不是已经知道$x_k$本身收敛多快？`

仍然没有。目标值误差与迭代点距离必须分开。

第六个问题于是变成：`在什么可验证条件下，图像迭代本身能线性收敛？`

Sinha 2024 把范围限制在线性逆问题和固定线性 denoiser，使一次 PnP 变成 affine map，再利用

$$\ker(A)\cap\operatorname{fix}(W)=\{0\}$$

排除共同盲区，并把带动量的二阶递推升到$2n$维状态空间，通过

$$\rho(R_\infty)<1$$

得到

$$\lVert x_k-x^\star\rVert\le C\beta^k.$$

所以 PnP-FISTA 最值得记住的不是一句“FISTA 的 prox 换成 denoiser”，而是完整的理论升级链：

$$\text{FISTA prox 接口}\rightarrow\text{denoiser 替换}\rightarrow\text{目标解释断裂}\rightarrow\text{scaled prox 找回正确几何}\rightarrow\text{梯度也同步 scaled}\rightarrow\text{恢复 objective }O(1/k^2)\rightarrow\text{再区分 iterate convergence}\rightarrow\text{线性系统 + 谱半径证明几何收敛}.$$

> `PnP-FISTA 真正训练的问题意识，是每当我们把一个理论算法里的数学算子换成黑盒模块时，都必须重新问：它还属于原来的数学类别吗？如果只在另一种几何里成立，整套算法是否也随之改变？最后声称的“收敛”究竟是目标值、迭代点，还是仅仅某个 fixed point？`
