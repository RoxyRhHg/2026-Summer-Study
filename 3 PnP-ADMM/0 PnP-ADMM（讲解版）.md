---
tags:
  - 暑假学习
  - 学习路径
  - PnP-ADMM
  - Plug-and-Play
  - 图像恢复
  - 逆问题
created: 2026-09-03
updated: 2026-09-04
---

# PnP-ADMM：从显式先验到可插拔 denoiser

> [!cite] Zotero 资料基线
> - [[Boyd 等 - 2011 - Distributed Optimization and Statistical Learning via the Alternating Direction Method of Multipliers|Boyd 等 2011：ADMM]] · [Zotero](zotero://select/library/items/B7AYLH4I)
> - [[Chan 等 - 2017 - Plug-and-Play ADMM for Image Restoration - Fixed Point Convergence and Applications|Chan 等 2017：PnP-ADMM]] · [Zotero](zotero://select/library/items/2HLFYW77)

> [!important] 先确认这一章为什么会出现在 FISTA 后面
> 这里不是把 [[1 FISTA（讲解版）|FISTA]] 继续升级成 ADMM。FISTA 与 ADMM 是从同一个 $f+g$ 复合问题出发的两种不同分裂框架：前者采用“梯度 + prox”，后者采用“变量分裂 + 共识协调”。如果对“为什么学完 FISTA 突然换成 ADMM”“ADMM 三步从哪里长出来”还不自然，先读 [[为什么从FISTA走向ADMM]]，再从本篇继续。

## 先抓住真正的矛盾：我们不是缺一个优化器，而是“先验接口”太死

图像逆问题常写成

$$\widehat x\in\operatorname*{argmin}_x f(x)+\lambda g(x).$$

这里两部分的职责很清楚：

- $f(x)$：`数据一致性`，它知道观测$y$和前向模型$A$；
- $g(x)$：`图像先验`，它规定什么样的$x$更像合理图像；
- $\lambda$：控制数据与先验之间的权衡。

如果使用总变差、$\ell_1$等显式正则项，这套写法没有问题。但图像处理中很快会遇到一个新的现实：BM3D、NLM、CNN denoiser 等往往比一个简单的显式$g$更会“修图”，可它们首先给我们的不是一个函数$g(x)$，而是一个映射

$$\mathcal D_\sigma:v\longmapsto\text{denoised image}.$$

这时不能直接把目标写成

$$f(x)+\mathcal D_\sigma(x),$$

因为$\mathcal D_\sigma$返回的是图像，不是一个可以和$f(x)$相加的标量代价。

所以真正的问题不是：

> 还能不能再设计一个更复杂的显式正则项？

而是：

> 有没有一种优化算法，能够把“数据一致性”和“先验处理”拆成两个独立接口，使第二个接口本来就长得像“输入一张图，输出一张处理后的图”？

这正是 ADMM 对 PnP 最关键的意义。PnP-ADMM 并不是先有“plug-and-play”这个口号，再随便挑一个 ADMM；而是 ADMM 的变量分裂恰好把先验步骤暴露成了一个可以被 denoiser 接管的模块。

---

## 一、为什么先把一个$x$拆成$x$和$v$：不是多造一个变量，而是把两种职责分开

原问题是

$$\min_x f(x)+\lambda g(x).$$

数据项$f$和先验项$g$都作用在同一个$x$上。数学上这没有错，但从算法接口看，它们被绑在了一起。

ADMM 先复制一个变量$v$：

$$\min_{x,v} f(x)+\lambda g(v)\qquad\text{s.t.}\qquad x=v.$$

乍看之下只是把一个变量变成两个变量，还多加了一条约束，好像把问题弄复杂了。真正的目的却恰好相反：

- $x$只交给$f$，负责解释观测；
- $v$只交给$g$，负责满足图像先验；
- $x=v$保证最后不能各做各的，二者必须达成一致。

因此变量分裂不是为了得到“两张最终图像”，而是为了让同一个最终图像在迭代过程中暂时拥有两个不同职责的副本。

对这个一致性约束构造增广拉格朗日，并使用 scaled dual variable $u$，ADMM 可以写成

$$x^{k+1}=\operatorname*{argmin}_x\left\{f(x)+\frac{\rho}{2}\left\lVert x-(v^k-u^k)\right\rVert_2^2\right\},$$

$$v^{k+1}=\operatorname*{argmin}_v\left\{\lambda g(v)+\frac{\rho}{2}\left\lVert v-(x^{k+1}+u^k)\right\rVert_2^2\right\},$$

$$u^{k+1}=u^k+x^{k+1}-v^{k+1}.$$

这三步先不要背成公式，而要看它们分别在解决什么冲突。

$x$步面对的是：

> 我当然想把$f(x)$降下来，但也不能完全不管$v^k-u^k$这个当前共识。

$v$步面对的是：

> 我当然想让图像更符合先验$g$，但也不能把刚才$x$步的数据修正全部推翻。

而$u$步记录的是

$$x^{k+1}-v^{k+1},$$

也就是两个模块还没有消除的分歧。因为

$$u^k=u^0+\sum_{j=1}^k(x^j-v^j),$$

所以$u$可以直观理解成一份“历史分歧账本”：如果$x$和$v$长期朝不同方向拉扯，后续子问题就会被这份累计偏差持续修正。

到这里，ADMM 已经完成了一件对 PnP 至关重要的事：`数据处理和先验处理不再混在同一个子问题里。`

---

## 二、为什么$v$步突然像“去噪”：prox 接口就是在这里被暴露出来的

现在专门看$v$步：

$$v^{k+1}=\operatorname*{argmin}_v\left\{\lambda g(v)+\frac{\rho}{2}\left\lVert v-(x^{k+1}+u^k)\right\rVert_2^2\right\}.$$

把目标整体除以$\lambda>0$不会改变最小点，于是

$$v^{k+1}=\operatorname*{argmin}_v\left\{g(v)+\frac{\rho}{2\lambda}\left\lVert v-(x^{k+1}+u^k)\right\rVert_2^2\right\}.$$

定义

$$\sigma=\sqrt{\frac{\lambda}{\rho}},$$

就有

$$v^{k+1}=\operatorname*{argmin}_v\left\{g(v)+\frac{1}{2\sigma^2}\left\lVert v-(x^{k+1}+u^k)\right\rVert_2^2\right\}.$$

这一步和前面学习[[近端梯度为什么自然出现]]时看到的 prox 结构本质相同：

$$\text{结构代价}\quad+\quad\text{不要离输入太远的二次代价}.$$

如果把

$$r^k=x^{k+1}+u^k$$

暂时看成一张“待处理图像”，那么$v$步就在回答：

> 在不离$r^k$太远的前提下，怎样找到一张更符合先验$g$的图像？

这已经非常像一个去噪问题。

例如$g$取 TV 时，它就是 TV denoising；如果$g$有易算的 prox，它就是一次 proximal denoising。

这正是 PnP 的入口。我们终于得到一个接口，它的计算语义不是“返回一个代价值”，而是：

$$\text{输入一张图像}\longrightarrow\text{输出一张更符合图像先验的图像}.$$

而这恰好就是现成 denoiser 擅长做的事情。

> [!important] $\sigma$不要机械理解成真实噪声标准差
> 在 Chan 2017 中，$\sigma_k$来自$\sqrt{\lambda/\rho_k}$，形式上对应 Gaussian denoising 的噪声参数，但迭代中的$x^{k+1}+u^k$并不真的等于“干净图像 + IID Gaussian 噪声”。因此这里更稳妥的理解是：$\sigma_k$是控制 denoiser 强度的旋钮。

---

## 三、PnP 真正替换了什么：只换先验模块，但不能假装什么都没变

既然$v$步已经长成一个去噪接口，最自然的动作就是把它直接交给一个 denoiser：

$$v^{k+1}=\mathcal D_{\sigma_k}(x^{k+1}+u^k).$$

于是 PnP-ADMM 的基本结构变成

$$x^{k+1}=\operatorname*{argmin}_x\left\{f(x)+\frac{\rho_k}{2}\left\lVert x-(v^k-u^k)\right\rVert_2^2\right\},$$

$$v^{k+1}=\mathcal D_{\sigma_k}(x^{k+1}+u^k),$$

$$u^{k+1}=u^k+x^{k+1}-v^{k+1},$$

其中

$$\sigma_k=\sqrt{\frac{\lambda}{\rho_k}}.$$

从代码结构上看，只改了第二行。但概念上要分清哪些东西保留了，哪些东西已经被换掉。

保留下来的有：

- $f$仍然负责数据一致性；
- $x$步仍然知道前向模型$A$和观测$y$；
- $u$仍然协调数据模块与先验模块的一致性；
- 整个算法仍然采用 ADMM 式的交替模块结构。

被换掉的是：

- 不再要求用户先显式写出一个$g$；
- 不再真的计算$\operatorname{prox}_g$；
- 先验知识直接封装在$\mathcal D_{\sigma}$里。

这给出 PnP 最吸引人的模块化性质：如果更换成像系统，只需要重写$f$对应的数据子问题；如果换一个更好的 denoiser，数据模型不需要跟着改。

但恰恰因为第二步不再明确来自$g$，新的理论问题马上出现了。

---

## 四、为什么“只替换一行”以后，原来 ADMM 的最优性结论不能直接搬过来

经典 ADMM 为什么能说“目标值趋于最优值”“原始残差趋于0”？因为它一开始就在解决一个明确的问题：

$$\min_{x,v}f(x)+\lambda g(v)\qquad\text{s.t.}\qquad x=v.$$

并且在 Boyd 2011 的经典条件下，$f$和$g$是 closed、proper、convex，未增广拉格朗日存在鞍点。

这时每一步都能追溯回同一个显式优化问题，经典结论包括：

- 原始残差趋于0；
- 目标值趋于$p^\star$；
- 对偶变量趋于对偶最优点。

而 PnP 替换后，我们只知道

$$v^{k+1}=\mathcal D_{\sigma_k}(x^{k+1}+u^k).$$

除非额外证明$\mathcal D_{\sigma}$确实是某个合适$g$的 proximal map，否则不能从“它看起来像 prox”推出

$$\exists g:\quad \mathcal D_\sigma=\operatorname{prox}_{\sigma^2 g}.$$

更不能进一步默认整个算法仍在最小化某个我们没写出来的$f+\lambda g$。

这就是 PnP 理论的第一个真正分叉：

> 如果 denoiser 不一定来自显式先验，那么原来的“优化最优点”还是不是正确的收敛目标？

早期一种思路是要求 denoiser 非扩张并具有对称梯度，从而尝试恢复 proximal-map 解释。但对实际非线性 denoiser，这些条件往往很难验证；Chan 2017 甚至给出 NLM 可能出现扩张的反例。

于是论文采取了另一条路线：既然暂时无法证明“最终点是某个未知目标的最优解”，那就先问一个更基础的问题。

---

## 五、如果“最优哪个目标”说不清，那至少能不能先保证迭代不会一直乱跑

这个新问题就是 `fixed-point convergence`。

经典优化的提问方式是：

> $x^k$是否越来越接近某个显式问题的最优解？

PnP 在一般 denoiser 下先把问题降一级：

> $(x^k,v^k,u^k)$这个整个算法状态会不会最终稳定下来？

如果存在

$$\theta^\star=(x^\star,v^\star,u^\star)$$

使

$$x^k\to x^\star,\qquad v^k\to v^\star,\qquad u^k\to u^\star,$$

那么算法至少不会永远振荡或漂移，而会进入稳定状态。

这个结论比“找到显式凸目标的全局最优点”弱，但它适用于更宽的 denoiser 范围。

接下来问题变成：

> 要让一个黑盒 denoiser 驱动的非线性迭代稳定，最少应该限制 denoiser 的什么行为？

如果要求它对任意两个输入都非扩张，条件太强；Chan 2017 改为限制一个更直接的量：`每次 denoiser 最多能把当前输入改多远。`

---

## 六、为什么 bounded denoiser 恰好控制“每轮会不会改得太猛”

Chan 2017 定义 bounded denoiser：存在与$n$、$\sigma$无关的常数$C$，使任意输入$x$满足

$$\frac1n\left\lVert\mathcal D_\sigma(x)-x\right\rVert_2^2\le\sigma^2C.$$

等价地，denoiser 对输入的改变量满足

$$\left\lVert\mathcal D_\sigma(x)-x\right\rVert_2\le\sigma\sqrt{nC}.$$

这个条件和`nonexpansive`不是一回事。

nonexpansive 比较两个输入：

$$\left\lVert\mathcal D(x)-\mathcal D(y)\right\rVert\le\left\lVert x-y\right\rVert.$$

bounded denoiser 只看一个输入和它自己的输出：

$$\left\lVert\mathcal D_\sigma(x)-x\right\rVert.$$

所以 bounded 条件没有要求 denoiser 在整个空间里都不放大任意两点之间的距离，它只要求：`给定去噪强度$\sigma$，一次去噪不能无限制地把图像推走。`

更关键的是，当

$$\sigma\to0$$

时，右边也趋于0，于是

$$\mathcal D_\sigma(x)\to x.$$

也就是说，denoiser 会渐近接近恒等映射。

但现在还差一步。如果$\sigma$始终固定，那么 bounded 条件只告诉我们“每轮最多改一个固定幅度”，这并不能保证无限多轮的总移动距离有限。

例如每一步都移动不超过0.1，并不意味着走一万步以后还待在一个有限小区域内。

所以要从“单步有界”走到“整个序列收敛”，还必须继续让单步影响逐渐减弱。这就是 continuation 为什么会自然出现。

---

## 七、仅有 bounded 还不够：为什么 continuation 必须和它配套出现

PnP-ADMM 令

$$\rho_{k+1}=\gamma_k\rho_k,\qquad\gamma_k\ge1,$$

同时

$$\sigma_k=\sqrt{\frac{\lambda}{\rho_k}}.$$

因此当$\rho_k$增大时，$\sigma_k$会减小：

$$\rho_k\uparrow\quad\Longrightarrow\quad\sigma_k\downarrow\quad\Longrightarrow\quad\mathcal D_{\sigma_k}\text{逐渐接近恒等映射}.$$

同时$x$步中的二次项

$$\frac{\rho_k}{2}\left\lVert x-(v^k-u^k)\right\rVert^2$$

也越来越强，使数据子问题不会轻易远离当前共识。

所以 continuation 的作用可以先理解成：

> 随着迭代推进，让“数据模块”和“去噪模块”每轮允许做的大动作逐渐变小。

Chan 2017 给出两种更新方式。

`单调 continuation`直接取固定$\gamma>1$：

$$\rho_{k+1}=\gamma\rho_k.$$

这时$\rho_k$几何增长，$\sigma_k$几何减小。

更有意思的是`自适应 continuation`。先定义整个状态的单步变化量

$$\Delta_{k+1}=\frac1{\sqrt n}\left(\left\lVert x^{k+1}-x^k\right\rVert_2+\left\lVert v^{k+1}-v^k\right\rVert_2+\left\lVert u^{k+1}-u^k\right\rVert_2\right).$$

然后检查这一轮是否已经收缩得足够快：

- 若$\Delta_{k+1}\ge\eta\Delta_k$，说明状态下降得不够快，于是增大$\rho_k$；
- 若$\Delta_{k+1}<\eta\Delta_k$，其中$0\le\eta<1$，说明状态本身已经按比例收缩，就暂时保持$\rho_k$不变。

这里非常容易误解成：

> 为了收敛，自适应 continuation 也一定要求$\rho_k\to\infty$。

完整证明并不是靠这一句话成立的。真正起作用的是两种不同的稳定机制：

1. 如果这一轮需要增大$\rho_k$，可以证明相邻状态的距离受到$O(1/\sqrt{\rho_k})$控制；反复增大$\rho_k$会让这个上界越来越小。
2. 如果这一轮不增大$\rho_k$，那是因为算法已经检测到$\Delta_{k+1}<\eta\Delta_k$，状态增量本身正在做几何收缩。

也就是说，不管当前走哪一个分支，都在逼着“下一步还能走多远”变小。

这才是 bounded denoiser 与 continuation 真正互相配合的地方：

$$\text{bounded 控制 denoiser 单次最大改变量}\quad+\quad\text{continuation / 自适应判据让改变量持续缩小}.$$

---

## 八、为什么证明最后会落到 Cauchy sequence，而不是重新找一个目标函数

既然现在的目标是证明整个算法状态稳定，最直接的数学对象就是

$$\theta^k=(x^k,v^k,u^k).$$

在这个状态空间里定义距离

$$D(\theta^k,\theta^j)=\frac1{\sqrt n}\left(\left\lVert x^k-x^j\right\rVert_2+\left\lVert v^k-v^j\right\rVert_2+\left\lVert u^k-u^j\right\rVert_2\right).$$

那么前面的$\Delta_{k+1}$就是

$$\Delta_{k+1}=D(\theta^{k+1},\theta^k).$$

证明的核心不再是构造某个隐藏目标函数，而是想办法得到

$$D(\theta^{k+1},\theta^k)\le C'\delta^k,\qquad0<\delta<1.$$

一旦相邻状态的移动量被几何序列控制，就可以把从第$k$步到第$m$步的总移动距离累加起来：

$$D(\theta^m,\theta^k)\le\sum_{j=k}^{m-1}C'\delta^j.$$

右边是几何级数的尾和，当$k\to\infty$时趋于0，因此

$$D(\theta^m,\theta^k)\to0.$$

这正是 Cauchy sequence 的含义。

而有限维欧氏空间是完备的，所以 Cauchy 序列一定收敛，于是得到

$$\theta^k\to\theta^\star.$$

这条证明链非常值得记住，因为它体现了 Chan 2017 和经典 ADMM 完全不同的理论策略：

$$\text{不先寻找隐藏目标}\rightarrow\text{直接控制每轮状态变化}\rightarrow\text{让变化量可求和}\rightarrow\text{Cauchy}\rightarrow\text{fixed point}.$$

---

## 九、收敛以后到底能说什么，又绝对不能多说什么

在论文条件下，可以得到

$$x^k\to x^\star,\qquad v^k\to v^\star,\qquad u^k\to u^\star.$$

并且因为

$$u^{k+1}-u^k=x^{k+1}-v^{k+1},$$

左边趋于0，所以

$$x^{k+1}-v^{k+1}\to0.$$

这说明数据模块和先验模块最终会达成一致。

但下面三件事必须严格分开：

| 问题 | 能否由 Chan 2017 的 fixed-point 结论直接推出 |
| --- | --- |
| 算法状态最终不再变化 | 可以 |
| $x^k-v^k\to0$，两个模块最终一致 | 可以 |
| 极限点一定是某个显式$f+\lambda g$的全局最优解 | 不可以无条件推出 |
| 极限点一定拥有最高 PSNR / 最好视觉质量 | 更不能由收敛定理推出 |

`稳定`、`优化最优`、`重建质量最好`是三个不同层次的问题。

这也是学习 PnP 后必须形成的习惯：每看到一句“算法收敛”，立刻继续问——

> 收敛的是什么量？收敛到什么对象？这个对象和某个显式目标的最优解有没有被真正证明为同一个东西？

---

## 十、为什么 PnP 的停止准则也不能直接照搬经典 ADMM

经典 ADMM 的停止准则与显式问题的最优性条件直接相连。

在一致性分裂$x=v$下，原始残差可写成

$$r^{k+1}=x^{k+1}-v^{k+1},$$

它衡量是否满足$x=v$。

对偶残差则反映相邻$v$更新对对偶可行性的影响；在这一简单分裂下常写成与

$$\rho(v^{k+1}-v^k)$$

成比例的量。

而 Chan 2017 研究的是 fixed point，因此更自然地检查整个状态是否停止变化：

$$\Delta_{k+1}=\frac1{\sqrt n}\left(\lVert\Delta x^k\rVert_2+\lVert\Delta v^k\rVert_2+\lVert\Delta u^k\rVert_2\right).$$

这两个判据回答的不是同一个问题：

- primal/dual residual：离显式 ADMM 最优性条件还有多远；
- $\Delta_{k+1}$：整个 PnP 动力系统离“不再变化”还有多远。

因此不能只因为$\Delta$很小，就把它翻译成“已经满足经典凸优化最优性条件”。

---

## 十一、放回逆问题：为什么换前向模型时主要改$x$步，而 denoiser 可以保持不变

现在回头看 PnP 的工程价值会更清楚。

假设数据项是最常见的线性 Gaussian 模型

$$f(x)=\frac12\lVert Ax-y\rVert_2^2.$$

那么$x$步为

$$x^{k+1}=\operatorname*{argmin}_x\left\{\frac12\lVert Ax-y\rVert_2^2+\frac{\rho_k}{2}\lVert x-(v^k-u^k)\rVert_2^2\right\}.$$

一阶条件给出

$$\left(A^TA+\rho_k I\right)x^{k+1}=A^Ty+\rho_k(v^k-u^k).$$

所以只要成像系统$A$改变，真正需要重新设计或加速的主要是这个数据反演步骤。

而先验模块仍然保持

$$v^{k+1}=\mathcal D_{\sigma_k}(x^{k+1}+u^k).$$

这就是“plug-and-play”比一句“把 prox 换成 denoiser”更深的一层含义：

> 它把`成像物理`和`图像先验`拆成两个算法模块，使二者可以相对独立地升级。

当然，模块独立不等于参数互不影响。$\rho_k$同时影响数据一致性约束强度和$\sigma_k=\sqrt{\lambda/\rho_k}$，因此求解器参数和 denoiser 强度仍然通过迭代机制耦合。

---

## 十二、PnP-ADMM 已经能插 denoiser 了，为什么下一步还会自然走向 PnP-FISTA

到这里 PnP-ADMM 已经解决了最关键的接口问题：

$$\text{显式 prior prox}\longrightarrow\text{任意合适的 denoiser}.$$

但 ADMM 的代价也已经暴露出来。

它每轮需要求解一个$x$子问题：

$$x^{k+1}=\operatorname*{argmin}_x\left\{f(x)+\frac{\rho_k}{2}\lVert x-c_k\rVert^2\right\}.$$

当$f$对应复杂 forward model 时，这一步可能需要线性方程求解、内层迭代或专门结构加速。Gavaskar 2021 也明确指出：PnP-FISTA 往往实现更简单，但要求$f$光滑；PnP-ADMM 对$f$更通用，却可能有更昂贵的数据子问题。

所以如果我们的数据项本来就是光滑的，例如

$$f(x)=\frac12\lVert Ax-y\rVert_2^2,$$

就会自然产生一个新问题：

> 既然[[1 FISTA（讲解版）|FISTA]]每轮只需要一次梯度、一次 prox 和一个动量外推，那么能不能把它的 prox 也换成 denoiser，从而得到一个更轻量的 PnP 骨架？

答案就是 [[0 PnP-FISTA（讲解版）]]。

但这次不能简单得出“ADMM 换成 FISTA 就更快”。因为 prox 被 denoiser 替换以后，FISTA 原来的$O(1/k^2)$究竟还在不在、非对称 denoiser 应该在哪种几何里理解、以及“目标值快”和“图像迭代点快”是不是一回事，都必须重新回答。

---

## 十三、把整条问题链收回来：PnP-ADMM 为什么会自然出现

现在可以从头把 PnP-ADMM 重新走一遍。

最开始，我们有显式逆问题

$$\min_x f(x)+\lambda g(x).$$

第一个问题是：`现代 denoiser 是一个图像到图像的映射，不能直接塞进标量目标函数，怎么办？`

于是需要一个能够暴露“先验处理接口”的算法骨架。

ADMM 通过

$$x=v$$

把数据职责和先验职责拆开：

$$\text{数据$x$步}\longleftrightarrow\text{先验$v$步}\longleftrightarrow\text{对偶$u$协调}.$$

第二个问题是：`为什么$v$步能直接接 denoiser？`

因为它经过重写后恰好是

$$g(v)+\frac1{2\sigma^2}\lVert v-r\rVert^2,$$

也就是标准的“结构 + 靠近输入”的 proximal denoising 形式。

第三个问题是：`把 prox 换成$\mathcal D_\sigma$以后，经典 ADMM 理论还能不能原封不动使用？`

一般不能，因为$\mathcal D_\sigma$未必是某个显式凸$g$的 prox。

第四个问题于是变成：`既然最优哪个目标暂时说不清，至少能不能保证迭代状态稳定？`

Chan 2017 用 bounded denoiser 控制单次去噪改变量，再用 continuation / 自适应收缩机制让相邻状态变化越来越小：

$$\text{bounded}\rightarrow\text{状态增量几何控制}\rightarrow\text{可求和}\rightarrow\text{Cauchy}\rightarrow\text{fixed point}.$$

最终真正应该记住的不是三条孤立迭代式，而是：

> `PnP-ADMM 的核心创新不是“ADMM + 一个去噪器”这么表面。ADMM 先通过变量分裂把显式先验的 proximal step 暴露成图像到图像的接口，PnP 才有位置把它替换成 denoiser；而一旦这样替换，理论目标就从经典显式优化最优性转向了必须重新证明的 fixed-point 稳定性。`
