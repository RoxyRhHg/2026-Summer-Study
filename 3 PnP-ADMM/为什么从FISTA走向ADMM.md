---
tags:
  - 暑假学习
  - 学习路径
  - FISTA
  - ADMM
  - 近端算法
  - PnP-ADMM
  - 逆问题
created: 2026-09-04
updated: 2026-09-04
aliases:
  - FISTA到ADMM的过渡
  - ADMM为什么自然出现
---

# 为什么从 FISTA 走向 ADMM

> [!cite] Zotero 资料基线
> - [[Beck 等 - 2009 - A Fast Iterative Shrinkage-Thresholding Algorithm for Linear Inverse Problems|Beck 与 Teboulle 2009：FISTA]] · [Zotero](zotero://select/library/items/EJXDP5IE)
> - [[Boyd 等 - 2011 - Distributed Optimization and Statistical Learning via the Alternating Direction Method of Multipliers|Boyd 等 2011：ADMM]] · [Zotero](zotero://select/library/items/B7AYLH4I)

## 先把最容易产生的误解说清楚：ADMM 不是 FISTA 的“下一代”

学完 [[1 FISTA（讲解版）]] 后，我们已经得到一条很完整的路线：

$$\min_x F(x)=f(x)+g(x),$$

其中 $f$ 光滑、$g$ 可以不可微但 prox 易算。FISTA 每轮做的核心事情是

$$y_k\longrightarrow y_k-\gamma\nabla f(y_k)\longrightarrow\operatorname{prox}_{\gamma g}(\cdot),$$

再利用历史点构造下一轮的外推位置。

所以刚进入 PnP-ADMM 时，一个非常自然的问题就是：

> FISTA 已经把“数据项”和“先验项”分开处理了，而且本身已经有一个 $\operatorname{prox}_g$ 接口。那为什么不直接把这个 prox 换成 denoiser？为什么还要突然学 ADMM？

这个疑问是对的，因为

$$\text{FISTA}\rightarrow\text{ADMM}$$

本来就不是“上一种算法继续升级成下一种算法”的关系。

更准确的知识结构是：FISTA 和 ADMM 都从同一个复合问题

$$f(x)+g(x)$$

出发，但采用了两种不同的`拆问题方式`。

- `FISTA / 近端梯度`：保留一个变量 $x$，把 $f$ 线性化，用梯度处理 $f$，用 prox 完整处理 $g$。
- `ADMM`：不要求一定把 $f$ 线性化，而是先复制变量，让 $f$ 和 $g$ 各自拥有自己的变量，再通过约束强迫两个副本最终一致。

所以从 FISTA 走向 ADMM，真正的问题不是：

> FISTA 哪里不够快，所以换一个更快算法？

而是：

> 除了“梯度 + prox”以外，还有没有另一种办法，把数据项和先验项拆成两个更独立的计算模块？

这才是 ADMM 应该出现的位置。

---

## 一、FISTA 已经在“拆” $f$ 和 $g$，ADMM 到底还想拆什么

先看 FISTA。

对于

$$\min_x f(x)+g(x),$$

FISTA 的近端梯度核心是

$$z_k=y_k-\gamma\nabla f(y_k),$$

$$x_k=\operatorname{prox}_{\gamma g}(z_k).$$

它实际上做了这样的职责划分：

$$\text{光滑数据项 }f\xrightarrow{\nabla f}\text{梯度步},$$

$$\text{结构项 }g\xrightarrow{\operatorname{prox}_g}\text{结构修正}.$$

但要注意：这是一种`不对称分工`。

$f$ 被局部线性化，只保留梯度信息；$g$ 则完整留在 proximal 子问题里。因此 FISTA 很适合

$$f\text{ 光滑},\qquad g\text{ 的 prox 容易求}$$

这一类问题。

ADMM 想做的是另一件事：

> 能不能让 $f$ 和 $g$ 都保留自己的完整形式，各自解决自己的子问题，而不是要求其中一边必须先被线性化？

这就需要比“把一个公式拆成两步计算”更彻底的分离：`把变量本身也拆开。`

两种思路可以先并排看：

| 思路 | FISTA / 近端梯度 | ADMM |
| --- | --- | --- |
| 原问题 | $f(x)+g(x)$ | $f(x)+g(x)$ |
| 怎样拆 | 一个变量，两种操作 | 复制变量，两个子问题 |
| $f$ 怎么处理 | 梯度 / 局部线性化 | 直接解 $f$ 子问题 |
| $g$ 怎么处理 | prox | 直接解 $g$ 子问题，常常也是 prox |
| 怎样保证二者一致 | 两步共同定义同一个 $x_{k+1}$ | 显式约束两个变量最终一致 |
| 额外状态 | 动量 / 外推状态 | 对偶变量 / 共识误差 |

所以 ADMM 不是把 FISTA 的某个步骤“加速”了，而是换了一种拆解整个优化问题的方式。

---

## 二、如果想让 $f$ 和 $g$ 各自算自己的，最直接的办法就是复制变量

仍然从

$$\min_x f(x)+g(x)$$

出发。

现在新建一个变量 $v$，把 $g$ 改成作用在 $v$ 上：

$$\min_{x,v} f(x)+g(v)\qquad\text{s.t.}\qquad x=v.$$

为什么这样做没有改变原问题？

因为只要约束 $x=v$ 最终满足，就有

$$f(x)+g(v)=f(x)+g(x).$$

所以这里不是在真的创造两个不同的最终答案，而是在迭代过程中让同一个未知图像暂时拥有两个“工作副本”：

- $x$：只面对数据项 $f$；
- $v$：只面对先验项 $g$；
- $x=v$：要求两个模块最后必须对同一个结果达成一致。

这一步叫做`variable splitting`。

它真正换来的东西是：原来紧紧绑在一个变量上的两个函数，现在被放到了两个不同变量上。

但问题也马上来了。

> 拆开当然容易，可拆开以后怎么保证 $x$ 和 $v$ 不会各做各的？

这就把我们带到约束优化和 Lagrange multiplier。

---

## 三、只写普通 Lagrangian 为什么还不够，于是出现 augmented Lagrangian

对于约束

$$x-v=0,$$

普通 Lagrangian 是

$$L_0(x,v,y)=f(x)+g(v)+y^T(x-v),$$

其中 $y$ 是对偶变量。

这里的线性项

$$y^T(x-v)$$

会根据当前 $x-v$ 的方向对两个变量施加修正，因此 $y$ 可以理解成一份不断更新的`一致性压力`。

但 Boyd 2011 在回顾 dual ascent 时指出，普通 dual ascent 的 $x$ 最小化在不少问题上会出现困难，例如子问题可能没有良好的有限最小点，因此其适用条件并不够稳健。

一个自然改进是在约束残差上再加入平方惩罚：

$$\frac{\rho}{2}\lVert x-v\rVert_2^2,\qquad \rho>0.$$

得到`增广拉格朗日`：

$$L_\rho(x,v,y)=f(x)+g(v)+y^T(x-v)+\frac{\rho}{2}\lVert x-v\rVert_2^2.$$

这个式子里有两种不同的一致性机制：

- $y^T(x-v)$：对偶变量根据分歧的方向不断累计修正；
- $(\rho/2)\lVert x-v\rVert^2$：当前这一轮只要两个副本相差太远，就立刻付出二次代价。

因此 augmented Lagrangian 不是单纯“再加一个正则项”，而是把

$$\text{乘子协调}\quad+\quad\text{当前残差惩罚}$$

放在了一起。

如果每轮都联合求

$$ (x^{k+1},v^{k+1})=\operatorname*{argmin}_{x,v}L_\rho(x,v,y^k),$$

再更新对偶变量

$$y^{k+1}=y^k+\rho(x^{k+1}-v^{k+1}),$$

这对应 `method of multipliers` 的思路。

它比普通 dual ascent 稳健，但又产生了一个新矛盾：我们前面费力把 $f$ 和 $g$ 拆开，结果现在却要求在一个联合最小化里同时求 $x$ 和 $v$，两个模块又重新耦在一起了。

Boyd 2011 对这一点说得很直接：augmented Lagrangian 改善了收敛性质，但联合最小化会破坏原本希望得到的 decomposition。

于是 ADMM 真正要解决的问题出现了：

> 能不能保留 augmented Lagrangian 带来的稳定性，却不要每轮把 $x$ 和 $v$ 联合求到底？

---

## 四、ADMM 的关键不是“多了三条公式”，而是把联合最小化改成一次交替更新

ADMM 的做法非常直接。

不再联合求

$$\min_{x,v}L_\rho(x,v,y^k),$$

而是一轮中依次做：

$$x^{k+1}=\operatorname*{argmin}_x L_\rho(x,v^k,y^k),$$

$$v^{k+1}=\operatorname*{argmin}_v L_\rho(x^{k+1},v,y^k),$$

$$y^{k+1}=y^k+\rho(x^{k+1}-v^{k+1}).$$

也就是：

1. 固定上一轮 $v^k$，只让数据变量 $x$ 更新；
2. 用刚得到的 $x^{k+1}$，只让先验变量 $v$ 更新；
3. 看这两个模块还差多少，用 $y$ 把这种分歧记下来并传给下一轮。

这就是名字里的 `Alternating Direction`：原来需要同时处理的两个 primal variable，现在采用顺序、交替的方式各走一步。

Boyd 2011 也可以把它理解成：method of multipliers 本来需要对 $x,v$ 做联合最小化，而 ADMM 只做一次顺序的 Gauss-Seidel 更新，从而重新获得 decomposition。

所以 ADMM 的核心问题意识可以压成一句话：

> `既想利用 augmented Lagrangian 强迫不同模块达成一致，又不想因为联合最小化把已经拆开的模块重新绑死，于是把一次联合最小化改成两个交替子问题。`

到这里，三步 ADMM 才不是凭空需要背的公式。

---

## 五、为什么又要引入 scaled dual variable $u$：只是把“分歧账本”写得更干净

PnP 文献里通常不直接使用 $y$，而是定义

$$u=\frac{y}{\rho}.$$

原因可以从增广项直接看出来。令

$$r=x-v,$$

那么

$$y^Tr+\frac{\rho}{2}\lVert r\rVert_2^2$$

配方得到

$$\frac{\rho}{2}\left\lVert r+\frac{y}{\rho}\right\rVert_2^2-\frac{1}{2\rho}\lVert y\rVert_2^2.$$

定义 $u=y/\rho$ 后，与当前优化变量无关的最后一项可以在各子问题中忽略，于是对于约束 $x=v$，scaled ADMM 写成

$$x^{k+1}=\operatorname*{argmin}_x\left\{f(x)+\frac{\rho}{2}\lVert x-v^k+u^k\rVert_2^2\right\},$$

$$v^{k+1}=\operatorname*{argmin}_v\left\{g(v)+\frac{\rho}{2}\lVert x^{k+1}-v+u^k\rVert_2^2\right\},$$

$$u^{k+1}=u^k+x^{k+1}-v^{k+1}.$$

最后一行尤其重要：

$$u^k=u^0+\sum_{j=1}^k(x^j-v^j).$$

所以 $u$ 可以先非常直观地理解成`历史分歧账本`。

- 如果 $x$ 和 $v$ 已经一致，$u$ 不再因为当前残差继续变化；
- 如果两个模块长期存在偏差，这个偏差会累积进 $u$，并改变后面的两个子问题。

这样就能把 ADMM 一轮看成：

$$\text{数据模块更新}\rightarrow\text{先验模块更新}\rightarrow\text{累计并纠正二者分歧}.$$

---

## 六、到这里才真正接回 FISTA 学过的 prox：ADMM 的 $v$ 步为什么像去噪

现在专门看

$$v^{k+1}=\operatorname*{argmin}_v\left\{g(v)+\frac{\rho}{2}\lVert v-(x^{k+1}+u^k)\rVert_2^2\right\}.$$

这和 [[近端梯度为什么自然出现]] 中已经反复见过的结构完全一样：

$$\text{结构代价}\quad+\quad\text{不要离某个输入点太远}.$$

按照 proximal operator 的定义，令 $\gamma=1/\rho$，就有

$$v^{k+1}=\operatorname{prox}_{\gamma g}(x^{k+1}+u^k).$$

这一步非常关键，因为它第一次把两条看似不同的路线重新接上了：

- FISTA 中，prox 出现在`梯度步之后`；
- ADMM 中，prox 出现在`变量分裂后的先验子问题`里。

两者都产生了同一种计算接口：

$$\text{输入一个当前图像}\longrightarrow\operatorname{prox}_g\longrightarrow\text{输出一个更符合先验的图像}.$$

也正因为这样，进入图像恢复以后，一个更大胆的问题才真正有了落脚点：

> 既然这一整步的计算语义已经如此像“输入一张图，输出一张更干净、更符合图像先验的图”，那它是否一定非得来自一个手工写出来的 $g$？能不能直接换成现成 denoiser？

这就是 [[0 PnP-ADMM（讲解版）]] 的真正起点。

---

## 七、那为什么学习路线先走 PnP-ADMM，而不是直接 PnP-FISTA

现在可以回答最开始那个疑问了。

从数学形式上说，当然可以直接从 FISTA 写出

$$\operatorname{prox}_{\gamma g}\longrightarrow\mathcal D$$

并得到 PnP-FISTA。PnP-FISTA 并不在逻辑上依赖“必须先运行一遍 PnP-ADMM”。

但在学习顺序上，PnP-ADMM 更适合作为第一次理解 Plug-and-Play 的入口，因为 ADMM 把三个角色暴露得非常彻底：

$$x\text{：数据模块},\qquad v\text{：先验模块},\qquad u\text{：模块间共识协调}.$$

于是替换

$$\operatorname{prox}_g\longrightarrow\mathcal D_\sigma$$

时，我们很容易看清楚：`真正被换掉的只有先验模块，forward model 仍然留在数据模块。`

等 PnP 这个思想本身建立以后，再回头看熟悉的 FISTA，就会自然产生下一问：

> 如果 Plug-and-Play 的本质不是“必须用 ADMM”，而是“把显式 proximal prior 换成 denoiser”，那么 FISTA 里的 prox 当然也可以被替换。这样是不是能保留更轻量的一阶数据步和历史外推？

这时才进入 [[0 PnP-FISTA（讲解版）]]。

所以合理的学习路线不是

$$\text{FISTA}\rightarrow\text{更高级的 ADMM}\rightarrow\text{更高级的 PnP-FISTA},$$

而是

$$\begin{aligned}
&\text{显式复合优化 }f+g\\
&\quad\downarrow\\
&\text{近端思想}\\
&\quad\swarrow\qquad\searrow\\
&\text{FISTA}\qquad\text{variable splitting / ADMM}\\
&\qquad\qquad\downarrow\\
&\qquad\text{先验子问题暴露为 prox 接口}\\
&\qquad\qquad\downarrow\\
&\qquad\text{PnP-ADMM：建立 Plug-and-Play 思想}\\
&\qquad\qquad\downarrow\\
&\qquad\text{回到 FISTA：同一个 prox 接口也能替换}\\
&\qquad\qquad\downarrow\\
&\qquad\text{PnP-FISTA}.
\end{aligned}$$

真正需要记住的是：

> `FISTA 和 ADMM 是处理同一类复合问题的两种不同分裂框架。学习 ADMM 不是因为 FISTA“失败了”，而是为了看到另一种更彻底的模块化方式：复制变量、分别处理数据与先验，再通过对偶变量维持共识。正是这种结构把先验步骤显式暴露成一个图像到图像的 prox 接口，从而让 PnP-ADMM 的 denoiser 替换变得自然。`
