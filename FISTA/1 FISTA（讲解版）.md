---
tags:
  - 暑假学习
  - 学习路径
  - FISTA
  - 近端梯度
  - 逆问题
created: 2026-08-30
updated: 2026-09-04
aliases:
  - FISTA讲解版
  - 近端梯度与FISTA讲解
---

# FISTA：把“梯度下降”一步一步改造成加速近端梯度

## 先抓住整件事：FISTA到底在解决什么

FISTA不是一个应该从三条递推公式开始背的算法。

它真正要解决的是这一类`复合凸优化`问题：

$$\min_x F(x)=f(x)+g(x).$$

这里通常有非常明确的分工：

- $f(x)$：`光滑项`，可以求梯度；
- $g(x)$：`结构项 / 正则项`，可能不可微，但它的 proximal operator 比较容易算。

例如图像重建或稀疏恢复里最经典的模型：

$$\min_x\frac12\lVert Ax-b\rVert_2^2+\lambda\lVert x\rVert_1.$$

其中

$$f(x)=\frac12\lVert Ax-b\rVert_2^2,\qquad g(x)=\lambda\lVert x\rVert_1.$$

第一项衡量“重建结果经过成像模型以后，和真实观测差多少”；第二项则要求$x$具有稀疏结构。

如果整个目标都光滑，我们会想到梯度下降：

$$x_{k+1}=x_k-\gamma\nabla F(x_k).$$

但现在$g(x)=\lambda\lVert x\rVert_1$在某些点不可微，因此$\nabla F$不一定存在。

于是问题变成：

> 能不能继续使用$f$的梯度信息，同时又完整保留$g$所表达的结构？

答案就是`近端梯度法`。而FISTA只是再往前走一步：

> 近端梯度已经能做了，但$O(1/k)$还是慢；能不能基本不增加单轮成本，把函数值收敛率提到$O(1/k^2)$？

因此整条理解链应该是：

$$\text{梯度下降}\longrightarrow\text{局部二次模型}\longrightarrow\text{近端梯度}\longrightarrow\text{ISTA}\longrightarrow\text{FISTA}. $$

---

## 一、先把梯度下降重新理解一遍

### 1. 为什么不直接最小化一阶 Taylor 展开

设$f$光滑，在当前点$y$附近，一阶 Taylor 近似为

$$f(x)\approx f(y)+\langle\nabla f(y),x-y\rangle.$$

这一项告诉我们当前最重要的局部信息：

- $\nabla f(y)$告诉我们函数增长最快的方向；
- $-\nabla f(y)$因此是最自然的下降方向。

但如果真的去最小化这个线性模型，会出问题。

令

$$x-y=-t\nabla f(y),\qquad t>0,$$

则线性项变成

$$\langle\nabla f(y),x-y\rangle=-t\lVert\nabla f(y)\rVert_2^2.$$

当$t\to\infty$时，它会一直跑向$-\infty$。

也就是说，一阶模型会给出一个荒谬建议：

> 既然负梯度方向能下降，那就沿这个方向一直冲到无穷远。

问题不在梯度方向错了，而在于：`一阶 Taylor 模型只在当前点附近可信。`

### 2. 所以给“离当前点太远”加成本

加入一个二次项：

$$f(y)+\langle\nabla f(y),x-y\rangle+\frac L2\lVert x-y\rVert_2^2.$$

现在局部模型里有两股力量：

- 一阶项：希望沿下降方向走；
- 二次项：走得越远，代价越大。

于是不会再无限跑远，而会找到一个有限的折中点。

为什么恰好是这个二次项，而不是随便拍脑袋加出来的，见[[为什么在线性化模型中加入二次惩罚]]。这里先保留最关键的一层：如果$\nabla f$是[[Lipschitz 条件]]，那么有

$$f(x)\le f(y)+\langle\nabla f(y),x-y\rangle+\frac L2\lVert x-y\rVert_2^2.$$

所以这个二次模型不只是“方便算”，还能够罩住真实$f(x)$。

### 3. 最小化它，梯度下降自己就出来了

忽略与$x$无关的$f(y)$，对$x$求导：

$$\nabla f(y)+L(x-y)=0.$$

解得

$$x=y-\frac1L\nabla f(y).$$

这就是步长为$1/L$的梯度下降。

所以梯度下降可以换一种更重要的理解：

> 每一轮不是简单“沿负梯度走一步”，而是在当前点构造一个容易求解的局部二次模型，再把这个局部模型的最小点当作下一步。

这个视角一旦建立，proximal gradient几乎是自然长出来的。

---

## 二、从“局部模型”到近端梯度：主推导放在哪里

第一部分已经得到一个更重要的梯度下降视角：

> `梯度下降 = 最小化“一阶 Taylor 模型 + 二次稳定项”。`

现在回到复合目标

$$F(x)=f(x)+g(x),$$

其中$f$光滑，$g$可能不可微。

这里最自然的设计不是对$f+g$一起求普通梯度，而是：

- $f$能求梯度，所以只对$f$做局部线性化；
- $g$可能不可微，而且往往正是我们想保留的结构，因此把$g(x)$完整留下。

于是得到局部模型

$$Q_L(x,y)=f(y)+\langle\nabla f(y),x-y\rangle+\frac L2\lVert x-y\rVert_2^2+g(x).$$

下一步仍然沿用第一部分的原则：`把局部模型的最小点当作下一步`。

$$x^+=\operatorname*{argmin}_x Q_L(x,y).$$

去掉与$x$无关的$f(y)$：

$$x^+=\operatorname*{argmin}_x\left\{g(x)+\langle\nabla f(y),x-y\rangle+\frac L2\lVert x-y\rVert_2^2\right\}.$$

把后两项配方，并记

$$z=y-\frac1L\nabla f(y),$$

就得到

$$x^+=\operatorname*{argmin}_x\left\{g(x)+\frac L2\lVert x-z\rVert_2^2\right\}.$$

讲解版在这里应该先停一下，因为真正容易卡住的地方已经不是配方，而是：

> 为什么这个`argmin`突然可以写成一个$\operatorname{prox}$？prox这个定义又为什么偏偏长成这样？[[近端梯度为什么自然出现]]。

这里只保留主线。

对于任意输入$v$，把这种反复出现的“结构项 + 不要离输入太远”的子问题定义成

$$\operatorname{prox}_{\gamma g}(v)=\operatorname*{argmin}_x\left\{g(x)+\frac1{2\gamma}\lVert x-v\rVert_2^2\right\}.$$

取$\gamma=1/L$、$v=z$，于是上面的局部模型最小点就是

$$x^+=\operatorname{prox}_{g/L}(z)=\operatorname{prox}_{g/L}\left(y-\frac1L\nabla f(y)\right).$$

到这里，“近端梯度”这个名字才真正有意义：

$$\underbrace{y-\frac1L\nabla f(y)}_{\text{gradient：处理光滑 }f}
\quad\longrightarrow\quad
\underbrace{\operatorname{prox}_{g/L}(\cdot)}_{\text{proximal：处理结构 }g}.$$

因此近端梯度不是突然多出来的一套公式，而是上一节局部二次模型在$g$不可微时的自然延伸。

---

## 三、prox不是额外的“魔法步骤”：它只是把结构重新放回梯度步

上一节已经得到

$$z=y-\frac1L\nabla f(y),$$

$$x^+=\operatorname{prox}_{g/L}(z).$$

理解这两行时，最容易产生一个错觉：先做完梯度下降以后，又额外加了一个神秘的prox步骤。

其实不是。

$z$只是光滑项$f$单独给出的`临时建议`；真正的下一步$x^+$还必须同时考虑$g$。prox做的事情就是：

> 以梯度步结果$z$为中心，在“别离$z$太远”和“让$g(x)$更符合结构要求”之间重新平衡。

详细的问题意识见 [[近端梯度为什么自然出现#10. 为什么叫“近端”：这个定义到底在做什么？|为什么叫近端]]。这里用两个极端情况建立直觉。

如果没有结构项，$g(x)\equiv0$，那么

$$\operatorname{prox}_{\gamma g}(z)=z.$$

也就是说，prox什么都不改，近端梯度自动退化成普通梯度下降。

如果$g$是凸集$C$的指示函数

$$g(x)=\begin{cases}0,&x\in C,\\+\infty,&x\notin C,\end{cases}$$

那么

$$\operatorname{prox}_{\gamma g}(z)=\operatorname*{argmin}_{x\in C}\frac12\lVert x-z\rVert_2^2,$$

它就是把$z$投影回可行集$C$。

所以prox可以先理解成一种更一般的`结构修正器 / 广义投影`：

- 没有结构要求时，不修正；
- 是硬约束时，变成投影；
- 是$\ell_1$稀疏正则时，下一节会看到它变成软阈值。

这正好把第二部分和第四部分接起来：

$$\text{局部二次模型}\rightarrow\text{prox抽象}\rightarrow\text{具体$g$决定prox长什么样}. $$

下面不再重复“prox为什么定义成这样”，而只研究最重要的具体情形：

$$g(x)=\lambda\lVert x\rVert_1.$$

---

## 四、当$g(x)=\lambda\lVert x\rVert_1$时，为什么出现软阈值

考虑

$$g(x)=\lambda\lVert x\rVert_1.$$

prox变成

$$\operatorname{prox}_{\gamma\lambda\lVert\cdot\rVert_1}(v)
=\operatorname*{argmin}_x\left\{\gamma\lambda\lVert x\rVert_1+\frac12\lVert x-v\rVert_2^2\right\}.$$

因为$\ell_1$范数可以拆成坐标之和

$$\lVert x\rVert_1=\sum_i\lvert x_i\rvert,$$

所以每个坐标都可以单独求解。只看一个标量：

$$\min_u\left\{\gamma\lambda\lvert u\rvert+\frac12(u-v)^2\right\}.$$

分三种情况。

若$u>0$，则$\lvert u\rvert=u$，求导得

$$\gamma\lambda+u-v=0,$$

所以

$$u=v-\gamma\lambda.$$

它只有在$v>\gamma\lambda$时才真的满足$u>0$。

若$u<0$，则$\lvert u\rvert=-u$，求导得

$$-\gamma\lambda+u-v=0,$$

所以

$$u=v+\gamma\lambda,$$

它要求$v<-\gamma\lambda$。

剩下的情况就是

$$\lvert v\rvert\le\gamma\lambda,$$

最优解落在$u=0$。

因此

$$S_{\gamma\lambda}(v)=\operatorname{sign}(v)\max\{\lvert v\rvert-\gamma\lambda,0\}.$$

这就是`软阈值`。

它的动作不是简单“把小数删掉”：

- $\lvert v\rvert\le\gamma\lambda$：直接压成0；
- $v>\gamma\lambda$：向0缩小$\gamma\lambda$；
- $v<-\gamma\lambda$：同样向0缩小$\gamma\lambda$。

因此$\ell_1$正则之所以会产生稀疏，是因为prox本身就在持续把系数往0拉，并把一部分系数直接压成0。

---

## 五、ISTA其实就是“反复做近端梯度”

现在把模型换成

$$\min_x\frac12\lVert Ax-b\rVert_2^2+\lambda\lVert x\rVert_1.$$

光滑项是

$$f(x)=\frac12\lVert Ax-b\rVert_2^2.$$

对它求梯度：

$$\nabla f(x)=A^T(Ax-b).$$

梯度Lipschitz常数可以取

$$L=\lVert A\rVert_2^2=\lambda_{\max}(A^TA).$$

于是一步近端梯度就是

$$x_{k+1}=S_{\lambda/L}\left(x_k-\frac1L A^T(Ax_k-b)\right).$$

这就是经典`ISTA`。

所以ISTA并不是另一套独立思想，而只是

$$\text{proximal gradient}+\ell_1\text{正则}. $$

一轮ISTA可以按物理意义拆成：

$$Ax_k\longrightarrow Ax_k-b\longrightarrow A^T(Ax_k-b)\longrightarrow x_k-\frac1L A^T(Ax_k-b)\longrightarrow S_{\lambda/L}(\cdot).$$

逐项解释：

- $Ax_k$：当前$x_k$会产生什么预测观测；
- $Ax_k-b$：预测和真实数据的残差；
- $A^T(Ax_k-b)$：把观测域误差传回未知量空间；
- 梯度步：先改善数据一致性；
- 软阈值：再把结果往稀疏先验方向拉。

所以ISTA一轮最核心的结构是：

$$\text{数据一致性更新}\longrightarrow\text{先验结构更新}. $$

---

## 六、步长$1/L$为什么不能随便取

前面提到，只要$\nabla f$是$L_f$-Lipschitz，且选

$$L\ge L_f,$$

就有二次上界

$$f(x)\le f(y)+\langle\nabla f(y),x-y\rangle+\frac L2\lVert x-y\rVert_2^2.$$

因为真正的梯度步长是

$$\gamma=\frac1L,$$

所以

$$L\uparrow\Longrightarrow\gamma\downarrow.$$

这很符合直觉：函数可能弯得越厉害，就越不能相信当前梯度走太远。

如果$L_f$不知道，可以用`backtracking`。

做法不是先精确求$L_f$，而是：

- 先猜一个$L$；
- 算出当前近端步$p_L(y)$；
- 检查这个$L$构造的二次模型是否真的够高；
- 不够就增大$L$，也就是减小步长，再试。

因此backtracking本质上是在动态寻找一个“足够安全的局部曲率上界”。

---

## 七、ISTA已经能求解了，为什么还需要FISTA

对一般凸复合问题，ISTA的经典函数值收敛率是

$$F(x_k)-F(x^*)=O(1/k).$$

这里一定要注意：说的是`目标函数值误差`，不是直接说

$$\lVert x_k-x^*\rVert_2=O(1/k).$$

两者不是一回事。

如果想让

$$F(x_k)-F(x^*)\le\varepsilon,$$

那么$O(1/k)$意味着大约需要

$$k=O(1/\varepsilon)$$

轮。

但每一轮ISTA本身已经很合理：

- 算一次梯度；
- 算一次prox。

真正昂贵的通常就是这些操作。

所以最吸引人的问题是：

> 能不能不增加第二次梯度、不增加第二次prox，只改变“下一轮从哪里开始算”，就让整体更快？

这就是FISTA的核心。

---

## 八、FISTA只改了一件关键的事：先外推，再做原来的近端梯度

ISTA是

$$x_k=p_L(x_{k-1}),$$

其中

$$p_L(y)=\operatorname{prox}_{g/L}\left(y-\frac1L\nabla f(y)\right).$$

也就是说，每一轮都从刚刚得到的$x_{k-1}$原地出发。

FISTA改成

$$x_k=p_L(y_k),$$

其中$y_k$不是上一轮结果，而是一个`外推点`。

标准固定步长FISTA写成：

$$x_k=\operatorname{prox}_{g/L}\left(y_k-\frac1L\nabla f(y_k)\right),$$

$$t_{k+1}=\frac{1+\sqrt{1+4t_k^2}}2,$$

$$y_{k+1}=x_k+\frac{t_k-1}{t_{k+1}}(x_k-x_{k-1}),$$

初始化

$$y_1=x_0,\qquad t_1=1.$$

先不要急着背。把最后一式改写成

$$y_{k+1}=x_k+\beta_k(x_k-x_{k-1}),$$

其中

$$\beta_k=\frac{t_k-1}{t_{k+1}}.$$

现在意义就很清楚了：

- $x_k$：当前已经走到的位置；
- $x_k-x_{k-1}$：最近一次移动方向；
- $\beta_k(x_k-x_{k-1})$：沿最近趋势再向前预测一点；
- $y_{k+1}$：下一轮真正计算梯度和prox的位置。

所以FISTA不是把prox改得更复杂，而是：

> 先利用历史方向猜一个更靠前的位置，再在那个位置做同样的一次近端梯度。

这就是`Nesterov型加速 / 惯性外推`最直观的一层。

但“惯性”只是帮助理解的类比。真正的加速并不是随便加一个动量系数，而是要非常特殊地设计$\beta_k$。

---

## 九、为什么需要$t_k$，那个根号公式到底从哪来

FISTA里最像“魔法”的公式是

$$t_{k+1}=\frac{1+\sqrt{1+4t_k^2}}2.$$

更好的学习顺序不是先背这个根号，而是先看它满足什么：

$$t_{k+1}^2-t_{k+1}=t_k^2.$$

验证很简单。由

$$2t_{k+1}-1=\sqrt{1+4t_k^2}$$

两边平方：

$$4t_{k+1}^2-4t_{k+1}+1=1+4t_k^2,$$

于是

$$t_{k+1}^2-t_{k+1}=t_k^2.$$

真正重要的是这个恒等式，而不是根号本身。

为什么？因为FISTA证明需要把第$k$轮和第$k+1$轮的误差按$t_k^2$加权后接起来。若权重之间没有这种关系，很多项无法恰好抵消。

所以更接近真实逻辑的理解是：

> 证明希望权重满足$t_{k+1}^2-t_{k+1}=t_k^2$，把它当一元二次方程求正根，才得到那个看起来奇怪的递推式。

换句话说，根号公式是`证明结构的结果`，不是一个凭直觉猜出来的经验参数。

---

## 十、FISTA为什么能从$O(1/k)$变成$O(1/k^2)$：先只看证明骨架

第一次学FISTA，不建议一上来就陷进完整长证明。先抓四层。

### 1. 一次近端步能给出一个“函数误差—距离变化”不等式

设

$$x^+=p_L(y).$$

近端步最优性条件是

$$0\in\partial g(x^+)+\nabla f(y)+L(x^+-y).$$

这意味着存在$s\in\partial g(x^+)$满足

$$s=-\nabla f(y)-L(x^+-y).$$

再结合$f$和$g$的凸性，以及二次模型上界，可以推出一个关键不等式：

$$F(x^+)-F(x)\le\frac L2\left(\lVert x-y\rVert_2^2-\lVert x-x^+\rVert_2^2\right).$$

先不用马上重推所有代数，先读懂它在说什么。

右边是两个距离平方之差：

- 更新前，比较点$x$离当前出发点$y$有多远；
- 更新后，比较点$x$离新点$x^+$有多远。

如果令$x=x^*$，就把当前的目标函数误差和“相对于最优点的几何变化”联系起来了。

这就是整个加速证明真正的地基。

### 2. ISTA基本是“一轮一轮”使用这个不等式

ISTA每次从$x_{k-1}$出发，因此每一轮都能得到类似的下降关系。

把这些关系累积起来，就得到$O(1/k)$。

### 3. FISTA不再逐轮孤立处理，而是把连续两轮按$t_k$加权

FISTA希望构造一个同时包含

- 函数值误差；
- 某个距离平方；

的`势函数`，并让它跨迭代递推。

一个典型结构是

$$\text{函数误差权重}\times\big(F(x_k)-F(x^*)\big)+\text{某个距离平方}.$$

为了让第$k$轮和第$k+1$轮拼接时杂项消掉，就需要前面那个关系

$$t_{k+1}^2-t_{k+1}=t_k^2.$$

所以$t_k$与外推系数的设计，本质上是在服务于这个势函数递推。

### 4. 最后$t_k$本身大约按$k$增长

可以证明

$$t_k\ge\frac{k+1}{2}.$$

而势函数分析最终会把函数误差压成大约

$$F(x_k)-F(x^*)\le\frac{C}{t_k^2}.$$

再代入$t_k\gtrsim k$，就得到

$$F(x_k)-F(x^*)=O(1/k^2).$$

经典固定步长结论可写成

$$F(x_k)-F(x^*)\le\frac{2L\lVert x_0-x^*\rVert_2^2}{(k+1)^2}.$$

因此证明主线其实只有：

$$\text{近端步基本不等式}\longrightarrow\text{构造加权势函数}\longrightarrow\text{用特殊$t_k$让项消掉}\longrightarrow t_k\sim k\longrightarrow O(1/k^2).$$

这条链比直接死背长证明更重要。

---

## 十一、$O(1/k^2)$到底快了多少

ISTA：

$$F(x_k)-F(x^*)\lesssim\frac Ck.$$

为了达到误差$\varepsilon$：

$$k=O(1/\varepsilon).$$

FISTA：

$$F(x_k)-F(x^*)\lesssim\frac C{k^2}.$$

于是

$$k=O(1/\sqrt\varepsilon).$$

比如只看数量级，如果希望误差从$1$压到$10^{-4}$：

- $1/k$量级大致需要$10^4$级别迭代；
- $1/k^2$量级大致需要$10^2$级别迭代。

实际常数、问题结构都会影响真实速度，但数量级上的差异已经说明为什么加速值得做。

更漂亮的是：FISTA每轮没有多算一次梯度或一次prox。

额外代价主要只有：

- 更新一个标量$t_k$；
- 做一次$x_k-x_{k-1}$；
- 做一次线性组合得到$y_{k+1}$。

所以它的核心价值是：

> 基本保留ISTA的单轮成本，却显著改善一般凸情形下的函数值复杂度。

---

## 十二、ISTA与FISTA放在一起看

| 对比 | ISTA | FISTA |
| --- | --- | --- |
| 基本问题 | $\min f+g$ | 同左 |
| 每轮梯度次数 | 1 | 1 |
| 每轮prox次数 | 1 | 1 |
| 近端梯度算子 | $p_L(\cdot)$ | 同一个$p_L(\cdot)$ |
| 近端步出发点 | $x_{k-1}$ | 外推点$y_k$ |
| 是否利用上一轮移动方向 | 否 | 是 |
| 额外状态 | 基本没有 | $t_k$与历史方向 |
| 函数值最坏收敛率 | $O(1/k)$ | $O(1/k^2)$ |
| 目标值是否必须逐轮下降 | 具有典型下降结构 | 标准FISTA不保证逐轮单调 |

最后一项很容易误解。

FISTA用了外推，因此某一轮可能“冲过头”，出现

$$F(x_{k+1})>F(x_k).$$

这并不自动意味着算法坏了，因为它的理论保证关注的是长期的函数值上界，而不是要求每一步都单调下降。

所以：

- `每轮下降`；
- `整体最坏情形收敛更快`；

是两件不同的事。

---

## 十三、把FISTA放回图像逆问题，一轮到底发生了什么

考虑

$$\min_x\frac12\lVert Ax-b\rVert_2^2+\lambda\lVert x\rVert_1.$$

第$k$轮先有一个外推点$y_k$。

第一步：预测观测

$$Ay_k.$$

第二步：算残差

$$Ay_k-b.$$

第三步：把误差反传回$x$空间

$$A^T(Ay_k-b).$$

第四步：做数据一致性梯度修正

$$z_k=y_k-\frac1L A^T(Ay_k-b).$$

第五步：做稀疏先验更新

$$x_k=S_{\lambda/L}(z_k).$$

第六步：利用最近运动趋势构造下一个外推点

$$y_{k+1}=x_k+\frac{t_k-1}{t_{k+1}}(x_k-x_{k-1}).$$

因此一轮FISTA可以压成四个角色：

$$\text{历史信息}\longrightarrow\text{数据一致性}\longrightarrow\text{先验约束}\longrightarrow\text{继续外推}. $$

如果以后进入Plug-and-Play，这里最值得记住的就是prox那一步：

$$\operatorname{prox}_{\gamma g}(z)$$

常常会被一个去噪器$D_\sigma(z)$替换。

于是

$$\text{显式正则的prox}\longrightarrow\text{隐式去噪先验}.$$

但要注意：普通去噪器未必真的是某个凸函数的prox，所以经典FISTA的收敛结论不能不加条件地直接搬过去。Zotero中的Sinha与Chaudhury 2024论文正是在特定线性去噪器和线性逆问题条件下进一步研究PnP-FISTA的迭代收敛。

---

## 十四、几个最容易混淆的点

### `FISTA是不是把ISTA的prox加速了？`

不是。

prox本身没有变化。FISTA主要改变的是`在哪个点做这一次近端梯度`：

$$x_{k-1}\quad\longrightarrow\quad y_k.$$

### `$O(1/k^2)$是不是说$x_k$到$x^*$的距离按$1/k^2$下降？`

经典FISTA结论首先说的是

$$F(x_k)-F(x^*)=O(1/k^2).$$

不是自动等价于

$$\lVert x_k-x^*\rVert=O(1/k^2).$$

### `为什么$g$不线性化？`

因为$g$可能不可微，而且往往恰好有容易计算的prox。完整保留它，反而比把它近似掉更有利。

### `为什么步长常写$1/L$？`

因为$L$控制$\nabla f$变化有多快，并给出二次上界；$L$越大，局部曲率可能越强，于是步长$1/L$越小。

### `FISTA的“惯性”是不是物理惯性？`

只是类比。真正的数学对象是外推

$$x_k+\beta_k(x_k-x_{k-1}),$$

而$\beta_k$不是随便选的，是为了满足加速证明的结构。

### `标准FISTA每一步都应该比前一步目标值更低吗？`

不要求。局部振荡可以出现，但仍可能满足整体$O(1/k^2)$复杂度。

---

## 十五、最后只保留一条真正的理解链

从最开始：

$$\min_x f(x)+g(x).$$

因为$f$光滑、$g$可能不可微，所以不能直接对整个目标做普通梯度下降。

先对$f$做一阶模型，并加入曲率安全垫：

$$f(y)+\langle\nabla f(y),x-y\rangle+\frac L2\lVert x-y\rVert_2^2.$$

再把$g(x)$完整放回来：

$$Q_L(x,y)=f(y)+\langle\nabla f(y),x-y\rangle+\frac L2\lVert x-y\rVert_2^2+g(x).$$

最小化这个局部模型得到

$$p_L(y)=\operatorname{prox}_{g/L}\left(y-\frac1L\nabla f(y)\right).$$

每次从上一个$x$做$p_L$，就是ISTA：

$$x_k=p_L(x_{k-1}),$$

函数值最坏收敛率为

$$O(1/k).$$

FISTA不改变$p_L$，而是先用历史方向构造

$$y_k=x_{k-1}+\beta_{k-1}(x_{k-1}-x_{k-2}),$$

再做

$$x_k=p_L(y_k).$$

特殊的$t_k$递推保证加权误差可以形成合适的势函数关系，最终得到

$$F(x_k)-F(x^*)=O(1/k^2).$$

所以FISTA真正该记成一句话：

> `FISTA = 近端梯度 + 精心设计的Nesterov型外推；它基本不改变每轮最贵的计算，却把一般凸复合问题的函数值收敛率从$O(1/k)$加速到$O(1/k^2)$。`

如果这句话里的每一部分都能解释出来，FISTA就已经不是三条需要硬背的公式，而是一套从梯度下降自然推出来的方法。

