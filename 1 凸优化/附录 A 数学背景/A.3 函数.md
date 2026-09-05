---
tags:
  - 暑假学习
  - 学习路径
  - 凸优化
  - 数学背景
  - 函数
created: 2026-08-31
updated: 2026-08-31
---

# A.3 函数

## 函数记号与定义域

- 本教材写$f:A\to B$时，只说明输入和输出的环境空间；函数真正允许的输入组成$$\operatorname{dom}f\subseteq A,$$不要求$\operatorname{dom}f=A$。
- 例如$f:\mathbb S^n\to\mathbb R$、$f(X)=\log\det X$的真实定义域是$$\operatorname{dom}f=\mathbb S_{++}^n.$$半正定但奇异的矩阵满足$\det X=0$，不能作为有限值$\log\det X$的输入。
- 因此，写出函数公式时必须同时核对输入维度和定义域；完整符号约定见[[基础知识/2026暑假学习/1 凸优化/记号]]。

## 有限值表示与扩展值表示

- `有限值表示`只在$x\in\operatorname{dom}f$时讨论$f(x)\in\mathbb R$，定义域外没有函数值。
- `扩展值表示`把函数写成$f:\mathbb R^n\to\mathbb R\cup\{+\infty\}$，并令定义域外$f(x)=+\infty$。这样，约束$x\in C$可以并入目标：$$\min_{x\in C}g(x)=\min_x\left(g(x)+I_C(x)\right),$$其中$I_C(x)=0$当$x\in C$，否则$I_C(x)=+\infty$。
- 采用扩展值表示后，$+\infty$不是普通实数函数值，而是“此点不允许参与有限目标比较”的编码。

## 连续性

- 函数$f:\mathbb R^n\to\mathbb R^m$在$x\in\operatorname{dom}f$处连续，是指对每个$\epsilon>0$，存在$\delta>0$使$$y\in\operatorname{dom}f,\ \lVert y-x\rVert_2\le\delta\ \Longrightarrow\ \lVert f(y)-f(x)\rVert_2\le\epsilon.$$
- 注意$y$只在$\operatorname{dom}f$内趋近$x$；这是一种相对于定义域的连续性。
- 等价的序列判定是：若$x_k\in\operatorname{dom}f$且$x_k\to x\in\operatorname{dom}f$，则$$f(x_k)\to f(x).$$
- 若函数在定义域的每一点连续，则称其为连续函数。

## 闭函数与下半连续性

- 对实值函数$f$，若每个下水平集$$\{x\in\operatorname{dom}f\mid f(x)\le\alpha\},\;\alpha\in\mathbb R,$$都是闭集，则称$f$是`闭函数`。
- 等价地，函数的上图$$\operatorname{epi}f=\{(x,t)\in\mathbb R^{n+1}\mid x\in\operatorname{dom}f,\ f(x)\le t\}$$是闭集。
- 对扩展实值函数，这一性质也称为`下半连续`。序列形式是$$f(x)\le\liminf_{k\to\infty}f(x_k)\;\text{只要 }x_k\to x.$$
- `连续`与`闭`不是同一概念：连续性比较定义域内邻近点的函数值；闭性还控制趋近定义域边界时上图是否丢失极限点。
- 若$f$在闭定义域上连续，则$f$是闭函数。
- 若$\operatorname{dom}f$开且$f$连续，则$f$闭当且仅当：每个从定义域内部趋近任意有限边界点的序列都满足$f(x_k)\to+\infty$。

## 三个边界例子

- $f(x)=x\log x$、$\operatorname{dom}f=\mathbb R_{++}$不是闭函数，因为$x\downarrow0$时$f(x)\to0$，但极限点$0$未包含在定义域内。
- 把它闭扩展为$$f(x)=\begin{cases}x\log x,&x>0,\\0,&x=0,\end{cases}\;\operatorname{dom}f=\mathbb R_+,$$所得函数闭。
- $f(x)=-\log x$、$\operatorname{dom}f=\mathbb R_{++}$虽然定义域开放，但$x\downarrow0$时$f(x)\to+\infty$，所以它是闭函数。

## 闭性在最小化中的作用

- 若一列点$x_k$逼近$x$，闭函数不允许极限处的函数值突然高于邻近值的下极限。
- 因而闭性适合最小化：它防止最优值通过趋近一个未被函数图像正确包含的边界点而“丢失”。
- 闭性本身仍不保证最小值存在；通常还需要某个下水平集非空且紧，或函数具有足够的强制增长。
