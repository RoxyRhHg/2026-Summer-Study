---
tags:
  - 暑假学习
  - 学习路径
  - 凸优化
  - 数学背景
  - 范数
created: 2026-08-31
updated: 2026-08-31
---

# A.1 范数

## 内积、欧几里得范数与角度

- 对$x,y\in\mathbb R^n$，标准内积为$$\langle x,y\rangle=x^Ty=\sum_{i=1}^n x_i y_i.$$本教材通常直接写$x^Ty$。
- `欧几里得范数`为$$\lVert x\rVert_2=(x^Tx)^{1/2}.$$其几何直觉与常用性质可回看[[基础知识/2026暑假学习/1 凸优化/欧几里得范数]]。
- `Cauchy–Schwarz 不等式`：$$\lvert x^Ty\rvert\le\lVert x\rVert_2\lVert y\rVert_2.$$它保证下式中的比值落在$[-1,1]$。
- 对非零向量$x,y$，无向夹角定义为$$\angle(x,y)=\cos^{-1}\left(\frac{x^Ty}{\lVert x\rVert_2\lVert y\rVert_2}\right)\in[0,\pi].$$当$x^Ty=0$时，称$x$与$y$正交。
- 对$X,Y\in\mathbb R^{m\times n}$，标准矩阵内积为$$\langle X,Y\rangle=\operatorname{tr}(X^TY)=\sum_{i=1}^m\sum_{j=1}^n X_{ij}Y_{ij}.$$
- 由该内积得到`Frobenius 范数`：$$\lVert X\rVert_F=\left(\operatorname{tr}(X^TX)\right)^{1/2}=\left(\sum_{i,j}X_{ij}^2\right)^{1/2}.$$它是把矩阵元素依次排成向量后得到的欧几里得范数。
- 若$X,Y\in\mathbb S^n$，则$X^T=X$，所以标准内积简化为$langle X,Y\rangle=\operatorname{tr}(XY)$。

## 一般范数、距离与单位球

- 函数$f:\mathbb R^n\to\mathbb R$在整个$\mathbb R^n$上有定义，并且满足以下四条时，称为`范数`，记作$f(x)=\lVert x\rVert$：
  - `非负性`：$\lVert x\rVert\ge0$；
  - `正定性`：$\lVert x\rVert=0$当且仅当$x=0$；
  - `绝对齐次性`：$\lVert tx\rVert=\lvert t\rvert\lVert x\rVert$；
  - `三角不等式`：$\lVert x+y\rVert\le\lVert x\rVert+\lVert y\rVert$。
- 范数诱导距离$$\operatorname{dist}(x,y)=\lVert x-y\rVert.$$因此“长度”和“距离”使用的是同一几何选择。不同范数的建模直觉见[[基础知识/2026暑假学习/1 凸优化/范数]]。
- `单位球`为$$\mathcal B=\{x\in\mathbb R^n\mid\lVert x\rVert\le1\}.$$它关于原点对称、凸、闭、有界，并且内部非空。
- 反过来，在有限维空间中，若集合$C$关于原点对称、凸、闭、有界且内部非空，则它是某个范数的单位球。对应范数是$C$的规范函数：$$\lVert x\rVert_C=\inf\{t>0\mid x\in tC\}.$$

## 常用范数

- 对$p\ge1$，`$\ell_p$范数`为$$\lVert x\rVert_p=\left(\sum_{i=1}^n\lvert x_i\rvert^p\right)^{1/p}.$$
- 两个端点和一个极限情形是$$\lVert x\rVert_1=\sum_i\lvert x_i\rvert,\;\lVert x\rVert_2=(x^Tx)^{1/2},\;\lVert x\rVert_\infty=\max_i\lvert x_i\rvert,$$并且$\lim_{p\to\infty}\lVert x\rVert_p=\lVert x\rVert_\infty$。
- 若$P\in\mathbb S_{++}^n$，则`二次范数`为$$\lVert x\rVert_P=(x^TPx)^{1/2}=\lVert P^{1/2}x\rVert_2.$$条件必须是$P\succ0$；若仅有$P\succeq0$且$P$奇异，则表达式可能在$x\ne0$时为零，只是半范数。
- 二次范数的单位球是椭球；反之，以原点为中心且内部非空的椭球可作为某个二次范数的单位球。
- 对矩阵$X\in\mathbb R^{m\times n}$，元素逐项范数还包括$$\lVert X\rVert_{\mathrm{sav}}=\sum_{i,j}\lvert X_{ij}\rvert,\;\lVert X\rVert_{\mathrm{mav}}=\max_{i,j}\lvert X_{ij}\rvert.$$

## 有限维中的范数等价

- 若$\lVert\cdot\rVert_a$和$\lVert\cdot\rVert_b$都是$\mathbb R^n$上的范数，则存在常数$\alpha,\beta>0$使$$\alpha\lVert x\rVert_a\le\lVert x\rVert_b\le\beta\lVert x\rVert_a,\;\forall x\in\mathbb R^n.$$
- 因此它们给出相同的开集、收敛序列和连续性概念，但常数$\alpha,\beta$会影响数值估计与算法界。
- `成立边界`：任意范数等价是有限维结论；在无限维函数空间中通常不能直接沿用。
- 更具体地，任意$\mathbb R^n$上的范数都存在一个二次范数$\lVert\cdot\rVert_P$满足$$\lVert x\rVert_P\le\lVert x\rVert\le\sqrt n\lVert x\rVert_P.$$这说明一般范数的单位球可由椭球在至多$\sqrt n$的因子内夹住。

## 算子范数

- 设$X\in\mathbb R^{m\times n}$，$\lVert\cdot\rVert_a$是$\mathbb R^m$上的范数，$\lVert\cdot\rVert_b$是$\mathbb R^n$上的范数。由它们诱导的`算子范数`为$$\lVert X\rVert_{a,b}=\sup_{\lVert u\rVert_b\le1}\lVert Xu\rVert_a=\sup_{u\ne0}\frac{\lVert Xu\rVert_a}{\lVert u\rVert_b}.$$
- 它度量线性映射从输入范数$\lVert\cdot\rVert_b$到输出范数$\lVert\cdot\rVert_a$的最大放大倍数，并满足$$\lVert Xu\rVert_a\le\lVert X\rVert_{a,b}\lVert u\rVert_b.$$
- 输入、输出都采用欧几里得范数时，$$\lVert X\rVert_2=\sigma_{\max}(X)=\sqrt{\lambda_{\max}(X^TX)}.$$这称为谱范数或$\ell_2$算子范数。
- 输入、输出都采用$\ell_\infty$范数时，诱导范数是最大行绝对值和：$$\lVert X\rVert_\infty=\max_i\sum_j\lvert X_{ij}\rvert.$$
- 输入、输出都采用$\ell_1$范数时，诱导范数是最大列绝对值和：$$\lVert X\rVert_1=\max_j\sum_i\lvert X_{ij}\rvert.$$
- `易混点`：$\lVert X\rVert_F$是元素向量化后的欧几里得范数；$\lVert X\rVert_2$是线性映射的最大放大倍数。两者一般不相等。

## 对偶范数

- 范数$\lVert\cdot\rVert$的`对偶范数`定义为$$\lVert z\rVert_*=\sup_{\lVert x\rVert\le1}z^Tx=\sup_{\lVert x\rVert\le1}\lvert z^Tx\rvert.$$
- 对所有$x,z$都有广义 Hölder 不等式$$z^Tx\le\lVert x\rVert\lVert z\rVert_*.$$对固定非零$x$或$z$，有限维中都能选到另一向量使等号成立。
- 在有限维空间中，对偶的对偶回到原范数：$$\lVert x\rVert_{**}=\lVert x\rVert.$$
- 若$1\le p\le\infty$且$1/p+1/q=1$，则$\ell_p$与$\ell_q$互为对偶；特别地，$\ell_2$自对偶，$\ell_1$与$\ell_\infty$互为对偶。
- 对矩阵使用迹内积$\operatorname{tr}(Z^TX)$时，谱范数的对偶是`核范数`：$$\lVert Z\rVert_*=\sup_{\lVert X\rVert_2\le1}\operatorname{tr}(Z^TX)=\sum_{i=1}^{\operatorname{rank}Z}\sigma_i(Z)=\operatorname{tr}\left((Z^TZ)^{1/2}\right).$$
- 对偶范数把“所有单位扰动下的最大线性响应”写成一个范数，因此会自然出现在鲁棒约束和等式约束范数最小化的对偶问题中。
