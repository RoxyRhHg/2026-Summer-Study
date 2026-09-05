---
tags:
  - 暑假学习
  - 学习路径
  - 凸优化
  - 数学背景
  - 思维导图
created: 2026-08-31
updated: 2026-08-31
---

# 附录 A 数学背景思维导图

## [[A.1 范数]]

### 内积几何

- 向量内积：$\langle x,y\rangle=x^Ty$；欧几里得范数：$\lVert x\rVert_2=(x^Tx)^{1/2}$。[[基础知识/2026暑假学习/1 凸优化/欧几里得范数]]
- Cauchy–Schwarz：$\lvert x^Ty\rvert\le\lVert x\rVert_2\lVert y\rVert_2$。
- 非零向量夹角：$\cos\angle(x,y)=x^Ty/(\lVert x\rVert_2\lVert y\rVert_2)$；$x^Ty=0\Leftrightarrow$正交。
- 矩阵内积：$\langle X,Y\rangle=\operatorname{tr}(X^TY)$；$\lVert X\rVert_F=(\operatorname{tr}(X^TX))^{1/2}$。

### 范数与单位球

- 四公理：非负、正定、绝对齐次、三角不等式。
- 距离：$\operatorname{dist}(x,y)=\lVert x-y\rVert$。
- 单位球：关于原点对称 + 凸 + 闭有界 + 内部非空。
- 反向构造：$\lVert x\rVert_C=\inf\{t>0\mid x\in tC\}$。[[基础知识/2026暑假学习/1 凸优化/范数]]
- $\ell_p$：$\lVert x\rVert_p=(\sum_i\lvert x_i\rvert^p)^{1/p}$，$p\ge1$；$p\to\infty$得到$\ell_\infty$。
- 二次范数：$\lVert x\rVert_P=(x^TPx)^{1/2}$必须有$P\succ0$。

### 等价、算子、对偶

- 有限维范数等价：$\alpha\lVert x\rVert_a\le\lVert x\rVert_b\le\beta\lVert x\rVert_a$；无限维不能直接沿用。
- 算子范数：$\lVert X\rVert_{a,b}=\sup_{u\ne0}\lVert Xu\rVert_a/\lVert u\rVert_b$。
- 谱范数：$\lVert X\rVert_2=\sigma_{\max}(X)$；$\lVert X\rVert_\infty=$最大行和；$\lVert X\rVert_1=$最大列和。
- 对偶范数：$\lVert z\rVert_*=\sup_{\lVert x\rVert\le1}z^Tx$；Hölder：$z^Tx\le\lVert z\rVert_*\lVert x\rVert$。
- 对偶对：$\ell_2\leftrightarrow\ell_2$，$\ell_1\leftrightarrow\ell_\infty$，$\ell_p\leftrightarrow\ell_q$且$1/p+1/q=1$。
- 矩阵：谱范数$\leftrightarrow$核范数$\sum_i\sigma_i$。

## [[A.2 分析基础]]

### 集合极限

- 内点：存在完全落在$C$内的小球；开集：$\operatorname{int}C=C$。
- 闭集：补集开；等价地，包含内部收敛序列的全部极限。
- 闭包：$x\in\operatorname{cl}C\Leftrightarrow$任意小邻域接触$C$。
- 边界：$\operatorname{bd}C=\operatorname{cl}C\setminus\operatorname{int}C$；任意小邻域同时接触集合内外。
- 有限维范数等价$\Rightarrow$开集、闭集、收敛概念与具体范数无关。

### 确界与达到

- $\sup C=$最小上界；$\inf C=$最大下界。
- $\sup\varnothing=-\infty$，$\inf\varnothing=+\infty$。
- $\sup C\in C\Rightarrow\max C=\sup C$；$\inf C\in C\Rightarrow\min C=\inf C$。
- 易混：有有限$\inf$不代表存在最小元。
- 有限维存在性入口：非空紧集 + 连续函数$\Rightarrow$最大、最小值均达到。

## [[A.3 函数]]

### 定义域

- $f:A\to B$只说明输入/输出环境；真实输入满足$\operatorname{dom}f\subseteq A$。[[基础知识/2026暑假学习/1 凸优化/记号]]
- 例：$f(X)=\log\det X$，输入类型$X\in\mathbb S^n$，真实定义域$\mathbb S_{++}^n$。
- 扩展值：定义域外置$+\infty$；$\min_{x\in C}g(x)=\min_x(g(x)+I_C(x))$。

### 连续与闭

- 连续：$x_k\to x$且$x_k,x\in\operatorname{dom}f\Rightarrow f(x_k)\to f(x)$。
- 闭函数：所有下水平集闭$\Leftrightarrow\operatorname{epi}f$闭$\Leftrightarrow$下半连续。
- 连续 + 闭定义域$\Rightarrow$闭函数。
- 连续 + 开定义域：闭函数需在每个有限边界点处趋于$+\infty$。
- $x\log x$在$(0,\infty)$上不闭；补上$f(0)=0$后闭；$-\log x$在$(0,\infty)$上闭。

## [[A.4 导数]]

### 一阶

- 可微点：$x\in\operatorname{int}(\operatorname{dom}f)$。
- $f:\mathbb R^n\to\mathbb R^m\Rightarrow Df(x)\in\mathbb R^{m\times n}$；$Df(x)_{ij}=\partial f_i/\partial x_j$。
- 标量函数：$\nabla f(x)=Df(x)^T\in\mathbb R^n$。
- 一阶近似：$f(x+d)\approx f(x)+Df(x)d$；标量情形为$f(x)+\nabla f(x)^Td$。
- 方向导数：$Df(x)[v]=\nabla f(x)^Tv$。
- 链式法则：$D(g\circ f)=Dg(f)Df$；若$h(x)=f(Ax+b)$，则$\nabla h(x)=A^T\nabla f(Ax+b)$。

### 二阶

- Hessian：$\nabla^2f=D(\nabla f)\in\mathbb R^{n\times n}$。
- 二阶近似：$f(x+d)\approx f(x)+\nabla f(x)^Td+\frac12d^T\nabla^2f(x)d$。
- 方向曲率：$v^T\nabla^2f(x)v$。
- 标量复合：$\nabla^2(g\circ f)=g'(f)\nabla^2f+g''(f)\nabla f\nabla f^T$。
- 仿射复合：$\nabla^2[f(Ax+b)]=A^T\nabla^2f(Ax+b)A$。
- $\log\det X$：$Df(X)[U]=\operatorname{tr}(X^{-1}U)$，$D^2f(X)[U,V]=-\operatorname{tr}(X^{-1}UX^{-1}V)$。

## [[A.5 线性代数]]

### 子空间

- 值域：$\mathcal R(A)=\{Ax\}$；零空间：$\mathcal N(A)=\{x\mid Ax=0\}$。
- $\mathcal N(A)=\mathcal R(A^T)^\perp$；$\mathcal N(A)\mathbin{\overset{\perp}{\oplus}}\mathcal R(A^T)=\mathbb R^n$。
- 满秩：$\operatorname{rank}A=\min\{m,n\}$；需结合形状区分满行秩/满列秩。

### 谱分解与正定性

- 对称分解：$A=Q\Lambda Q^T$，$Q^TQ=I$。
- Rayleigh 商：$\lambda_{\max}=\sup_{x\ne0}x^TAx/x^Tx$，$\lambda_{\min}=\inf_{x\ne0}x^TAx/x^Tx$。
- $A\succ0\Leftrightarrow\lambda_{\min}>0\Rightarrow A^{-1},A^{-1/2}$存在。
- $A\succeq0\Leftrightarrow\lambda_i\ge0$；允许零特征值和奇异性。
- 广义特征值：$B\succ0$时，$(A,B)$等价于$B^{-1/2}AB^{-1/2}$的普通特征值问题。

### SVD 与伪逆

- 紧 SVD：$A=U\Sigma V^T$；非零奇异值$\sigma_1\ge\cdots\ge\sigma_r>0$。
- $\sigma_i^2=$ $A^TA$与$AA^T$的非零特征值；$\sigma_{\max}=\lVert A\rVert_2$。
- 可逆方阵条件数：$\kappa(A)=\sigma_{\max}/\sigma_{\min}$。
- 伪逆：$A^\dagger=V\Sigma^{-1}U^T$。
- $A^\dagger b=$最小二乘问题的最小范数解。
- 投影：$AA^\dagger=\Pi_{\mathcal R(A)}$；$A^\dagger A=\Pi_{\mathcal R(A^T)}$。

### Schur 补

- $X=\begin{bmatrix}A&B\\B^T&C\end{bmatrix}$，$A$可逆：$S=C-B^TA^{-1}B$。
- $\det X=\det A\det S$；$A,S$均可逆时可写分块逆。
- $A\succ0$：$X\succ0\Leftrightarrow S\succ0$；$X\succeq0\Leftrightarrow S\succeq0$。
- $A\succeq0$可奇异：$X\succeq0\Leftrightarrow(I-AA^\dagger)B=0$且$C-B^TA^\dagger B\succeq0$。
- 值域条件：$(I-AA^\dagger)B=0\Leftrightarrow\mathcal R(B)\subseteq\mathcal R(A)$。
