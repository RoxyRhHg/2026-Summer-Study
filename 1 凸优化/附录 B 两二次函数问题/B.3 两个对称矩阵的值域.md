---
tags:
  - 暑假学习
  - 学习路径
  - 凸优化
  - 联合值域
  - 二次型
  - 隐藏凸性
created: 2026-08-31
updated: 2026-08-31
---

# B.3 两个对称矩阵的值域

## 秩一表示定理

- 设$A,B\in\mathbb S^n$。对任意$X\in\mathbb S_+^n$，都存在$x\in\mathbb R^n$使$$x^TAx=\operatorname{tr}(AX),\qquad x^TBx=\operatorname{tr}(BX).$$
- 左边来自秩一矩阵$xx^T$，右边来自任意半正定矩阵$X$。定理说明：若只观察两个线性测量$\operatorname{tr}(AX)$和$\operatorname{tr}(BX)$，任意半正定矩阵都能由某个秩一半正定矩阵给出完全相同的测量值。
- 这不是说$X=xx^T$，而是说$X$与$xx^T$在由$A,B$定义的两个方向上不可区分。
- 该结果正是[[B.1 单约束二次优化|单约束二次 SDP 松弛]]可以无损的核心：目标和唯一约束恰好只需要保留两个二次测量。

## 齐次联合值域

- 定义$$W(A,B)=\{(x^TAx,x^TBx)\mid x\in\mathbb R^n\}\subseteq\mathbb R^2.$$
- $W(A,B)$是锥，因为若$w$由$x$产生，则$tw$可由$\sqrt t,x$产生，其中$t\ge0$。
- 定义线性映射$$f(X)=(\operatorname{tr}(AX),\operatorname{tr}(BX)).$$
- 每个$x$给出$X=xx^T\succeq0$，所以$W(A,B)\subseteq f(\mathbb S_+^n)$；秩一表示定理给出反向包含。因此$$W(A,B)=f(\mathbb S_+^n).$$
- $\mathbb S_+^n$是凸锥，线性映射的像仍是凸锥，所以$W(A,B)$对任意$n$都是凸锥。

## 归一化值域

- 定义$$F(A,B)=\{(x^TAx,x^TBx)\mid\lVert x\rVert_2=1\}.$$
- $F(A,B)$称为矩阵对$(A,B)$的二维值域，而$W(A,B)$是由它生成的锥：$$W(A,B)=\{t z\mid t\ge0,\ z\in F(A,B)\}.$$
- `关键区别`：$W(A,B)$对任意$n$都凸；$F(A,B)$在实数情形下由教材引用的结果保证$n>2$时凸，$n=1$时退化为单点，但$n=2$时可能不凸。
- 一个$n=2$反例是$$A=\begin{bmatrix}1&0\\0&-1\end{bmatrix},\qquad B=\begin{bmatrix}0&1\\1&0\end{bmatrix}.$$
- 对$x=(\cos\theta,\sin\theta)$，有$$(x^TAx,x^TBx)=(\cos2\theta,\sin2\theta),$$所以$F(A,B)$是单位圆周，不包含圆盘内部，因而不是凸集；但其锥包$W(A,B)=\mathbb R^2$仍然凸。

> [!warning] 条件不能移位
> “$n>2$”针对单位球面上的归一化值域$F(A,B)$；附录证明 S-过程使用的是齐次锥值域$W(A,B)$，后者对任意$n$凸。把前者的维数条件误加到后者，会错误缩小 S-过程的适用范围。

## 为什么只有两个二次型特殊

- 对两个矩阵，半正定锥在映射$f(X)=(\operatorname{tr}(AX),\operatorname{tr}(BX))$下的每个像点都能由秩一矩阵实现。
- 若改为三个矩阵$A,B,C$，一般不存在$x$同时满足$$x^TAx=\operatorname{tr}(AX),\quad x^TBx=\operatorname{tr}(BX),\quad x^TCx=\operatorname{tr}(CX).$$
- 例如在$\mathbb R^2$中，三个独立二次测量可以重构$X\in\mathbb S^2$；取$X=I$时不可能有$xx^T=I$，因为前者秩为$2$而后者秩至多为$1$。
- 因此三个及以上二次型的联合值域一般不凸，这解释了为什么单约束 QCQP 的 SDP 松弛可以精确，而多约束非凸 QCQP 通常只有下界。

## 构造性证明：先降到秩二

- 证明按$\operatorname{rank}X$归纳。
- 秩$0$时取$x=0$；秩$1$时若$X=xx^T$，结论直接成立。
- 对秩$k+1$的$X\succeq0$，可写成$$X=yy^T+Z,$$其中$Z\succeq0$且$\operatorname{rank}Z=k$。
- 由归纳假设，存在$z$使$$z^TAz=\operatorname{tr}(AZ),\qquad z^TBz=\operatorname{tr}(BZ).$$
- 因而$X$的两个测量值与秩至多为$2$的矩阵$yy^T+zz^T$相同。若已经证明秩二情形，就能把它进一步压缩为一个秩一矩阵。
- 所以全体问题归结为：怎样把一个秩二半正定矩阵在两个二次测量下压缩为秩一矩阵。

## 秩二情形

- 设$X=VV^T$，其中$V=[v_1\ v_2]\in\mathbb R^{n\times2}$列满秩。
- 在二维列空间内作正交变换，可令$$V^TAV=\begin{bmatrix}\lambda_1&0\\0&\lambda_2\end{bmatrix},\qquad V^TBV=\begin{bmatrix}\sigma_1&\gamma\\\gamma&\sigma_2\end{bmatrix}.$$
- 需要寻找$x=\alpha v_1+\beta v_2$，使其二次测量等于$$w=(\lambda_1+\lambda_2,\sigma_1+\sigma_2).$$

### 交叉项可由对角测量表示

- 若$(0,\gamma)$可写成$$ (0,\gamma)=z_1(\lambda_1,\sigma_1)+z_2(\lambda_2,\sigma_2),$$则只需选择$\alpha,\beta$满足$$\alpha^2+2\alpha\beta z_1=1,\qquad\beta^2+2\alpha\beta z_2=1.$$
- 令$t=\beta/\alpha$，第二个条件化为$$t^2+2t(z_2-z_1)=1.$$
- 该二次方程有一正一负两个实根；至少可选一个根使$1+2tz_1>0$，再取$$\alpha=\pm\frac1{\sqrt{1+2tz_1}},\qquad\beta=t\alpha.$$
- 这样构造的$x$同时保留$A$和$B$的两个测量值。

### 交叉项不在张成空间中

- 若$(0,\gamma)$不在$(\lambda_1,\sigma_1)$和$(\lambda_2,\sigma_2)$的张成空间中，则后两个向量必线性相关。
- 它们的和$w$是其中至少一个向量的非负倍数。若$w=\alpha^2(\lambda_1,\sigma_1)$，取$x=\alpha v_1$；另一种情形类似。
- 两种情形共同完成秩二压缩，再配合秩归纳即得一般结论。

## 非齐次二次函数与齐次化

- 非齐次二次函数$q_i(x)=x^TA_ix+2b_i^Tx+c_i$可通过$$Q_i=\begin{bmatrix}A_i&b_i\\b_i^T&c_i\end{bmatrix},\qquad y=\begin{bmatrix}x\\1\end{bmatrix}$$写成$q_i(x)=y^TQ_i y$。
- 但秩一表示定理得到的一般向量是$y=(v,w)$，不保证$w=1$。
- 若$w\ne0$，可以除以$w$恢复$x=v/w$；若$w=0$，无法直接反齐次化。[[B.4 强对偶结果的证明]]正是利用严格可行点沿方向$v$作大步移动，处理这一“无穷远方向”异常情形。

## 核心边界

- $W(A,B)$凸：任意维数，齐次、允许任意长度$x$。
- $F(A,B)$凸：实对称矩阵对在$n>2$时有保证；$n=2$可能失败。
- 两个测量可秩一化：成立。
- 三个及以上测量可秩一化：一般不成立。
- 齐次向量可直接使用：成立；非齐次问题还需处理最后一个齐次坐标是否为零。
