---
tags:
  - 暑假学习
  - 学习路径
  - 凸优化
  - 二次优化
  - 强对偶
  - SDP松弛
created: 2026-08-31
updated: 2026-08-31
---

# B.1 单约束二次优化

## 问题模型

- 考虑只有一条二次不等式约束的问题：$$\begin{aligned}\text{minimize}\quad&q_0(x)=x^TA_0x+2b_0^Tx+c_0\\\text{subject to}\quad&q_1(x)=x^TA_1x+2b_1^Tx+c_1\le0,\end{aligned}$$其中$x\in\mathbb R^n$，$A_i\in\mathbb S^n$，$b_i\in\mathbb R^n$，$c_i\in\mathbb R$。
- 这里不要求$A_0\succeq0$或$A_1\succeq0$。因此目标可能非凸，约束的次水平集也可能非凸；不能把本节结论当作普通凸 QCQP 的直接推论。
- 本节的特殊性来自`只含一条二次不等式`。一般多约束非凸 QCQP 不具有同样的精确强对偶结论。
- 与[[4.4 二次优化问题|凸 QCQP]]相比，本节放弃了半正定条件，却用单约束结构换回了隐藏的精确性。

## Lagrangian 与对偶函数

- 对偶变量为$\lambda\ge0$，Lagrangian 为$$L(x,\lambda)=q_0(x)+\lambda q_1(x)=x^TA(\lambda)x+2b(\lambda)^Tx+c(\lambda),$$其中$$A(\lambda)=A_0+\lambda A_1,\qquad b(\lambda)=b_0+\lambda b_1,\qquad c(\lambda)=c_0+\lambda c_1.$$
- 对称二次函数关于$x$有有限下确界，当且仅当$$A(\lambda)\succeq0,\qquad b(\lambda)\in\mathcal R(A(\lambda)).$$
- 在上述条件下，全部极小点满足$A(\lambda)x=-b(\lambda)$，对偶函数为$$g(\lambda)=c(\lambda)-b(\lambda)^TA(\lambda)^\dagger b(\lambda).$$
- 若半正定条件或值域条件任一失败，则$L(\cdot,\lambda)$沿某个方向无界下降，故$g(\lambda)=-\infty$。
- 这与[[5.1 Lagrange 对偶函数]]的通用原则一致：原问题即使非凸，对偶函数仍给出下界；有限性条件会变成对偶问题中的隐式约束。

## SDP 形式的 Lagrange 对偶

- 标量$\gamma$满足$\gamma\le g(\lambda)$，等价于二次函数$L(x,\lambda)-\gamma$对所有$x$非负。
- 利用广义 Schur 补，这一条件等价于$$\begin{bmatrix}A_0+\lambda A_1&b_0+\lambda b_1\\(b_0+\lambda b_1)^T&c_0+\lambda c_1-\gamma\end{bmatrix}\succeq0.$$
- 因而 Lagrange 对偶可写成 SDP：$$\begin{aligned}\text{maximize}\quad&\gamma\\\text{subject to}\quad&\lambda\ge0,\\&\begin{bmatrix}A_0+\lambda A_1&b_0+\lambda b_1\\(b_0+\lambda b_1)^T&c_0+\lambda c_1-\gamma\end{bmatrix}\succeq0.\end{aligned}$$
- 对偶只有两个标量变量$\gamma,\lambda$，但矩阵不等式的阶数为$n+1$。

## 非凸问题的强对偶

- `严格可行条件`：存在$\hat x$使$$q_1(\hat x)<0.$$
- 在该条件下，原问题与上面的 Lagrange 对偶具有相同最优值：$$p^\star=d^\star.$$
- 这个结论与凸问题的 Slater 定理外形相似，但逻辑来源不同：这里原问题可以非凸，精确性依赖两个二次型联合值域的特殊凸性以及[[B.2 S-过程|S-过程]]。
- 若$p^\star$有限，S-过程还会为最优下界$p^\star$给出某个$\lambda\ge0$，因此对偶最优值达到。强对偶本身不应被误读为“每个原最优点都容易恢复”或“任意多条二次约束也成立”。
- 严格可行是本节定理的充分正则条件；没有严格可行时不能直接引用该定理，即使某个具体退化问题仍可能恰好零对偶间隙。

## 半正定提升与松弛

- 利用$x^TA_ix=\operatorname{tr}(A_i xx^T)$，引入矩阵变量$X=xx^T$，原问题等价于$$\begin{aligned}\text{minimize}\quad&\operatorname{tr}(A_0X)+2b_0^Tx+c_0\\\text{subject to}\quad&\operatorname{tr}(A_1X)+2b_1^Tx+c_1\le0,\\&X=xx^T.\end{aligned}$$
- 非凸性集中在秩一等式$X=xx^T$。把它放松为$$X\succeq xx^T$$就扩大了可行域。
- 因右下角标量为$1>0$，Schur 补给出等价 LMI：$$X\succeq xx^T\quad\Longleftrightarrow\quad\begin{bmatrix}X&x\\x^T&1\end{bmatrix}\succeq0.$$
- 得到 SDP 松弛：$$\begin{aligned}\text{minimize}\quad&\operatorname{tr}(A_0X)+2b_0^Tx+c_0\\\text{subject to}\quad&\operatorname{tr}(A_1X)+2b_1^Tx+c_1\le0,\\&\begin{bmatrix}X&x\\x^T&1\end{bmatrix}\succeq0.\end{aligned}$$
- 因可行域被扩大，松弛最优值不大于原最优值。若某个 SDP 最优解满足$X=xx^T$，则该$x$一定是原问题的全局最优解。
- 更强的是：原问题严格可行时，此 SDP 松弛与原问题的`最优值精确相同`。但“最优值精确”不等于每个 SDP 最优矩阵都秩一，不能省略恢复原变量的讨论。

> [!warning] 教材排印边界
> 松弛只替换$X=xx^T$这一条约束，因此线性项必须继续是$2b_i^Tx$。教材 PDF 的式 (B.5) 排印成$b_i^Tx$，与式 (B.3)、(B.4)及“只放松秩约束”的叙述不一致；本笔记保留数学上一致的系数$2$。

## 信赖域问题

- 标准`信赖域子问题`是本节最重要的特例：$$\begin{aligned}\text{minimize}\quad&x^TAx+2b^Tx+c\\\text{subject to}\quad&\lVert x\rVert_2^2\le\Delta^2,\end{aligned}$$其中$A$可以不定，$\Delta>0$。
- 约束写成$q_1(x)=x^Tx-\Delta^2\le0$。取$x=0$即有$q_1(0)=-\Delta^2<0$，所以严格可行条件自动满足。
- 因而全局最优解$x^\star$可由下列条件刻画：存在$\lambda^\star\ge0$使$$A+\lambda^\star I\succeq0,$$$$$(A+\lambda^\star I)x^\star=-b,$$$$\lVert x^\star\rVert_2\le\Delta,\qquad\lambda^\star(\lVert x^\star\rVert_2^2-\Delta^2)=0.$$
- 这些条件不仅是局部一阶条件。半正定条件保证$x^\star$全局最小化$L(\cdot,\lambda^\star)$，互补条件再把该下界与原目标在$x^\star$处对齐，因此它们对全局最优具有充分性。
- 若无约束二次函数已有内部极小点，则可有$\lambda^\star=0$；若最优点被信赖域边界截住，则通常$\lambda^\star>0$且$\lVert x^\star\rVert_2=\Delta$。

> [!warning] 不等式球与等式球
> 本节定理直接覆盖$\lVert x\rVert_2^2\le\Delta^2$，不直接覆盖$\lVert x\rVert_2^2=\Delta^2$。把等式拆成两个相反方向的不等式会产生两条二次约束，已经超出“单不等式”结论；即使等式信赖域另有隐藏凸性结果，也必须单独验证，不能由本节定理自动推出。

## 使用边界

- `可以推出`：单二次不等式、严格可行时，Lagrange 对偶与 SDP 松弛在最优值上精确。
- `不能推出`：多条非凸二次约束的一般 QCQP 也零对偶间隙。
- `不能推出`：严格可行是强对偶的必要条件。
- `不能推出`：SDP 最优解必定直接给出秩一矩阵。
- `不能忽略`：二次函数下确界有限时除了$A(\lambda)\succeq0$，还需要$b(\lambda)\in\mathcal R(A(\lambda))$。
