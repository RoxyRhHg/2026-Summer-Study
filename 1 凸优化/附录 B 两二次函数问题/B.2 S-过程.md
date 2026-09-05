---
tags:
  - 暑假学习
  - 学习路径
  - 凸优化
  - S过程
  - 二次蕴含
  - LMI
created: 2026-08-31
updated: 2026-08-31
---

# B.2 S-过程

## 二次函数的矩阵表示

- 对$i=1,2$，记$$q_i(x)=x^TF_ix+2g_i^Tx+h_i,$$并定义增广矩阵$$Q_i=\begin{bmatrix}F_i&g_i\\g_i^T&h_i\end{bmatrix}.$$
- 令$y=(x,1)$，则$$q_i(x)=y^TQ_i y.$$
- `问题`：怎样用一个有限维矩阵不等式判断二次蕴含$$q_1(x)\le0\quad\Longrightarrow\quad q_2(x)\le0$$是否对所有$x\in\mathbb R^n$成立？

## S-过程的充分方向

- 若存在$\lambda\ge0$使$$Q_2\preceq\lambda Q_1,$$则对任意满足$q_1(x)\le0$的$x$，都有$$q_2(x)=y^TQ_2y\le\lambda y^TQ_1y=\lambda q_1(x)\le0.$$
- 因而，`存在非负乘子`总是二次蕴含的充分条件。
- 这个方向只用到矩阵不等式与$\lambda\ge0$，不要求严格可行，也不要求$q_1,q_2$凸。
- 对多个前件$q_i(x)\le0$，若存在$\lambda_i\ge0$使$$Q_0\preceq\sum_i\lambda_iQ_i,$$同样可推出所有前件蕴含$q_0(x)\le0$。这仍只是普遍成立的充分方向。

## 单一前件下的必要方向

- 假设前件存在严格可行点：$$\exists\hat x,\qquad q_1(\hat x)<0.$$
- 此时，若$q_1(x)\le0\Rightarrow q_2(x)\le0$对所有$x$成立，则存在$\lambda\ge0$使$$Q_2\preceq\lambda Q_1.$$
- 因此在`一个二次前件 + 严格可行`的条件下，二次蕴含与一个 LMI 完全等价。这一无损结论常称为`S-过程`或`S-lemma`。
- 严格可行排除了前件可行集只有退化边界的情形。若不存在$q_1(x)<0$，乘子证书仍可能存在，但必要性不再由该定理保证。
- 对本节的非齐次形式，定理借助齐次锥值域$W(A,B)$的凸性，而该凸性对任意维数成立；不要误把[[B.3 两个对称矩阵的值域|归一化值域]]要求的$n>2$直接附加到此处。

## 二次不等式的强择一

- 令$$r_i(x)=x^TA_ix+2b_i^Tx+c_i,$$并假设存在$\hat x$使$r_2(\hat x)<0$。
- 以下两个系统恰有一个可行：
  - 存在$x$使$r_1(x)<0$且$r_2(x)\le0$；
  - 存在$\lambda\ge0$使$$\begin{bmatrix}A_1&b_1\\b_1^T&c_1\end{bmatrix}+\lambda\begin{bmatrix}A_2&b_2\\b_2^T&c_2\end{bmatrix}\succeq0.$$
- 二者不可能同时成立，因为若同时成立，则增广向量$y=(x,1)$会给出$$0\le y^T(Q_1+\lambda Q_2)y=r_1(x)+\lambda r_2(x)<0,$$产生矛盾。
- 严格可行条件保证反方向也成立：若矩阵证书不存在，就一定能找到满足两个二次不等式的$x$。证明见[[B.4 强对偶结果的证明]]。

## 椭球包含

- 设两个有非空内部的椭球写成$$\mathcal E=\{x\mid q(x)\le0\},\qquad\widetilde{\mathcal E}=\{x\mid\widetilde q(x)\le0\},$$其中$$q(x)=x^TFx+2g^Tx+h,\qquad F\succ0,$$$$\widetilde q(x)=x^T\widetilde Fx+2\widetilde g^Tx+\widetilde h,\qquad\widetilde F\succ0.$$
- $\mathcal E$有非空内部等价于存在$q(x)<0$；例如其中心$-F^{-1}g$处的函数值是$h-g^TF^{-1}g$，所以要求$$h-g^TF^{-1}g<0.$$
- 包含关系$\mathcal E\subseteq\widetilde{\mathcal E}$就是$q(x)\le0\Rightarrow\widetilde q(x)\le0$。
- S-过程给出精确 LMI：存在$\lambda\ge0$使$$\begin{bmatrix}\widetilde F&\widetilde g\\\widetilde g^T&\widetilde h\end{bmatrix}\preceq\lambda\begin{bmatrix}F&g\\g^T&h\end{bmatrix}.$$
- 因$\widetilde F\succ0$，$\lambda=0$不可能满足上式，故在两个非退化椭球的场景中可进一步写成$\lambda>0$。
- 这一变换是[[基础知识/2026暑假学习/1 凸优化/8 几何问题/8.4 极值体积椭球]]中把“对椭球内所有点成立”的半无限约束化成 SDP 的关键工具。

## 多约束为何通常有损

- 若前件是$q_i(x)\le0$，$i=1,\ldots,m$，找到$\lambda_i\ge0$使$$Q_0\preceq\sum_{i=1}^m\lambda_iQ_i$$仍然保证蕴含成立。
- 但当$m\ge2$时，即使前件存在严格可行点，蕴含成立通常也不能保证存在这样的标量乘子。
- 原因不是简单的技术缺口：两个二次型的联合齐次值域具有特殊凸性，而三个及以上二次型的联合值域一般不再凸，分离证书可能无法压缩为一组标量乘子。
- 因而，多个不确定二次约束下的 S-过程常作为`保守充分条件`使用；若得到 LMI 可行，就获得可靠保证，若 LMI 不可行，不能据此断言原蕴含不成立。

## 严格与非严格符号

| 位置 | 形式 | 作用 |
| --- | --- | --- |
| 前件可行集 | $q_1(x)\le0$ | 定义需要被蕴含覆盖的集合 |
| 正则条件 | $q_1(\hat x)<0$ | 保证乘子条件的必要性 |
| 结论 | $q_2(x)\le0$ | 被验证的二次不等式 |
| 乘子 | $\lambda\ge0$ | 保持不等式方向 |

- 不能把正则条件中的$<$改成$\le$；后者只表示可行，不表示存在内部余量。
- 也不能把结论中的$\le0$随意换成$<0$。严格结论通常需要额外裕量，不能由同一个半正定证书直接等价表达。

## 使用边界

- `无条件成立`：乘子 LMI$\Rightarrow$二次蕴含。
- `有条件等价`：一个二次前件且存在严格可行点时，二次蕴含$\Leftrightarrow$乘子 LMI。
- `一般不等价`：多个二次前件下，标量乘子法通常只有充分性。
- `不要求凸性`：$F_1,F_2$都可以不定。
- `不要混淆`：S-过程的严格可行条件与“归一化联合值域在$n>2$时凸”是不同层次的条件。
