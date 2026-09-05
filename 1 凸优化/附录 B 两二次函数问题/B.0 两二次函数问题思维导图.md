---
tags:
  - 暑假学习
  - 学习路径
  - 凸优化
  - 两二次函数问题
  - 思维导图
created: 2026-08-31
updated: 2026-08-31
---

# 附录 B 两二次函数问题思维导图

## 主线

- `隐藏凸性`：两个二次型的齐次联合值域是凸锥。
- `推理链`：联合值域秩一化 → S-过程无损 → 单约束二次强对偶 → SDP 松弛值精确。
- `关键限定`：两个二次函数、一个二次不等式、存在严格可行点。

## [[B.1 单约束二次优化]]

- 原问题：$\min q_0(x)$，满足$q_1(x)\le0$；$A_0,A_1$可不定。
- Lagrangian：$L=q_0+\lambda q_1$，$\lambda\ge0$。
- 有限下确界：$A_0+\lambda A_1\succeq0$且$b_0+\lambda b_1\in\mathcal R(A_0+\lambda A_1)$。
- 对偶 SDP：$\max\gamma$，满足$\lambda\ge0$及$\begin{bmatrix}A_0+\lambda A_1&b_0+\lambda b_1\\(b_0+\lambda b_1)^T&c_0+\lambda c_1-\gamma\end{bmatrix}\succeq0$。
- 严格可行$q_1(\hat x)<0$ → $p^\star=d^\star$。
- 提升：$X=xx^T$ → 放松为$X\succeq xx^T$ → $\begin{bmatrix}X&x\\x^T&1\end{bmatrix}\succeq0$。
- 严格可行 → SDP 松弛最优值精确；不保证每个最优$X$秩一。
- `信赖域`：$\lVert x\rVert_2^2\le\Delta^2$且$\Delta>0$自动严格可行。
- `边界`：等式球$\lVert x\rVert_2^2=\Delta^2$不由单不等式定理直接覆盖。

## [[B.2 S-过程]]

- 二次蕴含：$q_1(x)\le0\Rightarrow q_2(x)\le0$。
- 增广矩阵：$q_i(x)=(x,1)^TQ_i(x,1)$。
- 乘子证书：$Q_2\preceq\lambda Q_1$，$\lambda\ge0$。
- `充分性`：总成立。
- `必要性`：一个前件且存在$q_1(\hat x)<0$时成立。
- 椭球包含：$\mathcal E\subseteq\widetilde{\mathcal E}$ ↔ 单个乘子 LMI；连接[[基础知识/2026暑假学习/1 凸优化/8 几何问题/8.4 极值体积椭球]]。
- 多前件：$Q_0\preceq\sum_i\lambda_iQ_i$仍可靠，但一般有损。
- 严格/非严格：前件集合用$\le0$，正则条件必须有$<0$，乘子用$\ge0$。

## [[B.3 两个对称矩阵的值域]]

- 秩一化：对任意$X\succeq0$，存在$x$使$x^TAx=\operatorname{tr}(AX)$且$x^TBx=\operatorname{tr}(BX)$。
- 齐次值域：$W(A,B)=\{(x^TAx,x^TBx)\}$。
- 线性像：$W(A,B)=f(\mathbb S_+^n)$，故任意$n$下为凸锥。
- 归一化值域：$F(A,B)=\{(x^TAx,x^TBx)\mid\lVert x\rVert_2=1\}$。
- $F$在$n>2$时有凸性保证；$n=2$可能只是一条非凸圆周/椭圆周。
- `两个特殊`：两个二次测量可保留并秩一化；三个以上一般失败。
- 证明：高秩归纳降至秩二 → 在二维列空间显式构造$x=\alpha v_1+\beta v_2$。

## [[B.4 强对偶结果的证明]]

- LMI 择一：无乘子证书 ↔ 存在$X\succeq0$给出一严一弱的迹不等式。
- B.3：$X$ → 秩一向量$y=(v,w)$。
- $w\ne0$：取$x=v/w$。
- $w=0$：从严格可行点出发取$x(t)=\hat x+tv$，选大$\lvert t\rvert$与有利符号。
- 得到强择一 → S-过程。
- 下界蕴含：$q_1\le0\Rightarrow\gamma-q_0\le0$。
- S-过程把全部原问题下界精确变成对偶可行$\gamma$ → $p^\star=d^\star$。

## 边界速查

| 命题 | 结论 |
| --- | --- |
| 单前件乘子证书的充分性 | 无条件成立 |
| 单前件乘子证书的必要性 | 需要严格可行 |
| 多前件标量乘子法 | 通常只有充分性 |
| $W(A,B)$凸 | 任意维数 |
| $F(A,B)$凸 | $n>2$有保证，$n=2$可能失败 |
| 单约束非凸 QCQP 的值精确 SDP | 严格可行时成立 |
| SDP 最优矩阵必秩一 | 不成立 |
| 等式信赖域由 B.1 自动覆盖 | 不成立 |
