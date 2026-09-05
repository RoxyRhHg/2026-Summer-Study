---
tags:
  - 暑假学习
  - 学习路径
  - 凸优化
  - 强对偶证明
  - S过程
  - 分离证书
created: 2026-08-31
updated: 2026-08-31
---

# B.4 强对偶结果的证明

## 证明路线

- 附录中的两个强结果按以下链条建立：$$\text{LMI 择一定理}\longrightarrow\text{半正定证书 }X\longrightarrow\text{两个二次型的秩一化}\longrightarrow\text{S-过程}\longrightarrow\text{单约束二次强对偶}.$$
- [[B.3 两个对称矩阵的值域]]负责把半正定矩阵证书压缩为向量证书。
- 严格可行点负责排除 LMI 择一中的异常退化，并在齐次坐标为零时把“无穷远方向”转化为有限可行点。

## 从 LMI 择一到 S-过程

- 记$$Q_i=\begin{bmatrix}A_i&b_i\\b_i^T&c_i\end{bmatrix},\qquad i=1,2,$$并假设存在$\hat x$使$$\begin{bmatrix}\hat x\\1\end{bmatrix}^TQ_2\begin{bmatrix}\hat x\\1\end{bmatrix}<0.$$
- 因此$Q_2$至少有一个负特征值。于是$$\tau\ge0,\qquad\tau Q_2\succeq0\quad\Longrightarrow\quad\tau=0.$$
- 这个非退化条件允许应用[[5.8 例子|非严格 LMI 的择一定理]]：不存在$\lambda\ge0$使$$Q_1+\lambda Q_2\succeq0$$当且仅当存在$X\succeq0$使$$\operatorname{tr}(XQ_1)<0,qquad\operatorname{tr}(XQ_2)\le0.$$
- 这里一个迹不等式是严格的，另一个是非严格的；它们分别对应要找到的$r_1(x)<0$和$r_2(x)\le0$，不能互换或同时弱化。

## 用联合值域把矩阵证书秩一化

- 因证书只涉及$Q_1,Q_2$两个线性测量，B.3 的定理保证存在$y=(v,w)\in\mathbb R^{n+1}$使$$y^TQ_1y=\operatorname{tr}(XQ_1)<0,qquad y^TQ_2y=\operatorname{tr}(XQ_2)\le0.$$
- 现在需要把齐次向量$y=(v,w)$恢复为原变量$x\in\mathbb R^n$。

### 情形一：$w\ne0$

- 取$x=v/w$。由于二次齐次性，$$\begin{bmatrix}v\\w\end{bmatrix}=w\begin{bmatrix}x\\1\end{bmatrix},$$从而$$y^TQ_i y=w^2q_i(x).$$
- 因$w^2>0$，严格与非严格不等式方向都保持，所以$q_1(x)<0$且$q_2(x)\le0$。

### 情形二：$w=0$

- 此时证书只给出一个方向$v$：$$v^TA_1v<0,qquad v^TA_2v\le0.$$
- 不能令$x=v/w$。改从严格可行点出发，考察射线$$x(t)=\hat x+t v.$$
- 沿该射线有$$q_i(x(t))=q_i(\hat x)+2t(A_i\hat x+b_i)^Tv+t^2v^TA_i v.$$
- 因$v^TA_1v<0$，当$\lvert t\rvert$充分大时$q_1(x(t))<0$。
- 对$q_2$分三种情况：
  - 若$v^TA_2v<0$，二次项使$q_2(x(t))\to-\infty$，任一足够大的$\lvert t\rvert$都可；
  - 若$v^TA_2v=0$且$(A_2\hat x+b_2)^Tv\ne0$，选择$t$的符号使线性项趋向$-\infty$；
  - 若二次项和线性项都为零，则$q_2(x(t))=q_2(\hat x)<0$始终成立。
- 因而总能选到一个有限$t$，同时满足$q_1(x(t))<0$与$q_2(x(t))\le0$。
- 这完成了强择一：矩阵乘子证书不存在时，原二次不等式系统必有解；结合弱择一，就得到[[B.2 S-过程]]。

## 从 S-过程到单约束强对偶

- 对原问题$$p^\star=\inf\{q_0(x)\mid q_1(x)\le0\},$$一个标量$\gamma$是全局下界，当且仅当$$q_1(x)\le0\quad\Longrightarrow\quad q_0(x)\ge\gamma.$$
- 把结论写成$\gamma-q_0(x)\le0$，在$q_1$存在严格可行点时，S-过程说明上述蕴含等价于存在$\lambda\ge0$使$$\begin{bmatrix}-A_0&-b_0\\-b_0^T&\gamma-c_0\end{bmatrix}\preceq\lambda\begin{bmatrix}A_1&b_1\\b_1^T&c_1\end{bmatrix}.$$
- 移项得到$$\begin{bmatrix}A_0+\lambda A_1&b_0+\lambda b_1\\(b_0+\lambda b_1)^T&c_0+\lambda c_1-\gamma\end{bmatrix}\succeq0,$$这正是[[B.1 单约束二次优化|Lagrange 对偶 SDP]]的可行条件。
- 因此，“所有原问题下界$\gamma$”与“所有对偶可行目标值$\gamma$”完全相同。对两边取上确界即得$$d^\star=p^\star.$$
- 这个证明没有先假设$q_0$或$q_1$凸。强对偶来自下界蕴含能被一个非负乘子无损表示。

> [!warning] 教材排印边界
> 下界蕴含必须沿用原问题的$q_i(x)=x^TA_ix+2b_i^Tx+c_i$。教材 PDF 在 B.4 最后一处蕴含中漏排了系数$2$，但随后给出的增广矩阵仍对应$2b_i^Tx$；本笔记按全章一致的二次函数定义书写。

## 闭性、分离与达到性边界

- LMI 择一定理本质上是凸锥分离：若仿射射线$Q_1+\lambda Q_2$与半正定锥不相交，就用$X\succeq0$分离。
- 对一般凸集合，简单分离可能只得到非严格证书，或因像集不闭而无法得到所需证书。这里严格可行使$Q_2$含负方向，排除了$\tau Q_2\succeq0$的非零异常射线，从而得到所需的严格迹不等式。
- B.3 又保证两个迹测量的半正定像可以由秩一矩阵达到，而不仅是被秩一像的闭包逼近。因此证明需要的是`精确实现`，不是极限近似。
- 若$p^\star$有限，则$\gamma=p^\star$本身是原目标在可行集上的下界；S-过程因而给出达到$p^\star$的对偶乘子。原问题最优点是否达到是另一件事，不能只由“对偶达到”反推。
- 严格可行若失败，上述非退化步骤不再有保证；这并不等于强对偶必然失败，只表示当前证明与定理不能直接使用。

## 证明中最容易混淆的四点

- `半正定矩阵证书不是原变量`：必须经 B.3 秩一化，并处理齐次坐标$w$。
- `$w=0$不是无解`：它表示无穷远方向；严格可行点把该方向转成一族有限点$x(t)$。
- `强对偶不是凸性推论`：原问题可非凸，特殊性来自两个二次型的联合值域。
- `单约束不能推广为多约束`：三个以上二次测量时秩一化一般失败，标量乘子证书通常有损。
