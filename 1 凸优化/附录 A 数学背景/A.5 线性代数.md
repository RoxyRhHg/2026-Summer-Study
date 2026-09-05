---
tags:
  - 暑假学习
  - 学习路径
  - 凸优化
  - 数学背景
  - 线性代数
created: 2026-08-31
updated: 2026-08-31
---

# A.5 线性代数

## 值域、零空间与正交分解

- 对$A\in\mathbb R^{m\times n}$，`值域`是列向量的线性组合全体：$$\mathcal R(A)=\{Ax\mid x\in\mathbb R^n\}\subseteq\mathbb R^m.$$
- $\dim\mathcal R(A)=\operatorname{rank}A\le\min\{m,n\}$。当$\operatorname{rank}A=\min\{m,n\}$时，称$A$满秩；这不区分“满列秩”还是“满行秩”，需结合$m,n$判断。
- `零空间`为$$\mathcal N(A)=\{x\in\mathbb R^n\mid Ax=0\}.$$
- 子空间$\mathcal V\subseteq\mathbb R^n$的正交补为$$\mathcal V^\perp=\{x\mid z^Tx=0,\ \forall z\in\mathcal V\}.$$
- 基本关系是$$\mathcal N(A)=\mathcal R(A^T)^\perp,\;\mathcal R(A)=\mathcal N(A^T)^\perp,$$因此$$\mathcal N(A)\mathbin{\overset{\perp}{\oplus}}\mathcal R(A^T)=\mathbb R^n.$$
- 这意味着任意$x\in\mathbb R^n$都能唯一分解为零空间分量与行空间分量。

## 对称特征值分解

- 对$A\in\mathbb S^n$，存在正交矩阵$Q$和实对角矩阵$\Lambda$使$$A=Q\Lambda Q^T,\;Q^TQ=I,\;\Lambda=\operatorname{diag}(\lambda_1,\ldots,\lambda_n).$$
- 按$\lambda_1\ge\cdots\ge\lambda_n$排序，记$\lambda_{\max}(A)=\lambda_1$、$\lambda_{\min}(A)=\lambda_n$。
- 行列式、迹和两个常用范数可由特征值表示：$$\det A=\prod_i\lambda_i,\;\operatorname{tr}A=\sum_i\lambda_i,$$以及$$\lVert A\rVert_2=\max_i\lvert\lambda_i\rvert,\;\lVert A\rVert_F=\left(\sum_i\lambda_i^2\right)^{1/2}.$$
- `Rayleigh 商`给出极值特征值：$$\lambda_{\max}(A)=\sup_{x\ne0}\frac{x^TAx}{x^Tx},\;\lambda_{\min}(A)=\inf_{x\ne0}\frac{x^TAx}{x^Tx}.$$
- 因而对所有$x$，$$\lambda_{\min}(A)x^Tx\le x^TAx\le\lambda_{\max}(A)x^Tx.$$

## 正定、半正定与平方根

- 对称矩阵$A$满足$$A\succ0\Longleftrightarrow x^TAx>0\ \forall x\ne0\Longleftrightarrow\lambda_{\min}(A)>0.$$
- 半正定条件是$$A\succeq0\Longleftrightarrow x^TAx\ge0\ \forall x\Longleftrightarrow\lambda_i(A)\ge0\ \forall i.$$
- `条件差异`：$A\succ0$保证可逆；$A\succeq0$允许零特征值，因此可能奇异。
- 若$A\succeq0$，其唯一对称半正定平方根为$$A^{1/2}=Q\operatorname{diag}(\sqrt{\lambda_1},\ldots,\sqrt{\lambda_n})Q^T.$$只有$A\succ0$时，$A^{-1/2}$才存在。

## 广义特征值

- 对$A,B\in\mathbb S^n$，广义特征值由$$\det(sB-A)=0$$定义。
- 当且仅在这里采用的标准化前提$B\succ0$下，可构造对称矩阵$$\widetilde A=B^{-1/2}AB^{-1/2}.$$$(A,B)$的广义特征值就是$\widetilde A$的普通特征值，因此都是实数。
- 若$\widetilde A=Q\Lambda Q^T$，令$V=B^{1/2}Q$，则$$A=V\Lambda V^T,\;B=VV^T.$$这里$V$可逆但一般不正交。
- 极值广义特征值可写成$$\lambda_{\max}(A,B)=\sup_{x\ne0}\frac{x^TAx}{x^TBx},\;\lambda_{\min}(A,B)=\inf_{x\ne0}\frac{x^TAx}{x^TBx}.$$分母严格为正依赖$B\succ0$。

## 奇异值分解

- 若$A\in\mathbb R^{m\times n}$且$\operatorname{rank}A=r$，其紧奇异值分解为$$A=U\Sigma V^T,$$其中$U\in\mathbb R^{m\times r}$、$V\in\mathbb R^{n\times r}$满足$U^TU=V^TV=I_r$，且$$\Sigma=\operatorname{diag}(\sigma_1,\ldots,\sigma_r),\;\sigma_1\ge\cdots\ge\sigma_r>0.$$
- 等价的秩一展开为$$A=\sum_{i=1}^r\sigma_i u_i v_i^T.$$
- 非零奇异值的平方是$A^TA$和$AA^T$的非零特征值；$v_i$是右奇异向量，$u_i$是左奇异向量。
- 最大奇异值满足$$\sigma_{\max}(A)=\sup_{y\ne0}\frac{\lVert Ay\rVert_2}{\lVert y\rVert_2}=\lVert A\rVert_2.$$
- 教材把最小奇异值定义为$$\sigma_{\min}(A)=\begin{cases}\sigma_r(A),&r=\min\{m,n\},\\0,&r<\min\{m,n\}\end{cases}.$$因此$\sigma_{\min}(A)>0$当且仅当$A$满秩。
- 对可逆方阵$A\in\mathbb R^{n\times n}$，二范数条件数为$$\kappa(A)=\lVert A\rVert_2\lVert A^{-1}\rVert_2=\frac{\sigma_{\max}(A)}{\sigma_{\min}(A)}.$$

## Moore–Penrose 伪逆

- 由紧 SVD $A=U\Sigma V^T$定义$$A^\dagger=V\Sigma^{-1}U^T\in\mathbb R^{n\times m}.$$
- 两个极限表达式为$$A^\dagger=\lim_{\epsilon\downarrow0}(A^TA+\epsilon I)^{-1}A^T=\lim_{\epsilon\downarrow0}A^T(AA^T+\epsilon I)^{-1}.$$
- 特殊情形：
  - 满列秩时，$A^\dagger=(A^TA)^{-1}A^T$；
  - 满行秩时，$A^\dagger=A^T(AA^T)^{-1}$；
  - 方阵可逆时，$A^\dagger=A^{-1}$。
- 对任意$b$，$x=A^\dagger b$是$$\min_x\lVert Ax-b\rVert_2^2$$的最小欧几里得范数解。最小二乘解不唯一时，伪逆选择其中与$\mathcal N(A)$正交的一个。
- 两个投影矩阵是$$AA^\dagger=UU^T=\Pi_{\mathcal R(A)},\;A^\dagger A=VV^T=\Pi_{\mathcal R(A^T)}.$$
- 对一般二次函数$$\inf_x\left(\frac12x^TPx+q^Tx+r\right),\;P\in\mathbb S^n,$$最优值为$$p^\star=\begin{cases}-\frac12q^TP^\dagger q+r,&P\succeq0,\ q\in\mathcal R(P),\\-\infty,&\text{其他情形}.\end{cases}$$若$P$有负曲率，或$P\succeq0$但$q$在$\mathcal N(P)$方向仍有非零线性分量，目标都会向$-\infty$下降。

## 可逆主块的 Schur 补

- 设对称分块矩阵$$X=\begin{bmatrix}A&B\\B^T&C\end{bmatrix},\;A\in\mathbb S^k.$$若$A$可逆，则$A$在$X$中的`Schur 补`为$$S=C-B^TA^{-1}B.$$
- 行列式满足$$\det X=\det A\det S.$$因此在$A$可逆时，$X$可逆当且仅当$S$可逆。
- 当$A$和$S$都可逆时，$$X^{-1}=\begin{bmatrix}A^{-1}+A^{-1}BS^{-1}B^TA^{-1}&-A^{-1}BS^{-1}\\-S^{-1}B^TA^{-1}&S^{-1}\end{bmatrix}.$$
- 这个公式来自分块消元。实际计算时应通过解关于$A$和$S$的线性方程实现，而不是显式形成逆矩阵。

## Schur 补、部分最小化与矩阵正定性

- 当$A\succ0$且$v$固定时，关于$u$的二次函数$$u^TAu+2v^TB^Tu+v^TCv$$有唯一极小点$u^\star=-A^{-1}Bv$，最优值为$$v^T(C-B^TA^{-1}B)v=v^TSv.$$
- 因而$$X\succ0\Longleftrightarrow A\succ0\ \text{且}\ S\succ0.$$
- 在已经假设$A\succ0$时，还有$$X\succeq0\Longleftrightarrow S\succeq0.$$
- `易错条件`：上述半正定等价不能把$A\succ0$随意放宽成“$A$可逆”或“$A\succeq0$”；奇异情形需要额外值域条件。

## 奇异主块的广义 Schur 条件

- 若$A\succeq0$但可能奇异，则关于$u$的部分最小化只有在$Bv\in\mathcal R(A)$时有有限最优值，此时为$$v^T(C-B^TA^\dagger B)v.$$
- 若$Bv\notin\mathcal R(A)$，或者$A\not\succeq0$，该二次函数关于$u$无下界。
- 整个分块矩阵半正定的完整判定为$$X\succeq0\Longleftrightarrow A\succeq0,\ (I-AA^\dagger)B=0,\ C-B^TA^\dagger B\succeq0.$$
- $(I-AA^\dagger)B=0$等价于$\mathcal R(B)\subseteq\mathcal R(A)$。缺少这一条件时，即使$A\succeq0$和$C-B^TA^\dagger B\succeq0$，也不能推出$X\succeq0$。
