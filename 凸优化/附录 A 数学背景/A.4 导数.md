---
tags:
  - 暑假学习
  - 学习路径
  - 凸优化
  - 数学背景
  - 导数
created: 2026-08-31
updated: 2026-08-31
---

# A.4 导数

## 导数是最佳一阶线性近似

- 设$f:\mathbb R^n\to\mathbb R^m$，并取$x\in\operatorname{int}(\operatorname{dom}f)$。若存在矩阵$Df(x)\in\mathbb R^{m\times n}$使$$\lim_{z\to x,\ z\in\operatorname{dom}f,\ z\ne x}\frac{\lVert f(z)-f(x)-Df(x)(z-x)\rVert_2}{\lVert z-x\rVert_2}=0,$$则$f$在$x$处可微，$Df(x)$称为导数或`Jacobian`。
- 这里要求$x$位于定义域内部，使各个足够小方向都可用于局部比较。边界点可能有单侧方向导数，但不自动具有上述全空间导数。
- 一阶近似为$$f(z)\approx f(x)+Df(x)(z-x).$$误差相对于$\lVert z-x\rVert_2$是高阶小量。
- Jacobian 的维度必须与线性映射一致：$$Df(x)\in\mathbb R^{m\times n},\;Df(x)_{ij}=\frac{\partial f_i(x)}{\partial x_j}.$$每一行是一个输出分量对全部输入变量的导数。

## 梯度与方向导数

- 当$f:\mathbb R^n\to\mathbb R$为标量函数时，$Df(x)$是$1\times n$行向量，`梯度`定义为列向量$$\nabla f(x)=Df(x)^T\in\mathbb R^n.$$
- 标量函数的一阶近似是$$f(z)\approx f(x)+\nabla f(x)^T(z-x).$$
- 沿方向$v\in\mathbb R^n$限制到直线$\widetilde f(t)=f(x+tv)$，得到方向导数$$\widetilde f'(t)=\nabla f(x+tv)^Tv,$$特别地，$D f(x)[v]=\nabla f(x)^Tv$。

## 两个一阶例子

- 对$P\in\mathbb S^n$，二次函数$$f(x)=\frac12x^TPx+q^Tx+r$$的梯度为$$\nabla f(x)=Px+q.$$
- 对矩阵变量$X\in\mathbb S_{++}^n$，$f(X)=\log\det X$的方向微分为$$Df(X)[U]=\operatorname{tr}(X^{-1}U).$$以迹内积$\langle G,U\rangle=\operatorname{tr}(GU)$识别梯度，得到$$\nabla f(X)=X^{-1}.$$
- `维度提醒`：矩阵变量上的梯度仍属于同一个矩阵空间$\mathbb S^n$；不能不加说明地把它当成$n^2$维列向量。

## 一阶链式法则

- 设$f:\mathbb R^n\to\mathbb R^m$在$x$可微，$g:\mathbb R^m\to\mathbb R^p$在$f(x)$可微，并令$h=g\circ f$，则$$Dh(x)=Dg(f(x))Df(x).$$
- 维度核对为$$(p\times n)=(p\times m)(m\times n).$$矩阵乘法顺序由复合顺序决定，不能交换。
- 当$f:\mathbb R^n\to\mathbb R$、$g:\mathbb R\to\mathbb R$时，$$\nabla(g\circ f)(x)=g'(f(x))\nabla f(x).$$
- 若$g(x)=f(Ax+b)$，其中$A\in\mathbb R^{n\times p}$，则$$Dg(x)=Df(Ax+b)A.$$若$f$为标量函数，则$$\nabla g(x)=A^T\nabla f(Ax+b).$$
- 对$$F(x)=F_0+\sum_{i=1}^n x_iF_i\succ0,\;f(x)=\log\det F(x),$$有$$\frac{\partial f(x)}{\partial x_i}=\operatorname{tr}(F(x)^{-1}F_i).$$这里矩阵仿射映射与$\log\det$的链式法则共同决定结果。

## Hessian 与二阶近似

- 对二次可微的标量函数$f:\mathbb R^n\to\mathbb R$，`Hessian`为$$\nabla^2f(x)_{ij}=\frac{\partial^2 f(x)}{\partial x_i\partial x_j}\in\mathbb R^{n\times n}.$$
- 二阶近似为$$\widehat f(z)=f(x)+\nabla f(x)^T(z-x)+\frac12(z-x)^T\nabla^2f(x)(z-x).$$并满足$f(z)-\widehat f(z)=o(\lVert z-x\rVert_2^2)$。
- Hessian 是梯度映射的 Jacobian：$$D(\nabla f)(x)=\nabla^2f(x).$$在二阶偏导连续时，Hessian 对称。
- 对前述二次函数，$\nabla^2f(x)=P$，所以二阶近似就是函数本身。
- 对$f(X)=\log\det X$，二阶双线性型为$$D^2f(X)[U,V]=-\operatorname{tr}(X^{-1}UX^{-1}V).$$特别地，$D^2f(X)[U,U]\le0$，这揭示$\log\det$在$\mathbb S_{++}^n$上的凹性。

## 二阶链式法则

- 若$h(x)=g(f(x))$，其中$f:\mathbb R^n\to\mathbb R$、$g:\mathbb R\to\mathbb R$二次可微，则$$\nabla^2h(x)=g'(f(x))\nabla^2f(x)+g''(f(x))\nabla f(x)\nabla f(x)^T.$$
- 第一项继承内层函数曲率，第二项来自外层函数曲率；遗漏任一项都会得到错误 Hessian。
- 若$g(x)=f(Ax+b)$且$f$为标量函数，则$$\nabla^2g(x)=A^T\nabla^2f(Ax+b)A.$$
- 在线方向$\widetilde f(t)=f(x+tv)$上，$$\widetilde f''(t)=v^T\nabla^2f(x+tv)v.$$这把多维曲率压缩为方向$v$上的一维曲率。

## Log-sum-exp 的导数结构

- 设$$f(x)=\log\sum_{i=1}^m\exp(a_i^Tx+b_i),$$令$A\in\mathbb R^{m\times n}$的第$i$行为$a_i^T$，$z_i=\exp(a_i^Tx+b_i)$，并定义$p=z/(\mathbf1^Tz)$。
- 梯度为$$\nabla f(x)=A^Tp.$$
- Hessian 为$$\nabla^2f(x)=A^T\left(\operatorname{diag}(p)-pp^T\right)A.$$
- 中间矩阵$\operatorname{diag}(p)-pp^T$是概率向量$p$对应的协方差矩阵，因而半正定；这直接给出 log-sum-exp 的凸性。
