---
type: 原子笔记
title: "[[Lipschitz 条件]]"
tags:
  - 暑假学习
  - 学习路径
  - 凸优化
  - FISTA
  - Lipschitz
  - 光滑性
  - 梯度下降
created: 2026-08-30
updated: 2026-08-30
aliases:
  - Lipschitz 连续
  - L-smooth
  - 梯度 Lipschitz
  - Lipschitz
---

# Lipschitz 条件

## 1. 它到底想限制什么？

`Lipschitz 条件`的核心思想不是“函数必须很平”，而是：

> 两个输入离得不远时，函数值或梯度也不能突然变化得离谱。

最基本的 Lipschitz 连续定义是：若存在常数$L>0$，使任意$x,y$满足

$$\lVert f(x)-f(y)\rVert\le L\lVert x-y\rVert,$$

则称$f$是$L$-Lipschitz 的。

直觉上，$L$控制函数的最大变化速度：

- $L$小：函数变化比较温和；
- $L$大：允许变化得更快；
- 若根本找不到有限$L$，说明函数可能在某些位置变化得非常剧烈。

---

## 2. 为什么在 FISTA 里更关心“梯度 Lipschitz”，而不是函数本身 Lipschitz？

在梯度下降、[[1 FISTA（讲解版）]] 中，更常见的条件是：

$$\lVert\nabla f(x)-\nabla f(y)\rVert_2\le L\lVert x-y\rVert_2.$$

这叫`梯度 Lipschitz`，也常称$f$是`L-smooth`。

它限制的不是$f$本身变化多快，而是：

> `梯度变化不能太快。`

因为梯度描述局部斜率，所以梯度变化速度本质上就在控制函数的`弯曲程度`。

因此在优化里可以把$L$先直观理解成：

$$L\approx\text{函数允许的最大曲率尺度}.$$

---

## 3. 一个一维例子

考虑

$$f(x)=\frac{a}{2}x^2.$$

梯度为

$$f'(x)=ax.$$

于是

$$\lvert f'(x)-f'(y)\rvert=\lvert a\rvert\lvert x-y\rvert.$$

因此梯度 Lipschitz 常数可以取

$$L=\lvert a\rvert.$$

如果$a$很大，抛物线很陡、弯曲得更厉害；如果$a$很小，函数更平缓。

这说明$L$确实能够反映函数曲率的大小。

---

## 4. 它和 Hessian 有什么关系？

如果$f$二阶可微，那么梯度变化由 Hessian 控制。

在一维中：

$$f'(x)-f'(y)=f''(\xi)(x-y)$$

其中$\xi$位于$x,y$之间。

所以若

$$\lvert f''(x)\rvert\le L,$$

就有

$$\lvert f'(x)-f'(y)\rvert\le L\lvert x-y\rvert.$$

多维情况下，若

$$\lVert\nabla^2 f(x)\rVert_2\le L,$$

则可推出

$$\lVert\nabla f(x)-\nabla f(y)\rVert_2\le L\lVert x-y\rVert_2.$$

对于凸且二阶可微的$f$，常写成

$$0\preceq\nabla^2 f(x)\preceq LI.$$

所以：

> `梯度 Lipschitz = Hessian 不会无限大 = 函数曲率有统一上界。`

---

## 5. 为什么这个条件会导出二次上界？

这是 Lipschitz 条件在优化里最关键的用途。

若$\nabla f$是$L$-Lipschitz，则有 `Descent Lemma`：

$$f(x)\le f(y)+\langle\nabla f(y),x-y\rangle+\frac{L}{2}\lVert x-y\rVert_2^2.$$

这意味着：

- $f(y)$：当前函数值；
- $\langle\nabla f(y),x-y\rangle$：一阶 Taylor 给出的线性预测；
- $\frac{L}{2}\lVert x-y\rVert_2^2$：为忽略的曲率提供安全补偿。

所以加入二次项并不是随意正则化，而是在用$L$给一阶 Taylor 误差罩一个可靠上界。

详细见 [[为什么在线性化模型中加入二次惩罚]]。

---

## 6. 为什么$L$越大，梯度下降步长反而越小？

最小化局部二次模型

$$f(y)+\langle\nabla f(y),x-y\rangle+\frac{L}{2}\lVert x-y\rVert_2^2$$

得到

$$x^+=y-\frac1L\nabla f(y).$$

因此典型步长为

$$\alpha=\frac1L.$$

直觉上：

- $L$大：函数可能弯得很厉害，当前梯度只能在较小邻域内相信，因此步长应小；
- $L$小：函数较平，局部线性近似在更大范围内仍较可靠，因此可以走得更远。

所以

$$L\uparrow\Rightarrow\text{曲率上界更大}\Rightarrow\text{更保守}\Rightarrow\alpha=1/L\downarrow.$$

---

## 7. 为什么它能够保证一步下降？

令

$$x^+=y-\frac1L\nabla f(y).$$

把它代入二次上界：

$$f(x^+)\le f(y)+\left\langle\nabla f(y),-\frac1L\nabla f(y)\right\rangle+\frac{L}{2}\left\lVert\frac1L\nabla f(y)\right\rVert_2^2.$$

整理：

$$f(x^+)\le f(y)-\frac1L\lVert\nabla f(y)\rVert_2^2+\frac1{2L}\lVert\nabla f(y)\rVert_2^2.$$

所以

$$f(x^+)\le f(y)-\frac1{2L}\lVert\nabla f(y)\rVert_2^2.$$

只要$\nabla f(y)\neq0$，右边严格小于$f(y)$。

因此 Lipschitz 梯度条件给出了一个很重要的保证：

> 步长取$1/L$时，梯度步不会因为走得过头而把函数值反而抬高。

---

## 8. 在线性最小二乘里，$L$怎么求？

考虑

$$f(x)=\frac12\lVert Ax-b\rVert_2^2.$$

梯度为

$$\nabla f(x)=A^T(Ax-b).$$

因此

$$\nabla f(x)-\nabla f(y)=A^TA(x-y).$$

利用谱范数：

$$\lVert A^TA(x-y)\rVert_2\le\lVert A^TA\rVert_2\lVert x-y\rVert_2.$$

所以可以取

$$L=\lVert A^TA\rVert_2=\lambda_{\max}(A^TA)=\lVert A\rVert_2^2.$$

如果目标写成

$$f(x)=\lVert Ax-b\rVert_2^2,$$

因为梯度多一个系数$2$，则

$$L=2\lVert A\rVert_2^2.$$

因此一定要注意目标函数前面有没有$\frac12$，否则 Lipschitz 常数会差一个$2$。

---

## 9. 不知道$L$怎么办？

实际问题中，精确的$L$可能难以计算。

这时可以使用`backtracking line search`：

1. 先猜一个$L$；
2. 用这个$L$做一步梯度或 proximal gradient；
3. 检查二次上界条件是否成立；
4. 若不成立，就增大$L$，也就是减小步长$1/L$；
5. 直到局部二次模型足够可靠。

所以 backtracking 的本质是：

> 自动寻找一个足够大的曲率上界$L$。

---

## 10. 最容易混淆：函数 Lipschitz 和梯度 Lipschitz 不是一回事

### 函数本身 Lipschitz

$$\lvert f(x)-f(y)\rvert\le L\lVert x-y\rVert.$$

它控制的是：`函数值变化速度`。

例如$f(x)=\lvert x\rvert$是 Lipschitz 连续的，但在$x=0$不可微。

### 梯度 Lipschitz

$$\lVert\nabla f(x)-\nabla f(y)\rVert\le L\lVert x-y\rVert.$$

它控制的是：`斜率变化速度 / 曲率`。

梯度法、ISTA、FISTA 中通常真正需要的是第二种。

因此看到“$f$是$L$-smooth”时，应自动翻译成：

$$\nabla f\text{ 是 }L\text{-Lipschitz}.$$

而不是$f$本身是$L$-Lipschitz。

---

## 11. 在 FISTA 中它到底负责什么？

FISTA 解决

$$\min_x F(x)=f(x)+g(x),$$

其中$f$光滑，$g$可以不可微。

对$f$使用

$$f(x)\le f(y)+\langle\nabla f(y),x-y\rangle+\frac L2\lVert x-y\rVert_2^2,$$

得到局部模型

$$Q_L(x,y)=f(y)+\langle\nabla f(y),x-y\rangle+\frac L2\lVert x-y\rVert_2^2+g(x).$$

然后最小化$Q_L$得到 proximal gradient：

$$p_L(y)=\operatorname{prox}_{g/L}\left(y-\frac1L\nabla f(y)\right).$$

所以 Lipschitz 条件在 FISTA 中的作用可以概括为：

> 它保证我们能够用“线性化$f$ + 二次项”构造一个可靠的局部模型，而$L$同时决定这个模型的曲率和梯度步长。

---

## 12. 最终应该怎样记？

不要只背：

$$\lVert\nabla f(x)-\nabla f(y)\rVert\le L\lVert x-y\rVert.$$

更应该形成下面这条因果链：

$$\text{梯度 Lipschitz}\Rightarrow\text{曲率有上界}\Rightarrow\text{一阶 Taylor 误差可被二次项控制}\Rightarrow\text{得到二次上界}\Rightarrow\text{安全步长约为 }1/L.$$

最核心的一句话是：

> `Lipschitz 梯度条件告诉我们：函数虽然会弯，但不会突然弯得无限厉害。`

因此我们才敢用当前梯度做局部预测，并用$\frac L2\lVert x-y\rVert^2$给这个预测加一个足够安全的曲率缓冲。 
