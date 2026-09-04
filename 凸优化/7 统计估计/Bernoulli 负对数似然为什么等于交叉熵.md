---
type: 原子笔记
title: "[[Bernoulli 负对数似然为什么等于交叉熵]]"
aliases:
  - Logistic 回归为什么产生交叉熵
  - Bernoulli 最大似然与二元交叉熵
tags:
  - 暑假学习
  - 学习路径
  - 凸优化
  - 统计估计
  - Logistic回归
  - 交叉熵
created: 2026-08-29
updated: 2026-08-30
---

## 结论成立需要哪些假设

对每个样本 $i$，假设

$$
Y_i\mid u_i\sim\operatorname{Bernoulli}(p_i),
\qquad
p_i=\sigma(a^Tu_i+b),
$$

并且给定解释变量后，各 $Y_i$ 条件独立。此时最大化联合似然等价于最小化二元交叉熵。

## 从两种结果合成一个公式

Bernoulli 分布满足

$$
\Pr(Y_i=1\mid u_i)=p_i,
\qquad
\Pr(Y_i=0\mid u_i)=1-p_i.
$$

由于 $y_i\in\{0,1\}$，两种情况可以合写为

$$
\Pr(Y_i=y_i\mid u_i)=p_i^{y_i}(1-p_i)^{1-y_i}.
$$

- $y_i=1$ 时，右侧为 $p_i$；
- $y_i=0$ 时，右侧为 $1-p_i$。

条件独立性给出

$$
L(a,b)=\prod_{i=1}^m p_i^{y_i}(1-p_i)^{1-y_i}.
$$

取负对数：

$$
\boxed{
-\log L(a,b)
=-\sum_{i=1}^m[y_i\log p_i+(1-y_i)\log(1-p_i)]
}.
$$

右侧就是未取平均的二元交叉熵。取总和还是样本平均值只相差正常数 $1/m$，不改变最优参数。

## 数值核对

若真实标签为 $y=1$，则 $\mathcal L=-\log p$：

- 预测 $p=0.9$ 时，$\mathcal L\approx0.105$；
- 预测 $p=0.1$ 时，$\mathcal L\approx2.303$。

若真实标签为 $y=0$，则 $\mathcal L=-\log(1-p)$，上面两个数值恰好交换。模型给真实类别的概率越接近零，损失越大；因此交叉熵会严厉处罚“高置信度的错误预测”。

## 不能混淆的两层结论

- Bernoulli 负对数似然等于二元交叉熵，与 $p_i$ 是否由 sigmoid 参数化无关；只要输出是合法的 Bernoulli 概率，这个等式就成立。
- sigmoid 规定 $p_i$ 如何依赖特征，并使 log-odds 关于特征仿射。
- 交叉熵关于参数是否凸，还要检查 $p_i$ 的参数化。Logistic 参数化下的凸性见 [[Logistic 损失为什么是凸函数]]。

因此 [[7.1 参数分布估计#Logistic 回归|Logistic 回归]] 中，“概率模型产生交叉熵”和“参数化保证凸性”是两个相邻但不同的结论。
