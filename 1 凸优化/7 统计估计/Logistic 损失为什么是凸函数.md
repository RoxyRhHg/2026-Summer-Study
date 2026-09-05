---
type: 原子笔记
title: "[[Logistic 损失为什么是凸函数]]"
tags:
  - 暑假学习
  - 学习路径
  - 凸优化
  - 统计估计
  - Logistic回归
  - 凸函数
created: 2026-08-29
updated: 2026-08-30
---

## 标量损失的二阶导

令增广特征与参数为

$$
\tilde u_i=\begin{bmatrix}u_i\\1\end{bmatrix},
\qquad
\theta=\begin{bmatrix}a\\b\end{bmatrix},
\qquad
z_i=\tilde u_i^T\theta.
$$

Bernoulli 负对数似然代入 $p_i=\sigma(z_i)$ 后为

$$
\ell_i(z_i)=\log(1+e^{z_i})-y_iz_i.
$$

一、二阶导数分别为

$$
\ell_i'(z_i)=\sigma(z_i)-y_i=p_i-y_i,
$$

$$
\ell_i''(z_i)=\sigma(z_i)[1-\sigma(z_i)]=p_i(1-p_i)>0.
$$

所以每个有限 $z_i$ 上的标量损失严格凸。

## 参数空间中的 Hessian

总损失为 $\mathcal L(\theta)=\sum_i\ell_i(\tilde u_i^T\theta)$，因此

$$
\nabla\mathcal L(\theta)=\sum_i(p_i-y_i)\tilde u_i,
$$

$$
\boxed{
\nabla^2\mathcal L(\theta)
=\sum_i p_i(1-p_i)\tilde u_i\tilde u_i^T
=\tilde U^TW\tilde U\succeq0
},
$$

其中 $\tilde U$ 的第 $i$ 行为 $\tilde u_i^T$，$W=\operatorname{diag}(p_i(1-p_i))\succ0$。

任取方向 $d$，

$$
d^T\nabla^2\mathcal L(\theta)d
=\sum_i p_i(1-p_i)(\tilde u_i^Td)^2\ge0,
$$

故 $\mathcal L$ 关于 $(a,b)$ 凸。

## 什么时候严格凸

因为有限参数下 $W\succ0$，

$$
d^T\nabla^2\mathcal Ld=0
\Longleftrightarrow
\tilde Ud=0.
$$

所以

$$
\boxed{
\mathcal L\text{ 关于参数严格凸}
\Longleftrightarrow
\tilde U\text{ 满列秩}
}
$$

（这里讨论有限参数处的 Hessian）。若增广特征矩阵不满列秩，则存在非零方向 $d$ 使所有 $z_i$ 不变，参数自然不能唯一识别。

严格凸至多保证“若最优解存在则唯一”，并不保证最优解一定存在。完全可分数据的不存在性见 [[线性可分时 Logistic 最大似然为什么可能不存在]]。

本结论为 [[基础知识/2026暑假学习/1 凸优化/7 统计估计/7.1 参数分布估计#凸性|Logistic 回归的凸性判断]] 提供了矩阵形式的证明。
