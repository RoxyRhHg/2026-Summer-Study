---
type: 原子笔记
title: "[[线性可分时 Logistic 最大似然为什么可能不存在]]"
tags:
  - 暑假学习
  - 学习路径
  - 凸优化
  - 统计估计
  - Logistic回归
  - 线性可分
created: 2026-08-29
updated: 2026-08-30
---

## 完全线性可分

把标签改写成 $s_i=2y_i-1\in\{-1,1\}$，并使用增广特征 $\tilde u_i=(u_i^T,1)^T$。若存在向量 $\bar\theta$ 使

$$
s_i\tilde u_i^T\bar\theta>0,
\qquad i=1,\ldots,m,
$$

则数据完全线性可分：同一个超平面把两类样本严格分开。

## 沿分离方向放大参数

Logistic 单样本损失可写为

$$
\ell_i(\theta)=\log\left(1+e^{-s_i\tilde u_i^T\theta}\right).
$$

取 $\theta=t\bar\theta$，$t>0$，由于 $s_i\tilde u_i^T\bar\theta>0$，有

$$
\ell_i(t\bar\theta)
=\log\left(1+e^{-t s_i\tilde u_i^T\bar\theta}\right)
\longrightarrow0
\qquad(t\to\infty).
$$

所以总损失满足 $\inf_\theta\mathcal L(\theta)=0$。但对任意有限 $\theta$，每一项 $\ell_i(\theta)>0$，从而 $\mathcal L(\theta)>0$。

因此下确界 $0$ 只能在 $\|\theta\|\to\infty$ 时逼近，不能由任何有限参数取得。等价地，似然的上确界为 $1$，但没有有限的 ML 参数。

## 它说明什么

- 问题仍然是凸优化，因为 [[Logistic 损失为什么是凸函数|目标函数仍然凸]]。
- 失败的是最优解存在性，而不是凸性。
- 分类边界由 $\theta$ 的方向决定；把 $\theta$ 乘以正数基本不改变边界，却让预测概率越来越接近 $0$ 或 $1$。
- 对全部参数 $\theta=(a^T,b)^T$ 加入 $\lambda\|\theta\|_2^2$（$\lambda>0$），目标强凸且强制增长，因此存在唯一有限解。若只正则化 $a$ 而不正则化截距 $b$，还需排除沿截距方向发散的情形。

这就是 [[基础知识/2026暑假学习/1 凸优化/7 统计估计/7.1 参数分布估计#有限最优解不一定存在|Logistic 最大似然的边界条件]]。
