---
type: 原子笔记
title: "[[最大熵与最小 KL 散度的关系]]"
aliases:
  - 最大熵是均匀先验下的最小KL投影
tags:
  - 暑假学习
  - 学习路径
  - 凸优化
  - 统计估计
  - 最大熵
  - KL散度
created: 2026-08-29
updated: 2026-08-30
---

## 两个问题看似不同

设 $\mathcal P$ 是满足已知矩、概率或其他凸约束的分布集合。

最大熵在 $\mathcal P$ 中选择

$$
\max_{p\in\mathcal P}H(p),
\qquad
H(p)=-\sum_i p_i\log p_i.
$$

它适合“只有约束，没有偏好的参考分布”的情形：在已知信息之外，不额外把概率质量集中到少数结果。

最小 KL 散度则从参考分布 $q$ 出发：

$$
\min_{p\in\mathcal P}D_{\mathrm{KL}}(p\|q),
\qquad
D_{\mathrm{KL}}(p\|q)=\sum_i p_i\log\frac{p_i}{q_i}.
$$

它适合“已有一个基准分布，但新约束要求修正”的情形：在满足新信息的前提下，尽量少偏离 $q$。

## 均匀参考分布把二者统一起来

若 $q_i=1/n$，则

$$
\begin{aligned}
D_{\mathrm{KL}}(p\|q)
&=\sum_i p_i\log\frac{p_i}{1/n}\\
&=\sum_i p_i\log p_i+\log n\sum_i p_i\\
&=\sum_i p_i\log p_i+\log n.
\end{aligned}
$$

因为 $\log n$ 与 $p$ 无关，

$$
\arg\min_{p\in\mathcal P}D_{\mathrm{KL}}(p\|q)
=\arg\min_{p\in\mathcal P}\sum_i p_i\log p_i
=\arg\max_{p\in\mathcal P}H(p).
$$

所以最大熵并非独立于最小 KL 的另一套方法，而是==参考分布为均匀分布时的最小 KL 投影==。

## 选择时问什么

- 只有约束，不希望引入额外偏好：使用最大熵；
- 已有可信参考分布 $q$，希望最小幅度更新：使用最小 KL；
- 有观测频数并相信 IID 采样：优先从最大似然出发。

三者都在选择概率向量，但依据不同：最大似然服从数据，最大熵避免无依据的集中，最小 KL 保留参考分布。它们在 [[7.2 非参数分布估计]] 中共享同一个凸可行集 $\mathcal P$，却回答不同的信息融合问题。

