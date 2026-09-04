---
tags:
  - 暑假学习
  - 学习路径
  - 凸优化
  - 数学背景
  - 分析基础
created: 2026-08-31
updated: 2026-08-31
---

# A.2 分析基础

## 内点、开集与闭集

- 对$C\subseteq\mathbb R^n$，若$x\in C$且存在$\epsilon>0$使$$\{y\mid\lVert y-x\rVert_2\le\epsilon\}\subseteq C,$$则$x$是$C$的`内点`。
- 所有内点组成$\operatorname{int}C$。当$\operatorname{int}C=C$时，$C$是`开集`。
- `有限维条件`：由[[A.1 范数#有限维中的范数等价|范数等价性]]，把定义中的欧几里得球换成任意范数球，不会改变$\mathbb R^n$中的内点和开集。
- 若补集$\mathbb R^n\setminus C$是开集，则$C$是`闭集`。
- 低维仿射集合在环境空间中通常没有内点。例如，$\mathbb R^2$中的直线内部为空，但它仍可相对于自己的仿射包讨论相对内部；二者不能混用。

## 闭包、序列与边界

- `闭包`定义为$$\operatorname{cl}C=\mathbb R^n\setminus\operatorname{int}(\mathbb R^n\setminus C).$$
- 等价地，$x\in\operatorname{cl}C$当且仅当对每个$\epsilon>0$，都存在$y\in C$满足$\lVert x-y\rVert_2\le\epsilon$。
- 在$\mathbb R^n$中，$C$闭当且仅当：任何满足$x_k\in C$且$x_k\to x$的收敛序列，其极限$x$仍属于$C$。
- `边界`为$$\operatorname{bd}C=\operatorname{cl}C\setminus\operatorname{int}C.$$
- 若$x\in\operatorname{bd}C$，则任意小邻域同时包含$C$中的点和$C$外的点。
- 等价判定：$C$闭当且仅当$\operatorname{bd}C\subseteq C$；$C$开当且仅当$C\cap\operatorname{bd}C=\varnothing$。

## 上确界与达到

- 对$C\subseteq\mathbb R$，若$a$满足$x\le a$对所有$x\in C$成立，则$a$是$C$的`上界`。
- 若$C$有上界，则所有上界中最小者称为`上确界`，记作$\sup C$。
- 约定$$\sup\varnothing=-\infty,$$若$C$无上界，则$\sup C=+\infty$。
- `达到`与`存在上确界`不同：只有当$\sup C$有限且$\sup C\in C$时，才称上确界达到，并写作$$\max C=\sup C.$$
- 例如$C=(0,1)$满足$\sup C=1$，但$1\notin C$，所以$C$没有最大元。

## 下确界与达到

- 若$a\le x$对所有$x\in C$成立，则$a$是$C$的`下界`。
- `下确界`为所有下界中最大者，并可写成$$\inf C=-\sup(-C).$$
- 约定$$\inf\varnothing=+\infty,$$若$C$无下界，则$\inf C=-\infty$。
- 只有当$\inf C$有限且$\inf C\in C$时，才称下确界达到，并写作$$\min C=\inf C.$$
- 优化问题写成$\inf_x f(x)$只陈述最优值；写成$\min_x f(x)$还声称存在某个可行$x$取得该值。

## 为什么闭性与有界性影响最优解存在

- 在有限维空间中，闭且有界的集合是紧集。
- 若$C\subseteq\mathbb R^n$非空且紧，$f:C\to\mathbb R$连续，则$f$在$C$上同时取得最大值和最小值。
- 闭性防止极限点“逃出”可行集，有界性防止最优序列向无穷远逃逸，连续性保证函数值随极限传递。
- 这些条件是常用充分条件，不是必要条件；开放可行域上的问题仍可能有最优解，无界可行域上的强制函数也可能达到最小值。
