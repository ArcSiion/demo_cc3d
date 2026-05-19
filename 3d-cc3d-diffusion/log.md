# 3d-cc3d-diffusion Log

## 协作方式

- 后续围绕本文件推进：每次开始前先看 `log.md` 的“下一步”，完成实现后更新“已完成”和“验证结果”。
- “下一步”部分由助手先写草案，用户再修改确认。
- 如果用户修改了“下一步”，后续实现以用户修改后的版本为准。
- 每次迭代尽量保持目标小而清楚，避免一次性把 CC3D 特征、向量化、blocking 全部混在一起。

## 2026-05-19
### ====================================第一次迭代====================================
### 已完成

- 新建 `3d-cc3d-diffusion` 目录，作为 CompuCell3D PDE diffusion kernel 的补充实验入口。
- 从仓库根目录复制 `cc3d.md` 到本目录，作为该实验的源码依据、公式推导和论文表述边界说明。
- 创建最小可运行骨架：`Makefile`、`define.h`、`main.c`、`naive_scalar.c`。
- 初始阶段曾将 `3d-cc3d-diffusion` 加入 TemporalStencil 根目录 `Makefile`；后续为避免影响原项目，已迁移到 `/home/zhuhl/dataspace/demo_cc3d/` 并从原项目根 `Makefile` 移除。
- 在 `/home/zhuhl/dataspace/demo_cc3d/` 中创建独立根 `Makefile`，后续 demo 测试在该目录内执行。
- 从 TemporalStencil 中复制参考目录：`2d-jacobi-life`、`3d-jacobi-7p`、`3d-jacobi-27p`。
- 复制 TemporalStencil 公共构建依赖：`common.h`、`common.make`。
- 第一版只实现 `float` 精度的 `naive_scalar`，不引入 vector、vectime、blocking、cell type、boundary refresh 或 secretion。
- 当前公式采用 CC3D 常系数 diffusion + decay 的内部点特例：

```text
C_next = D * neighbor_sum + (1 - lambda - 6D) * center
```

- 当前参数为 `D = 0.1f`、`lambda = 0.0f`，因此 `center` 系数为 `0.4f`，与现有 `3d-jacobi-7p` 默认公式保持可对齐。

### 验证结果

- 在 TemporalStencil 原项目内，`make -C 3d-cc3d-diffusion clean && make -C 3d-cc3d-diffusion` 编译通过。
- 在 TemporalStencil 原项目内，`./exe-3d-cc3d-diffusion 16 16 16 4` 运行通过，输出：

```text
naive_scalar, NX = 16, NY = 16, NZ = 16, T = 4, checksum = 2.079554e+06, GStencil/s = 0.327680
```

- 在 demo 工作区执行 `make clean`，清理复制参考目录时带来的旧可执行文件。
- 在 demo 工作区执行 `make -C 3d-cc3d-diffusion && ./3d-cc3d-diffusion/exe-3d-cc3d-diffusion 16 16 16 4`，迁移后编译运行通过，输出：

```text
naive_scalar, NX = 16, NY = 16, NZ = 16, T = 4, checksum = 2.079554e+06, GStencil/s = 0.321255
```

### 下一步

- 先定实验分层，再继续写代码。当前目标是从 `naive_scalar` 出发，尽可能提取 CC3D 的真实特征，增加实验的生物背景；同时控制计算模式复杂度，避免后续 `vectime` 难以取得提升。

### 下一次迭代计划草案

#### 实验设计原则

- 保留 CC3D 的真实来源：每个版本都要能对应到 `cc3d.md` 中的源码依据或公式推导。
- 从标量版本开始：先把计算模式定义清楚，再做 `naive_vector`，最后再做 `vectime`。
- 特征逐步加入：每次只引入一个主要 CC3D 特征，方便解释性能变化。
- 不先做完整 solver：暂时不加入 secretion、uptake、variable diffusion coefficient、复杂 boundary refresh、BoxWatcher 和 CC3D interleaved layout。
- 后续能向量化：新增特征不能让数据依赖变成不规则图计算；尽量保持规则 3D stencil + 少量规则数组查表。

#### 建议版本序列

1. `V0_constant_scalar`
   - 当前已实现。
   - 常系数 `D` 和 `lambda`。
   - 公式为：

```text
C_next = D * neighbor_sum + (1 - lambda - 6D) * center
```

2. `V1_decay_scalar`
   - 显式把 `D`、`lambda` 作为实验参数保留下来。
   - 仍然是常系数 3D 7 点 stencil。
   - 目的：说明我们不是普通 3D7P，而是来自 CC3D diffusion + decay 公式。
   - 实现建议：先用宏控制 `CC3D_DIFFUSION_COEF` 和 `CC3D_DECAY_COEF`，暂时不做命令行参数。

3. `V2_cell_type_scalar`
   - 加入 CC3D 的核心生物建模特征：每个 lattice site 有 `cell_type`。
   - 增加：

```text
unsigned char cell_type[x][y][z]
float diff_coef[256]
float decay_coef[256]
```

   - 每个点按当前点类型查表：

```text
D_p = diff_coef[cell_type[p]]
lambda_p = decay_coef[cell_type[p]]
C_next = D_p * neighbor_sum + (1 - lambda_p - 6D_p) * center
```

   - 这是下一阶段最值得做的 CC3D 特征，因为它既有真实背景，又仍然是规则 stencil + 查表。

4. `V3_avoid_type_scalar`
   - 加入 non-diffusive / avoid type。
   - 邻居如果 `diff_coef[cell_type[q]] == 0`，则不加入 `neighbor_sum`，并减少有效邻居数 `N_p`。
   - 公式为：

```text
N_p = count(q in neighbors where diff_coef[cell_type[q]] != 0)
C_next = D_p * (sum_active_neighbors - N_p * center) + (1 - lambda_p) * center
```

   - 这个版本更接近 CC3D，但会引入分支或 mask。建议作为 scalar 生物背景版本先做，是否进入 `vectime` 视性能和实现复杂度决定。

5. `V4_naive_vector`
   - 在 `V0` 或 `V2` 稳定后再做。
   - 如果 `V2_cell_type_scalar` 的查表模式仍然清晰，优先给 `V2` 做 `naive_vector`；否则先给 `V0` 做。

6. `V5_vectime`
   - 不急着开始。
   - 等 `V0/V1/V2` 的标量语义和测试稳定后，再选择最合适的版本做时间向量化。

#### 下一次实际操作建议

- 暂时不要创建 `vectime`。
- 先把 `main.c` 和 `define.h` 整理成便于后续加版本的结构：
  - 保留当前 `naive_scalar`。
  - 增加一个独立 reference/check 机制，避免后续加入 cell type 后只有 checksum。
  - 在 `log.md` 记录固定测试命令和 checksum。
- 然后实现 `V1_decay_scalar` 或直接进入 `V2_cell_type_scalar`。如果优先增加 CC3D 背景，建议下一步直接做 `V2_cell_type_scalar`。
