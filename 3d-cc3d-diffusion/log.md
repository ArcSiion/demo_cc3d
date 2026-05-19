# 3d-cc3d-diffusion Log

## 记录规则

- 后续围绕本文件推进：每次开始前先看上一轮“下一步规划”，完成后追加新的更新块。
- 从第三次更新开始，只追加新内容，不修改历史更新块。
- 每次更新固定包含四大块：
  - `分析`
  - `做了什么`
  - `结果`
  - `下一步规划`
- “下一步规划”由助手先写草案，用户再修改确认。
- 如果用户修改了“下一步规划”，下一轮实现以用户修改后的版本为准。
- 每次迭代尽量保持目标小而清楚，避免一次性把 CC3D 特征、向量化、blocking 全部混在一起。
- 避免增加大块的代码或公式污染本文件。
- 每次按`下一步规划`执行更新后，在`log.md`中更新补充这是第i次更新

## 2026-05-19

### ====================================第一次更新====================================

#### 分析

- 目标是从 CompuCell3D CPU PDE solver 中抽取一个补充实验 kernel，证明工作不是凭空构造的 stencil。
- 第一阶段不直接做完整 CompuCell3D solver，也不马上做 `vectime`，而是先建立一个可运行、可解释、可逐步扩展的 `float` 标量 baseline。
- 为避免影响 TemporalStencil 原项目，后续实验需要迁移到独立 demo 工作区。
- 第一版采用 CC3D diffusion + decay 的内部点常系数特例：

```text
C_next = D * neighbor_sum + (1 - lambda - 6D) * center
```

- 当前参数为 `D = 0.1f`、`lambda = 0.0f`，因此 `center` 系数为 `0.4f`，与现有 `3d-jacobi-7p` 默认公式保持可对齐。
W
#### 做了什么

- 新建 `3d-cc3d-diffusion` 目录，作为 CompuCell3D PDE diffusion kernel 的补充实验入口。
- 从仓库根目录复制 `cc3d.md` 到本目录，作为该实验的源码依据、公式推导和论文表述边界说明。
- 创建最小可运行骨架：
  - `Makefile`
  - `define.h`
  - `main.c`
  - `naive_scalar.c`
- 第一版只实现 `float` 精度的 `naive_scalar`，不引入 vector、vectime、blocking、cell type、boundary refresh 或 secretion。
- 初始阶段曾将 `3d-cc3d-diffusion` 加入 TemporalStencil 根目录 `Makefile`；后续为避免影响原项目，已迁移到 `/home/zhuhl/dataspace/demo_cc3d/` 并从原项目根 `Makefile` 移除。
- 在 `/home/zhuhl/dataspace/demo_cc3d/` 中创建独立根 `Makefile`，后续 demo 测试在该目录内执行。
- 从 TemporalStencil 中复制参考目录：
  - `2d-jacobi-life`
  - `3d-jacobi-7p`
  - `3d-jacobi-27p`
- 复制 TemporalStencil 公共构建依赖：
  - `common.h`
  - `common.make`

#### 结果

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

#### 下一步规划

- 先定实验分层，再继续写代码。
- 当前目标是从 `naive_scalar` 出发，尽可能提取 CC3D 的真实特征，增加实验的生物背景。
- 同时控制计算模式复杂度，避免后续 `vectime` 难以取得提升。

### ====================================第二次更新====================================

#### 分析

- 从项目推进角度看，当前重点仍然不是马上写 `vectime`，而是先把实验分层定清楚，并保证每个版本都能对应到 CC3D 的真实背景。

#### 做了什么

- 将 `log.md` 整理为固定格式。
- 增加“记录规则”，明确从第三次更新开始只追加新内容，不修改历史更新块。
- 将第一次记录重排为四大块：
  - `分析`
  - `做了什么`
  - `结果`
  - `下一步规划`
- 将本次日志整理记为“第二次更新”。
- 保留第一次更新中的关键信息：
  - `float` 标量 baseline 已实现；
  - CC3D 常系数 diffusion 公式已落地；
  - demo 工作区已创建；
  - 参考目录已复制；
  - 迁移后编译运行通过。

#### 结果

- `log.md` 现在具备后续协作模板。
- 后续每轮只需要追加新的更新块，并按四大块填写即可。
- 当前代码没有变化，本次只整理日志。

#### 下一步规划

- 用户先检查并修改本节“下一步规划”。
- 如果用户确认继续按当前方向推进，第三次更新建议做以下事情：
  - 暂时不创建 `vectime`。
  - 验证机制可以先不做，执行exe-3d-cc3d-diffusion能返回GStencil/s等信息即可，因为我们要继续迭代不同的scalar版本。
  - 整理 `main.c` 和 `define.h`，让后续可以更容易添加不同 CC3D 特征版本。
  - 下一阶段优先考虑 `V2_cell_type_scalar`，因为它最能体现 CC3D 的生物建模背景，同时仍然保持规则 3D stencil + 查表模式。另外在做V2_cell_type_scalar时候，补充下我们做了什么，防止我不知道我们增加了什么模块。
- 建议版本序列如下：

    ```text
    V0_constant_scalar       已完成
    V1_decay_scalar          可作为参数语义补充
    V2_cell_type_scalar      优先建议下一步做
    V3_avoid_type_scalar     后续做，接近 CC3D 但会引入分支/maskW
    V4_naive_vector          在 scalar 语义稳定后做
    V5_vectime               最后选择合适版本做
    ```

### ====================================第三次更新====================================

#### 分析

- 本轮按照第二次更新的“下一步规划”执行。
- 暂时不创建 `vectime`，也不增加额外 reference/check 机制。
- 目标是在 scalar 阶段加入最有 CC3D 生物背景的特征：`cell_type + diff_coef/decay_coef` 查表。
- 该特征对应 CC3D 中每个 lattice site 有 cell type，并根据当前点类型选择 diffusion 和 decay 参数。
- 这个模式仍然保持规则 3D stencil，只是在每个点多了一次 cell type 查表和两次系数查表，后续仍有机会做 `naive_vector` 或 `vectime`。

#### 做了什么

- 整理 `define.h`，增加 CC3D cell type 相关定义：
  - `cc3d_cell_type_t`
  - `CC3D_TYPE_MEDIUM`
  - `CC3D_TYPE_CELL_A`
  - `CC3D_TYPE_CELL_B`
  - `CC3D_NUM_CELL_TYPES`
- 新增 `Compute_scalar_cell_type`，用于按当前点的 `cell_type` 查 `diffCoef` 和 `decayCoef`。
- 保留原有 `naive_scalar` 作为 `V0_constant_scalar`。
- 新增 `naive_scalar_cell_type.c`，实现 `V2_cell_type_scalar`。
- 在 `main.c` 中新增 cell type 初始化：
  - 中心区域设为 `CellA`；
  - 规则分布的一部分点设为 `CellB`；
  - 其余点设为 `Medium`。
- 在 `main.c` 中新增系数表初始化：
  - `Medium`: `D=0.10`, `decay=0.00`
  - `CellA`: `D=0.06`, `decay=0.02`
  - `CellB`: `D=0.12`, `decay=0.01`
- 程序现在会依次运行：
  - `naive_scalar`
  - `naive_scalar_cell_type`
- 程序会输出内部区域的 cell type 数量，便于确认 `V2` 确实启用了 cell type 分布。

#### 结果

- 编译命令通过：

```text
make clean && make
```

- 小规模运行命令通过：

```text
./exe-3d-cc3d-diffusion 16 16 16 4
```

- 输出结果：

```text
cell_type_counts, Medium = 2944, CellA = 512, CellB = 640
naive_scalar, NX = 16, NY = 16, NZ = 16, T = 4, checksum = 2.079554e+06, GStencil/s = 0.292571
naive_scalar_cell_type, NX = 16, NY = 16, NZ = 16, T = 4, checksum = 2.046478e+06, GStencil/s = 0.150312
```

- `naive_scalar_cell_type` 的 checksum 与 `naive_scalar` 不同是预期行为，因为不同 cell type 使用了不同 diffusion/decay 系数。
- `naive_scalar_cell_type` 的性能低于常系数版本也是预期行为，因为它增加了 `cell_type` 读取和系数查表。

#### 下一步规划

- 如果用户确认继续按当前方向推进，第四次更新建议做以下事情：
  - 暂时仍不创建 `vectime`。
  - 我先汇总一下我们接下来还没有模拟的功能，但不可能全都去做，为了temporal vectorization的效果好，我们要取舍：1. secretion / uptake 2. boundaryConditionInit 3. avoid type，也就是某些邻居不参与 diffusion 4. variableDiffusionCoefficientFlag，即邻居点 diffusion coefficient 也参与计算 5. cell type 随时间变化 6. OpenMP 并行 7. SIMD / temporal vectorization
  - 你接下来做：1.secretion / uptake我们不做，把扩散和衰减做好，我认为真实性就足够了。2.boundaryConditionInit边界条件初始化我们还是要做的，我认为这部分我们做无流出边界simple no-flux halo boundary就可以了，应该是保证化学物质不会从模拟区域边界流出去。这里我不太懂，你看着做吧，我们不做CC3D 里完整的 boundaryConditionInit，我希望对边界的处理我们有一个可解释的逻辑就行。3.当遇到不能参与diffusion的点时，我们不做跳过，而是对每个方向计算m_dir * (neighbor - center)，也就是：
                        new =
                            (1 - decay_center) * center
                          + D_center * Σ_dir [ m_dir * (neighbor_dir - center) ]，m_dir = 1  表示这个邻居参与扩散
  m_dir = 0  表示这个邻居是 avoid type，不参与扩散，也可以m_dir = 0.5   半透过，来体现不同 cell type 对扩散的阻挡程度不同。我们先做到上面这些，然后我们再看看效果。如果我的规划有错误、不合理的地方，写在第四次更新的分析部分

### ====================================第四次更新====================================

#### 分析

- 本轮按照第三次更新后用户修改的“下一步规划”执行。
- `secretion / uptake` 暂时不做是合理的：它们会在 diffusion 子步之间改写 field，容易破坏连续 diffusion 时间步的规则性，不利于后续 temporal vectorization。
- `boundaryConditionInit` 选择 simple no-flux halo boundary 是合理折中：它不是 CC3D 的完整边界系统，但语义清楚，表示化学物质不会从模拟区域边界流出。
- 对 avoid type 不采用“跳过邻居并减少邻居数”的写法，而采用方向 permeability `m_dir * (neighbor - center)`，这个抽象更适合后续向量化。
- 当前 permeability 语义为：`m=1.0` 正常扩散，`m=0.5` 半透，`m=0.0` 阻断。
- 这个模型是对 CC3D avoid/non-diffusive cell type 的可解释抽象，但不是 CC3D 原始代码逐行等价实现。它保留了规则 3D stencil 结构，同时加入了 cell type 对扩散通量的影响。

#### 做了什么

- 新增 `boundary.c`，实现 simple no-flux halo refresh。
- 在每个 scalar 版本的每个时间步开始前刷新当前时间层的 halo。
- 新增 cell type：
  - `CC3D_TYPE_BLOCKED`
  - `CC3D_TRACKED_CELL_TYPES`
- 扩展 cell type 初始化，当前内部点包含 `Medium`、`CellA`、`CellB`、`Blocked`。
- 新增 permeability 系数表：
  - `Medium`: `1.0`
  - `CellA`: `1.0`
  - `CellB`: `0.5`
  - `Blocked`: `0.0`
- 新增 `naive_scalar_permeability.c`，实现 `V3_permeability_scalar`。
- 当前程序依次运行：
  - `naive_scalar`
  - `naive_scalar_cell_type`
  - `naive_scalar_permeability`
- `naive_scalar_permeability` 使用当前点类型决定 `D_center` 和 `decay_center`，使用邻居点类型决定该方向 permeability。

#### 结果

- 编译命令通过：

```text
make clean && make
```

- 小规模运行命令通过：

```text
./exe-3d-cc3d-diffusion 16 16 16 4
```

- 输出结果：

```text
cell_type_counts, Medium = 2887, CellA = 512, CellB = 592, Blocked = 105
naive_scalar, NX = 16, NY = 16, NZ = 16, T = 4, checksum = 2.076110e+06, GStencil/s = 0.287439
naive_scalar_cell_type, NX = 16, NY = 16, NZ = 16, T = 4, checksum = 2.040889e+06, GStencil/s = 0.138847
naive_scalar_permeability, NX = 16, NY = 16, NZ = 16, T = 4, checksum = 2.040990e+06, GStencil/s = 0.055165
```

- `naive_scalar` 的 checksum 相比第三次更新发生变化是预期行为，因为现在每个时间步都会执行 no-flux halo refresh。
- `naive_scalar_permeability` 明显慢于 `naive_scalar_cell_type` 是预期行为，因为每个点额外读取 6 个邻居 cell type 并执行 6 个 permeability 乘法。

#### 下一步规划

- 用户先检查并修改本节“下一步规划”。
- 如果用户确认继续按当前方向推进，第五次更新建议做以下事情：
  - 暂时仍不创建 `vectime`。
  - 先做一轮 scalar 尺寸扫描，例如 `32^3`、`64^3`、`128^3`，观察三个版本的 GStencil/s 和性能差距。
  - 根据尺寸扫描结果决定后续向量化目标：
    - 如果 `naive_scalar_permeability` 过慢且复杂度过高，后续先对 `naive_scalar_cell_type` 做 `naive_vector`。
    - 如果 `naive_scalar_permeability` 性能仍可接受，再考虑给它做 vector 版本。
  - 暂时不加入 variable diffusion coefficient 和 cell type 随时间变化；这两个特征会进一步增加复杂度，可能不利于 temporal vectorization。
  - 建议下一步先补一个简单运行脚本，固定输入规模并记录三种 scalar 的输出，方便后续对比。
