# CompuCell3D PDE Solver 抽取背景

本文档记录从 CompuCell3D CPU PDE solver 中抽取补充实验 kernel 的依据。目标不是声称已经加速完整 CompuCell3D 系统，而是从真实应用代码中抽取一个有明确来源的 diffusion stencil，用于验证时间维度向量化在应用型 stencil 上的效果。

## 1. 源码背景

候选源码来自当前本地 CompuCell3D checkout：

- `/home/zhuhl/dataspace/cc3d_study/CompuCell3D/CompuCell3D/core/CompuCell3D/steppables/PDESolvers/DiffusionSolverFE_CPU.cpp`
- `/home/zhuhl/dataspace/cc3d_study/CompuCell3D/CompuCell3D/core/CompuCell3D/Field3D/Array3D.h`

关键函数是 `DiffusionSolverFE_CPU::diffuseSingleField`。它在 `stepImpl` 中被反复调用：每个 field 会根据 `scalingExtraMCSVec[i]` 执行若干 PDE diffusion 子步，每个子步先刷新边界条件，再执行一次 diffusion。这提供了时间维度 `T`，因此适合抽取为时间向量化实验。

### 源码中 t、x、y、z 的直接证据

在 CompuCell3D 这段 solver 中，PDE 子步的时间维度不是显式命名为 `t`，而是由 `extraMCS` 循环体现。也就是说，TemporalStencil 实验中的 `t = 0..T-1` 可以对应到这里连续执行的 diffusion 子步。

文件：`CompuCell3D/CompuCell3D/core/CompuCell3D/steppables/PDESolvers/DiffusionSolverFE_CPU.cpp`，行 `801-807`：

```cpp
for (int extraMCS = 0; extraMCS < scalingExtraMCSVec[i]; extraMCS++) {
    boundaryConditionInit(i);//initializing boundary conditions
    diffuseSingleField(i);
    for (unsigned int j = 0; j < diffSecrFieldTuppleVec[i].secrData.secretionFcnPtrVec.size(); ++j) {
        (this->*diffSecrFieldTuppleVec[i].secrData.secretionFcnPtrVec[j])(i);
    }
}
```

同一个函数的另一个分支也有同样的 diffusion 子步循环。文件：`DiffusionSolverFE_CPU.cpp`，行 `815-818`：

```cpp
for (int extraMCS = 0; extraMCS < scalingExtraMCSVec[i]; extraMCS++) {
    boundaryConditionInit(i);//initializing boundary conditions
    diffuseSingleField(i);
}
```

空间维度 `x, y, z` 直接体现在 `diffuseSingleField` 的三重循环中。文件：`DiffusionSolverFE_CPU.cpp`，行 `944-948`：

```cpp
for (int z = minDim.z; z < maxDim.z; z++)
    for (int y = minDim.y; y < maxDim.y; y++)
        for (int x = minDim.x; x < maxDim.x; x++) {

            currentConcentration = concentrationField.getDirect(x, y, z);
```

邻居访问也是围绕当前 `(x, y, z)` 加 offset 得到 `(offX, offY, offZ)`。文件：`DiffusionSolverFE_CPU.cpp`，行 `966-987`：

```cpp
for (register int i = 0; i <= maxNeighborIndex /*offsetVec.size()*/ ; ++i) {
    const Point3D &offset = offsetVecRef[i];

    int offX = x + offset.x;
    int offY = y + offset.y;
    int offZ = z + offset.z;

    float diffCoef_offset = diffCoef[cellTypeArray.getDirect(offX, offY, offZ)];

    // No diffusion where not diffusive
    if (abs(diffCoef_offset) < FLT_EPSILON) {
        --numNeighbors;
        continue;
    }

    concentrationSum += concentrationField.getDirect(offX, offY, offZ);

}

concentrationSum -= numNeighbors * currentConcentration;

concentrationSum *= currentDiffCoef;
```

最终更新值写入 scratch，再 swap 到下一时间层。文件：`DiffusionSolverFE_CPU.cpp`，关键位置在行 `1133-1134`、`1148`、`1155`：

```cpp
updatedConcentration = (concentrationSum + varDiffSumTerm) +
                       (1 - decayCoef[currentCellType]) * currentConcentration;

concentrationField.setDirectSwap(x, y, z, updatedConcentration);//updating scratch

concentrationField.swapArrays();
```

因此，抽取到 TemporalStencil 后可以建立如下对应关系：

```text
TemporalStencil 的 t 循环       <-> CC3D 的 extraMCS diffusion 子步循环
TemporalStencil 的 x/y/z 循环   <-> CC3D diffuseSingleField 中的 x/y/z 三重循环
TemporalStencil 的邻居访问      <-> CC3D offsetVecRef 生成的 offX/offY/offZ 邻居访问
TemporalStencil 的双缓冲切换    <-> CC3D setDirectSwap(...) + swapArrays()
```

源码中的有限差分注释给出了基本形式：
```text
T_{0, delta tau} = F * sum_i T_i + (1 - N * F) * T_0

```

其中 `N` 是邻居数，`F` 对应 diffusion coefficient。实现中，内部网格点先累加邻居浓度，再减去 `numNeighbors * currentConcentration`，乘以当前 cell type 的 diffusion coefficient，最后加上 decay 后的中心值。

## 2. 最小版本的公式依据

先考虑最简单、最适合与现有 `3d-jacobi-7p` 对接的情况：

* 3D square lattice；
* 只看内部点；
* 一阶最近邻，即 6 个邻居；
* 所有点都可扩散；
* diffusion coefficient 为常数 `D`；
* decay coefficient 为常数 `lambda`；
* 不启用 variable diffusion coefficient；
* 不启用 threshold、secretion、复杂边界条件。

此时 CC3D 的内部点更新可写成：

```text
C^{t+1}_{i,j,k}
  = D * (C^t_{i-1,j,k} + C^t_{i+1,j,k}
       + C^t_{i,j-1,k} + C^t_{i,j+1,k}
       + C^t_{i,j,k-1} + C^t_{i,j,k+1}
       - 6 * C^t_{i,j,k})
    + (1 - lambda) * C^t_{i,j,k}

```

整理后：

```text
C^{t+1}_{i,j,k}
  = D * (C^t_{i-1,j,k} + C^t_{i+1,j,k}
       + C^t_{i,j-1,k} + C^t_{i,j+1,k}
       + C^t_{i,j,k-1} + C^t_{i,j,k+1})
    + (1 - lambda - 6D) * C^t_{i,j,k}

```

因此它就是一个 3D 7 点 Jacobi stencil：

```text
C^{t+1}_{i,j,k} = C0 * C^t_{i,j,k} + C1 * neighbor_sum

```

其中：

```text
C1 = D
C0 = 1 - lambda - 6D

```

这说明最小版本可以直接基于现有 `3d-jacobi-7p` 框架实现。当前 TemporalStencil 里的默认 3D 7 点公式：

```text
C^{t+1} = 0.4 * center + 0.1 * neighbor_sum

```

可以解释为 CC3D diffusion kernel 的一个特例：`D = 0.1`，`lambda = 0`，`C0 = 1 - 6 * 0.1 = 0.4`。

## 3. 为什么先做最小版本

完整 CompuCell3D solver 包含很多与 stencil 核心计算正交的应用逻辑：

* boundary condition refresh；
* cell type dependent diffusion/decay；
* avoid type，也就是某些 cell type 不参与 diffusion；
* variable diffusion coefficient；
* secretion 和 uptake；
* concentration threshold；
* `Array3DContiguous` 的 current/scratch 交错布局；
* OpenMP 分区与 BoxWatcher。

这些特征有真实应用意义，但会掩盖时间向量化本身的收益。最小版本先固定为常系数 3D 7 点 diffusion，可以作为从源码公式到 TemporalStencil kernel 的可验证起点。后续再逐步加入 CC3D 特征，能清楚说明每个特征对性能和向量化难度的影响。

## 4. 建议实验分层

### V0: CC3D-basic diffusion

这是第一版补充实验，目标是建立源码到公式的最短闭环。

特征：

* 3D 7 点 Jacobi；
* 常系数 `D` 和 `lambda`；
* `float` 优先，因为 CC3D 的 concentration field 是 `float`；
* 双缓冲；
* 不包含 cell type、secretion、threshold、复杂边界。

公式：

```text
C^{t+1}_{i,j,k}
  = D * neighbor_sum + (1 - lambda - 6D) * C^t_{i,j,k}

```

论文中可描述为：

```text
We first extract the internal-region, constant-coefficient case of the CompuCell3D forward-Euler diffusion solver. This reduces to a 3D 7-point Jacobi stencil with decay.

```

### V1: decay-aware diffusion

如果 V0 默认 `lambda = 0`，V1 可以显式扫不同 decay coefficient。它仍然是常系数 3D 7 点 stencil，但更贴近 CC3D 的 `DiffusionData::decayCoef` 语义。

公式：

```text
C0 = 1 - lambda - 6D
C1 = D

```

价值：

* 证明不是任意选的 7 点 stencil，而是保留了 CC3D diffusion + decay 的参数结构；
* 可以测试不同中心系数对性能几乎无影响，但对数值语义有应用背景。

### V2: cell-type dependent coefficients

这一版开始加入 CC3D 特征。CC3D 中每个 lattice site 有 cell type，diffusion 和 decay 通过 type 查表：

```text
D_p = diffCoef[type(p)]
lambda_p = decayCoef[type(p)]

```

内部点公式变成：

```text
C^{t+1}_p
  = D_p * (sum_{q in N(p)} C^t_q - N_p * C^t_p)
    + (1 - lambda_p) * C^t_p

```

如果所有邻居都可扩散，`N_p = 6`。整理后：

```text
C^{t+1}_p
  = D_p * sum_{q in N(p)} C^t_q
    + (1 - lambda_p - 6D_p) * C^t_p

```

价值：

* 这是区别于普通 `3d-jacobi-7p` 的主要应用特征；
* 多一个 `cellType` 数组和 coefficient lookup，能反映真实 CC3D kernel 的访存与分支开销；
* 可以讨论时间向量化在非完全常系数 stencil 上的适用性。

### V3: avoid type / non-diffusive sites

CC3D 中如果邻居 cell type 的 diffusion coefficient 接近 0，会跳过该邻居，并减少 `numNeighbors`。抽象公式为：

```text
N_p = count(q in N(p) where D_q != 0)
C^{t+1}_p
  = D_p * (sum_{q in active N(p)} C^t_q - N_p * C^t_p)
    + (1 - lambda_p) * C^t_p

```

价值：

* 这是 CC3D 的真实生物建模特征之一：某些 cell type 可以阻止 diffusion；
* 但它会引入分支或 mask，可能降低 SIMD 效率；
* 适合作为扩展实验，不适合作为第一版。

### V4: boundary condition

CC3D 每个 diffusion 子步前会执行 boundary refresh。最小实验可以先只统计内部区域，或者在 timing 外刷新 halo。后续可加入：

* periodic boundary；
* no-flux boundary；
* constant value boundary；
* constant derivative boundary。

建议论文中把边界作为应用约束说明，不要让边界逻辑成为第一版性能结论的核心。

### V5: CC3D memory layout

CC3D 的 `Array3DContiguous` 把 current 和 scratch 交错存在同一块数组中，通过 `shiftArray` 和 `shiftSwap` 切换当前层。等价地说，它不是普通的 `A[2][x][y][z]`，而是把两个时间层按偏移交错到一个 1D buffer 中。

价值：

* 更接近真实 CC3D；
* 可作为 layout sensitivity experiment；
* 实现复杂度高于普通双数组，因此建议放在 cell-type 版本之后。

## 5. 不建议第一阶段加入的部分

以下功能真实存在，但不适合作为第一阶段补充实验：

* secretion / uptake：它们在 diffusion 子步之间修改 field，会破坏多个纯 diffusion 时间步连续合并的假设；
* variable diffusion coefficient：需要额外计算 `varDiffSumTerm`，公式复杂，适合作为后续扩展；
* threshold：本质是 min/max clamp，会引入非线性和分支；
* BoxWatcher：改变计算区域，不改变 stencil 公式本身；
* hexagonal lattice：邻居结构与 3D square 7 点不同，会分散实验重点。

## 6. 论文表述边界

建议使用如下表述：

```text
We extract an application-derived diffusion kernel from the CPU PDE solver of CompuCell3D. The first version focuses on the internal-region, constant-coefficient forward-Euler update, which reduces to a 3D 7-point stencil with decay. We then incrementally add CompuCell3D-specific features such as cell-type dependent coefficients and non-diffusive sites.

```

不建议使用如下表述：

```text
We accelerate the full CompuCell3D PDE solver.

```

原因是完整 solver 还包括 boundary refresh、secretion、uptake、threshold、variable coefficient、cell field preparation 等逻辑。当前计划更准确地说是：

```text
application-derived kernel extraction and evaluation

```

而不是：

```text
end-to-end CompuCell3D acceleration

```

## 7. 实验价值总结

这个补充实验的价值在于：

1. 它不是凭空构造的 stencil，而是来自 CompuCell3D 的真实 PDE solver。
2. 最小版本能和已有 `3d-jacobi-7p` 直接对齐，便于验证正确性。
3. 后续加入 cell type、avoid type、CC3D layout 后，可以逐步展示真实应用特征对时间向量化的影响。
4. 论文叙述可以从简单公式推导到应用特征扩展，逻辑链条清楚。

建议实施顺序：

```text
V0 constant-coefficient diffusion
-> V1 explicit decay parameter
-> V2 cell-type dependent coefficients
-> V3 avoid type / masked neighbor
-> V4 boundary handling
-> V5 CC3D interleaved layout

```

```

```
