# 轨迹居中与整链重组模块 (`01_center`)

本模块作为整个物理分析管线的**第一阶段（Stage-1）数据预处理核心**，负责从原始粗粒化分子动力学轨迹中提取稳态平衡帧、将 slab 凝聚相平移至周期性模拟盒子几何中心，并修复高分子链跨越周期性边界（PBC）导致的分子割裂。

---

## 一、 整体结构与关键文件

```text
01_center/
├── run.py                 # 批量执行入口与全局截取参数配置
├── dependencies/
│   └── center_local.py    # 质心平移、最大连续核心检测与 PBC 整链重组算法
├── output/                # 各模拟体系标准化输出目录
│   └── <simulation_id>/
│       ├── centered.dcd               # 居中并恢复完整分子后的二进制轨迹 (单位: Å)
│       ├── topology.pdb              # 标准居中拓扑结构文件
│       ├── chains.json               # 各单链残基原子索引映射表
│       ├── frame_index.csv           # 居中分析帧与原始轨迹帧号及物理时间对照表
│       ├── centered.diagnostics.csv   # 居中漂移量与质心偏移诊断日志
│       ├── metadata.json             # 上游元数据、盒尺寸与版本记录
│       └── COMPLETE                  # 完成校验标志文件
└── CENTER_README.md
```

---

## 二、 核心算法与物理处理机制

### 1. 稳态区间取样策略
- 原始分子动力学模拟由 `simulations.yaml` 统一定义，原始轨迹只读不修改。
- 默认采样参数（在 `run.py` 顶部定义）：
  - `start = 750`（跳过前 750 帧升温与初始相分离不稳定期，确保进入热力学稳态）。
  - `step = 3`（每隔 3 帧取样一次，消除高频时间自相关，兼顾相空间采样密度）。
  - `expected_frames = 750`（标准生产段共提取 750 帧有效构象用于统计系综）。

### 2. 双阶段空间居中定位
- **阶段一（密度粗平移）**：沿 slab  법向（$z$ 轴）对残基密度进行空间分箱，计算密度峰值位置，并将整体坐标沿 $z$ 轴平移，使凝聚相主峰粗略对齐到盒中心 $z = L_z / 2$。
- **阶段二（最大连续占据核心加权微调）**：检测处于浓相内部连续占据区域的质心，进行精细微调，彻底消除 slab 整体在布朗运动中的随机漂移。

### 3. PBC 周期性边界跳跃修复与整链重建
- 在 OpenMM / LAMMPS 等引擎输出中，分子跨越模拟盒子边界时会被空间折叠（Wrap），导致高分子链被割裂在盒子两端。
- 算法以每条链的第一个残基为空间锚点，沿多肽共价主链依次应用最小镜像条件（Minimum Image Condition）：
  $$\Delta \mathbf{r}_{i, i-1} = (\mathbf{r}_i - \mathbf{r}_{i-1}) - \mathbf{L} \cdot \left\lfloor \frac{\mathbf{r}_i - \mathbf{r}_{i-1}}{\mathbf{L}} + 0.5 \right\rfloor$$
  保证相邻残基虚假键长不超过 $1.0\text{ nm}$，完全恢复整条共价链在空间中的真实物理连续构象。

---

## 三、 输出规约与字段对照

| 产物文件名 | 格式 / 单位 | 科学用途与数据定义 |
| :--- | :--- | :--- |
| `centered.dcd` | DCD 二进制 / 坐标 $\text{\AA}$ | 居中后的连续轨迹。下游 MDTraj 读取时自动转换为国际标准单位 $\text{nm}$。 |
| `chains.json` | JSON 数组 | 记录 $M$ 条链在全体系原子列表中的起始与结束索引（如 $[[0, \dots, 128], [129, \dots]]]$）。 |
| `frame_index.csv` | CSV 映射表 | 包含三列：`frame`（居中连续分析帧索引 $0 \sim 749$）、`source_frame`（原始模拟帧号）、`time_ps`（真实物理纳秒时间戳）。 |
| `centered.diagnostics.csv` | CSV 诊断日志 | 记录每帧 $z$ 轴移动距离与残余偏移，验证是否成功居中。 |

---

## 四、 命令行操作

```bash
# 处理 simulations.yaml 中定义的全部模拟体系
python -u 01_center/run.py

# 仅处理单个模拟体系
python -u 01_center/run.py --simulation B01_T260
```
