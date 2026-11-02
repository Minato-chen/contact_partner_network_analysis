# 基础计算模块 (`02_basic_calculate`)

本模块计算与具体接触图谱无关的**基础热力学、空间形态与序列物理化学参数**，作为整个分析流水线的公共上游支撑。每个基础指标独立维护在一个 Python 脚本中，支持单指标独立运行或全自动批量执行。

---

## 一、 整体结构与流水线定位

```text
02_basic_calculate/
├── run.py                 # 基础计算批量入口，优先执行相界，自动调度下游指标
├── phase_boundary.py      # 时间平均密度剖面、对称中心拟合与浓相核心半宽
├── csat.py                # 稀相临界饱和浓度与时间块误差分析
├── rg.py                  # 全体系全链质量加权回转半径
├── kappa.py               # 基于 localCIDER 计算电荷不对称分布参数 κ
├── scd.py                 # 序列电荷退化度参数 (Sequence Charge Decoration)
├── _sequence.py           # 居中拓扑单链氨基酸序列读取与校验工具
├── residue_masses.csv     # CALVADOS2 力场标准氨基酸质量表
└── BASIC_CALCULATE_README.md
```

---

## 二、 严谨指标命名与物理定义

| 指标命名 | 英文标量 Key | 单位 | 物理定义与科学含义 |
| :--- | :--- | :--- | :--- |
| **浓相核心半宽** | `phase_boundary.dense_half_width_nm` | `nm` | 时间平均密度剖面中，凝聚相主体半高宽的 $0.6$ 倍（$0.6 w$）。作为下游划分浓相单链系综（`dense_chains`）的客观几何边界。 |
| **临界饱和浓度** | `csat.csat_mM` | `mM` (链等效浓度) | 选取远离中心凝聚体（$|z| > 2w$）的两侧稀相储库区域，统计游离单链珠数并换算为摩尔浓度。直接量化相分离驱动力（$C_{\mathrm{sat}}$ 越小，相分离能力越强）。 |
| **全系综回转半径** | `rg.rg_nm` | `nm` | 对模拟盒子中全部链（包含浓相与稀相）计算质量加权均方根回转半径 $R_g$ 并取时间系综平均。表征体系宏观平均构象尺寸。 |
| **电荷图案参数 $\kappa$** | `kappa` | 无量纲 (`dimensionless`) | 基于 localCIDER 严格计算的正负电荷交替/聚集度参数。$\kappa \to 0$ 表示正负交替最均匀；$\kappa \to 1$ 表示强电荷区块化分离。 |
| **序列电荷退化度 (SCD)** | `scd` | 无量纲 (`dimensionless`) | $\text{SCD} = \frac{1}{N}\sum_{i < j} q_i q_j \sqrt{j - i}$。表征沿序列一维空间电荷排布对多肽链长程静电吸引/排斥贡献的理论标度。 |

---

## 三、 核心算法与计算原理

### 1. 密度分布与自适应相界探测 (`phase_boundary.py`)
- **空间离散化**：沿垂直于 slab 界面的 $z$ 轴以 $\Delta z = 0.1\text{ nm}$ 进行精细分箱，统计残基珠数密度并归一化为链等效浓度（$\text{mM}$）。
- **中心与峰高**：以 slab 几何中心为原点，取中心及左右各 5 箱（共 11 箱，跨度 $1.1\text{ nm}$）的时间平均浓度均值作为凝聚相峰值浓度 $C_{\mathrm{peak}}$。
- **自适应半高探测机制**：
  - **常规体系**：向右扫描，首个浓度低于峰高一半（$C \le C_{\mathrm{peak}} / 2$）的箱中心位置记为半高截断距离 $w$。
  - **弱相分离体系自适应跨越**：对于稀相基线浓度较高的高溶解弱体系（如 C22，稀相基准 $\sim 4.6\text{ mM}$，导致峰半高低于基线），自动采用峰顶与稀相基线的中点作为过零截断，保证对任意极端弱相分离体系物理稳健不回退。
- **浓相核心划分**：定义浓相核心边界为 $|z| < 0.6 w$。

### 2. 临界饱和浓度分析 (`csat.py`)
- **储库取样区间**：严格选取 $|z| > 2w$ 的两端平坦区域（而非 $2 \times 0.6w$），完全避开相界面过度毛刺区。
- **逐帧等权统计**：每帧累加储库内残基珠数，除以链长（129残基）与储库体积（$L_x L_y (L_z - 4w)$），换算为 $\text{mM}$。
- **误差分析（Block SEM）**：以 50 帧为一个独立时间块，计算块间标准差并除以 $\sqrt{N_{\mathrm{blocks}}}$ 得到时间自相关修正后的统计标准误。

### 3. 全链质量加权回转半径 (`rg.py`)
- **定义**：
  $$R_g = \sqrt{\frac{\sum_i m_i \|\mathbf{r}_i - \mathbf{r}_{\mathrm{com}}\|^2}{\sum_i m_i}}$$
  其中 $m_i$ 来自 CALVADOS2 力场残基质量表，N端与C端分别补偿末端原子质量（$+2\text{ Da}$ / $+16\text{ Da}$）。整链坐标基于周期性最小镜像整链重组，消除跨盒断链假象。

### 4. 纯序列指标 (`kappa.py` & `scd.py`)
- 仅依赖多肽一级序列，从居中拓扑中校验所有链的序列一致性后计算一次。
- 属于严格静态序列固有属性，无时间误差 SEM。

---

## 四、 命令行执行指南

在工程根目录下执行：

```bash
# 一键运行全部基础指标
python 02_basic_calculate/run.py

# 指定单个模拟执行相界与 csat
python 02_basic_calculate/phase_boundary.py --simulation B01_T260
python 02_basic_calculate/csat.py --simulation B01_T260

# 运行全系综回转半径
python 02_basic_calculate/rg.py
```

---

## 五、 输出规范与数据流转

- **明细输出路径**：`03_results/<simulation_id>/01_basic/`
  - `phase_boundary/`：包含 `density_profile.csv`（全 $z$ 轴剖面）、`boundary.json`（核心半宽与采样参数）。
  - `csat/`：包含 `concentration_frames.csv`（逐帧稀相浓度）、`summary.csv`（系综均值与块SEM）。
  - `rg/`：包含全链各帧 $R_g$ 记录与系综汇总。
  - `kappa.scalars.json`、`scd.scalars.json`：单值标量标准记录。
- **标量自动汇总**：运行 `python 04_summary/run.py` 即可自动提取生成 `04_summary/output/basic.csv`，包含全模拟的对比基础参数。
