# 结果汇总与对比表聚合模块 (`04_summary`)

本模块作为整个分析流水线的**第四阶段（Stage-4）聚合引擎**，专门用于从 `03_results/` 中收集已完成的标量记录，自动化生成跨体系横向对比的宽表（Wide Tables）。本模块不扫描轨迹，不执行任何数值拟合，纯粹作为严格的数据汇聚与 Schema 校验层。

---

## 一、 整体结构与输出产物

```text
04_summary/
├── run.py                 # 汇总收集主脚本，自动遍历并严格对齐标量
├── output/                # 最终汇总交付表
│   ├── basic.csv                         # 基础热力学与序列指标横向对比表 (1行/模拟)
│   ├── contacts_dense_chains.csv         # 浓相全系综接触指标对比表 (含第一/二/三层指标及 SEM)
│   ├── contacts_central_chain.csv        # 中心单链探针接触指标对比表
│   ├── status.csv                        # 77 项标量收集状态诊断报告 (SUCCESS/MISSING/STALE)
│   └── sources.json                      # 各指标溯源记录、计算参数与元数据快照
└── SUMMARY_README.md
```

---

## 二、 核心聚合机制与数据流转

### 1. 分类映射与表名派生
- 遍历 `03_results/<simulation_id>/` 下所有 `scalars.json` 与 `*.scalars.json`。
- **表名自动生成规则**：剥离目录首部的阶段序号（`re.sub(r"^\d+_", "", category)`）：
  - `01_basic` $\longrightarrow$ `basic.csv`（基础指标，无区分对象）。
  - `01_contacts` + `object=dense_chains` $\longrightarrow$ `contacts_dense_chains.csv`。
  - `01_contacts` + `object=central_chain` $\longrightarrow$ `contacts_central_chain.csv`。

### 2. 严格的 Schema 强一致性防御机制
- **单位与定义锁（Definition Lock）**：同一指标在不同模拟体系中，必须拥有完全一致的物理单位（`unit`）与参数定义字典（`definition`，如截断距离 `0.8 nm`、分箱宽度等）。若定义不一致，立即熔断终止，杜绝不可比数据被拼接为同一列。
- **条件参数标签化**：对于多阈值指标（如 `threshold=1, 5, 10`），自动在列名中以中括号展开，如：
  `inter.partner_count [threshold=1] (chains)` 与 `inter.partner_count [threshold=1] (chains) SEM`。
- **缺失与异常处理**：若某模拟由于未完成（缺少 `COMPLETE` 标记）或上游过期（`STALE`），该模拟指标在对比表中以 `NaN` 占位，并在 `status.csv` 中明确标出。

---

## 三、 输出宽表字段速查

### 1. `basic.csv`
- `phase_boundary.dense_half_width_nm`：浓相核心半宽（nm）。
- `csat.csat_mM` 及 SEM：稀相饱和浓度（mM）。
- `rg.rg_nm` 及 SEM：全盒子质量加权回转半径（nm）。
- `kappa`：序列电荷不对称分布参数（无量纲）。
- `scd`：序列电荷退化度（无量纲）。

### 2. `contacts_dense_chains.csv` 与 `contacts_central_chain.csv`
- **分子内维度 (`01_intra`)**：
  - `intra.contact_count`、`intra.rg_nm`、`intra.long_range_fraction`、`intra.global_efficiency`、`intra.largest_component_fraction`。
- **分子间维度 (`02_inter`)**：
  - **接触数**：`inter.total_contacts`、`inter.attract_contacts`、`inter.repulse_contacts`、`inter.neutral_contacts`、`inter.attract_fraction`。
  - **伙伴数**（按 `threshold=1, 5, 10` 展开）：`inter.partner_count`、`inter.contacts_per_partner`、`inter.attract_contacts_per_partner`、`inter.max_attract_contacts_per_partner`。
  - **网络**（按 `threshold=1, 5, 10, 20, 30, 35` 展开）：`network.global_efficiency`、`network.largest_component_fraction`。
- **联合维度 (`03_joint`)**：
  - `joint.intra_contacts_per_residue`、`joint.inter_contacts_per_residue`、`joint.inter_fraction`、`joint.inter_intra_ratio`。

---

## 四、 命令行操作指南

```bash
# 全自动收集并更新 04_summary/output/ 下的所有表格
python 04_summary/run.py

# 快速验证收集状态（要求全部 77/77 个标量呈现 SUCCESS）
python -c '
import pandas as pd
df = pd.read_csv("04_summary/output/status.csv")
print(df["status"].value_counts().to_dict())
'
```
