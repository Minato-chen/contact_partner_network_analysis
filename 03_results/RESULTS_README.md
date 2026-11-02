# 数值结果目录 (`03_results`)

本目录保存各个独立模拟体系（7 个变体）经过上游计算流水线生成的全量明细数据、构象时序、关联图谱与轻量标量指标。

---

## 一、 体系与目录层级架构

第一层为各模拟系统的唯一标识目录，第二层严格遵循两阶段标准计算结构：

```text
03_results/
├── <simulation_id>/               # 7 个系统：A18_T260, B01_T260, B10_T260, A20_T260, C04_T260, B09_T260, C22_T260
│   ├── 01_basic/                  # 基础构象与宏观相行为结果
│   │   ├── phase_boundary/        # 平板对称拟合、界面宽度、密度剖面及 Csat (含 COMPLETE)
│   │   ├── rg/                    # 全体系与浓相链回转半径时间序列 (含 COMPLETE)
│   │   ├── charge_hydropathy/     # 序列电荷与疏水性参数 (含 COMPLETE)
│   │   ├── sequence_properties/   # 严格生化参数 kappa, SCD 等 (含 COMPLETE)
│   │   └── basic.scalars.json     # 汇总结算的轻量标量键值字典
│   └── 02_contacts/               # 分子内外接触、多价伙伴与拓扑网络结果 (0.8 nm 严格定义)
│       ├── 01_intra/              # 分子内构象紧实度与残基接触谱 (含 dense_chains 与 central_chain)
│       ├── 02_inter/              # 分子间电荷拆解、多价伙伴与网络逾渗
│       ├── 03_joint/              # 分子内外接触平衡与联合多价景观
│       ├── metadata.json          # 接触计算元数据
│       ├── COMPLETE               # 状态标记锁
│       └── contacts.scalars.json  # 汇总结算的轻量标量键值字典 (含各对象 77 个严谨指标)
├── verify_results.py              # 自动化结果完整性检验工具
└── RESULTS_README.md              # 本技术文档
```

---

## 二、 结果验证工具 (`verify_results.py`)

为确保下游 `04_summary` 收集与 `05_plot` 绘图万无一失，提供了一键自动化校验脚本：

```bash
# 验证全部 7 个变体的 basic 与 contacts 输出完整性
python 03_results/verify_results.py
```

### 校验项：
1. **01_basic 阶段**：检验 `phase_boundary/COMPLETE`、`rg/COMPLETE`、`charge_hydropathy/COMPLETE`、`sequence_properties/COMPLETE` 与 `basic.scalars.json`。
2. **02_contacts 阶段**：检验 `COMPLETE` 锁、`metadata.json`、`contacts.scalars.json`，以及 `01_intra/`、`02_inter/`、`03_joint/` 下所有明细时序与矩阵。

---

## 三、 使用与维护规范

1. **上游写保护**：明细结果均由 `02_basic_calculate` 与 `02_calculate/01_contacts` 自动产出，严禁手工修改任何 `.csv`、`.json` 或 `.npy` 文件。
2. **状态锁守则**：只有在对应任务完全成功后系统才会写入 `COMPLETE` 标记。下游模块如果发现缺少 `COMPLETE`，会自动发出警报并阻止脏数据污染汇总表。
3. **成果交付分流**：详细明细数据保存在本目录；全局横向横切汇总表在 `04_summary/output/`；科学论文级别可视化成果在 `05_plot/`。
