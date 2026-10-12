# Tau K18 T260：接触数、伙伴链数与稠相网络分析

本文件夹用于从 A18、B01、B10 三组原始轨迹复现 PPT 对应的三组图。流程为：

```text
原始 DCD + top.pdb
  → 01_center.py → centered/*_centered.dcd
  → 02_compute.py → results/*.csv.gz、*.csv、*.npz
  → notebook/ 中的 01、02、03 → figures/ 中的 PPT 用图
```

原始 DCD 和 `top.pdb` 不复制到本文件夹。两个脚本中的 `DATA_ROOT` 固定为 `/Users/minato/Desktop/slab_idr`。若原始数据搬到别处，需要同时修改 `01_center.py` 和 `02_compute.py` 的 `DATA_ROOT`。

## 1. 环境与输入文件

需要 Python 3.11，以及 NumPy、pandas、SciPy、MDTraj、MDAnalysis、Matplotlib、Jupyter/IPython。先进入分析文件夹：

```bash
cd /Users/minato/Desktop/slab_idr/contact_partner_network_analysis
```

每个变体的原始文件应位于以下位置（将 `A18` 换为 `B01`、`B10` 即是另两组）：

```text
/Users/minato/Desktop/slab_idr/tau_k18_129aa_A18_T260_5/top.pdb
/Users/minato/Desktop/slab_idr/tau_k18_129aa_A18_T260_5/tau_k18_129aa_A18_T260_5.dcd
```

缺少输入时，脚本会报 `FileNotFoundError` 并指出缺少的路径。

## 2. 第一步：生成居中轨迹

```bash
python 01_center.py
```

脚本按 A18、B01、B10 的顺序处理。它从原轨迹的**零起始索引 750** 开始，每隔 3 帧取一帧：750、753、…、2997；预期每组输出 750 帧，保存为 `centered/{变体}_T260_centered.dcd`。具体居中由 `dependencies/center_local.py` 执行。

- `[ok] A18: ...`：A18 的居中 DCD 已写入。
- `[skip] A18: ...`：目标文件已存在，因此跳过；这里只检查是否存在，不重新验证文件内容。
- `Step 1 complete`：三组都已处理或跳过，可以进入第二步。

只处理 B01：`python 01_center.py --variant B01`。重新生成已有文件：`python 01_center.py --force`。两者也可组合，例如 `python 01_center.py --variant B01 --force`。

若重新生成任何居中 DCD，随后应运行 `python 02_compute.py --force`，避免新轨迹配旧计算结果。

MDAnalysis 可能打印有关 `xdrlib` 或 `DCDReader` 的 `DeprecationWarning`。它们是依赖库的未来兼容性提示；请以各变体的 `[ok]`／`[skip]`、最后的 `Step 1 complete` 和进程退出状态判断本次是否完成。警告出现后脚本可能仍在处理下一组轨迹。

## 3. 第二步：计算 CSV 与 NPZ

确认三组居中 DCD 已存在，然后运行：

```bash
python 02_compute.py
```

脚本读取 `centered/` 的 DCD 和原始目录中的 `top.pdb`。链间接触距离阈值为 **1.0 nm**，只将稠相内部链之间的接触计入分析。A18、B01、B10 的稠相中心半宽分别为 **11.25、8.91、9.81 nm**。逐链伙伴数分别按链对至少有 **1、5、10** 个接触计算；汇总使用连续 **50 帧**为一个时间块。

`results/` 中的产物：

| 文件 | 内容 | 用途 |
|---|---|---|
| `{变体}_T260_chain_partner_rows.csv.gz` | 逐帧、逐稠相链、逐伙伴阈值的接触数与伙伴数 | notebook 01；生成分布 CSV |
| `chain_partner_summary.csv` | 时间块汇总表 | 数值核对；当前三张 PPT 图不直接读取 |
| `partner_count_distribution.csv` | 逐帧归一化后平均的伙伴数分布 | notebook 02 |
| `partner_count_distribution_by_block.csv` | 各时间块的伙伴数分布 | notebook 02 的波动带 |
| `partner_count_time_blocks.csv` | 稠相链数与平均伙伴数的时间块诊断数据 | 辅助检查；当前 PPT 图不读取 |
| `{变体}_edge_arrays.npz` | 每帧每对链的接触数，以及双方是否均在稠相的标记 | notebook 03 |

CSV 文件名的 `{变体}` 后缀为 `_T260`；NPZ 文件名没有该后缀。CSV 和 NPZ 是同级计算产物；不能只用逐链 CSV 重建 notebook 03 所需的链对数组。

`[compute]` 表示正在计算，`[reuse]` 表示对应文件已存在并被复用；最后的 `[ok] 750 frames per variant` 表示第二步完成。复用逐链 CSV 时，脚本仍会更新汇总与分布 CSV。`[reuse]` 也只根据文件存在判断。完整重算三组结果：

```bash
python 02_compute.py --force
```

第二步目前没有只重算单个变体的选项。

## 4. 从 notebook 生成 PPT 图

**必须从 `notebook/` 目录启动 Jupyter**。三个 notebook 通过 `../results/` 读入数据，并把图写到 `../figures/`：

```bash
cd /Users/minato/Desktop/slab_idr/contact_partner_network_analysis/notebook
jupyter lab
```

| Notebook | 输入 | 最终 PNG |
|---|---|---|
| `01_接触数とパートナー鎖数_V26.ipynb` | 三组逐链 CSV | `figures/01_partner/24_chain_partner_test.png` |
| `02_パートナー鎖数の分布_V26.ipynb` | 分布及时间块分布 CSV | `figures/02_distribution/partner_count_distribution_threshold1_5.png` |
| `03_稠相ネットワーク閾値解析_V26.ipynb` | 三组链对 NPZ | `figures/03_network/graph_threshold_for_slide.png` |

三份 notebook 相互独立，无需依次运行；每张图还会导出 SVG、PDF。`figures/` 是可重新生成的绘图输出，不是第二步的输入。表中的 `figures/` 路径以分析文件夹为起点。

## 与 PPT 的数值定义

- notebook 01：a 图的伙伴阈值是 **1**；b 图比较阈值 **1、5、10**。现有 V35 PPT 的旧表格及旧图可能采用阈值 5；替换图片时应同步核对标题、图注和讲稿。
- notebook 02：展示阈值 **1 和 5** 的伙伴数分布；波动带来自时间块之间的标准差。
- notebook 03：从 NPZ 计算网络全局效率和最大连通分量比例，绘制 PPT 用双面板图。波动带描述时间块间差异，不是独立模拟重复的误差。

## 文件夹位置与依赖文件

分析文件夹可以整体移动，只需保持 `centered/`、`results/`、`figures/` 与 `notebook/` 的相对关系。原始数据路径不会自动改变，因为 `DATA_ROOT` 是绝对路径。

`dependencies/center_local.py` 由第一步调用；`dependencies/partner_analysis.py` 和 `dependencies/contact_persistence_test.py` 由第二步调用。正常流程只运行顶层的两个脚本和所需 notebook，无需单独运行依赖文件。脚本使用临时英文路径处理 DCD，以兼容含日文字符的目录名。
