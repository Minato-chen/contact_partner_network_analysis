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

脚本按 A18、B01、B10 的顺序处理。它从原轨迹的**零起始索引 750** 开始，每隔 3 帧取一帧：750、753、…、2997；预期每组输出 750 帧，保存为 `centered/{变体}_T260_centered.dcd`。具体居中由 `dependencies/center_local.py` 执行。程序检查取帧数是否为750，并在写出后重新读取全部帧，检查盒子、有限坐标、相邻骨架距离和链中心位置。

- `[ok] A18: ...`：A18 的居中 DCD 已写入。
- `[skip] verified trajectory and settings: ...`：输入记录、参数、算法源码指纹、输出DCD及诊断文件的校验均匹配，才复用。旧版没有记录的DCD会报错，需选择新输出目录或显式使用`--force`。
- `Step 1 complete`：三组都已处理或跳过，可以进入第二步。

新版本建议先用`python 01_center.py --output-dir centered_verified`保存到独立目录，保留旧`centered/`。后续使用`python 02_compute.py --centered-dir centered_verified --results-dir results_verified`读取新轨迹并单独输出结果。`notebook/`中若要画新版结果，需将`RESULTS`指向`../results_verified`，避免仍读取旧图数据。

只处理 B01：`python 01_center.py --variant B01`。重新生成已有文件：`python 01_center.py --force`。两者也可组合，例如 `python 01_center.py --variant B01 --force`。

若重新生成任何居中 DCD，随后应运行 `python 02_compute.py --force`，避免新轨迹配旧计算结果。

MDAnalysis 可能打印有关 `xdrlib` 或 `DCDReader` 的 `DeprecationWarning`。它们是依赖库的未来兼容性提示；请以各变体的 `[ok]`／`[skip]`、最后的 `Step 1 complete` 和进程退出状态判断本次是否完成。警告出现后脚本可能仍在处理下一组轨迹。

### 居中到底做了什么（易懂说明）

可以把模拟盒子想象成一根长管，蛋白质较密集的一团会沿长方向移动。居中的目的，是把每帧中这团蛋白质放回盒子中部，方便后续在相同位置统计；不是把某一条中心链固定住。

1. **选画面**：对3000帧原轨迹，跳过零索引0–749，之后每3帧取一次，共750帧。逐帧读取原轨迹，写出时只保留这些帧。第一步会核对实际取帧数与`EXPECTED_FRAMES=750`一致。
2. **粗对齐**：沿z方向把盒子切成厚1 Å（0.1 nm）的薄片，数每片里的残基珠。把珠最多的一片移到盒子中部；所有珠一起平移。
3. **精对齐**：重新数薄片中的珠，把连续非空的薄片看作一个区域。选择珠数总和最多的区域，计算其中珠的平均z位置，再把这个位置移到盒子中部。这里是按珠数加权，不是按质量加权，也不是根据两侧界面求中点。
4. **处理周期边界**：定位密度前以及每次平移后，将珠的周期位置放回盒内。最后沿每条链的残基顺序，将相邻骨架位移换成周期最短位移，逐珠恢复连续链；再整链平移，使几何中心位于主盒内。完整链的个别珠允许落在盒外，不能再逐珠wrap。
5. **保存**：先写临时DCD，完成后再替换到`centered/`。原始轨迹不改动。

当前盒长150 nm，所以居中目标是z=75 nm，不是z=0；后续选链时才减去半盒长。整个过程只沿z平移，没有旋转体系，也没有在此阶段决定哪些链属于浓相。

**适用范围**：固定正交盒、每个segment对应一条按残基顺序排列的线性链、每残基一珠。盒长变化、非正交盒、不合理骨架距离或帧数不符时拒绝输出。居中不证明相分离，也不保证多个片层时一直选择同一片层；不应用居中轨迹直接计算绝对扩散。

**2026-10-22优化**：保留原有密度峰及连续区域的两次对齐，改用沿骨架的周期重建绕过本地unwrap提前返回问题。原有文件不自动覆盖。`max_bond_nm`默认1.0 nm仅作宽松的骨架/拓扑检查，和残基接触阈值无关；不应为掩盖异常任意放宽。

每个DCD另有两份文件：
- `*.diagnostics.csv`：原始帧号、两次平移量、区域范围/宽度、珠数占比、第二大区域占比、是否触及盒边、写出前重新分箱后识别的主要区域中心偏移、最大相邻珠距离。逐帧查看这些数据可检查区域切换或稀区牵引；不把一次对齐成功当作相分离证明。`final_patch_center_offset_nm`是重新分箱后再选主要非空区域的偏移，用于检查区域识别是否稳定；不是已独立确定的物理浓相中心误差。A18/B10部分帧此项偏大，需结合密度剖面判断，不能据此自动删帧或重定选区。
- `*.metadata.json`：输入绝对路径、大小、修改时间、拓扑SHA256、取帧参数、算法源码SHA256、输出与诊断文件SHA256。原始DCD使用文件大小/修改时间识别，不是完整内容哈希；如果人为保留大小和时间修改原DCD，此记录不能识别。metadata最后发布，作为完成记录。

程序逐帧读取而非将整条DCD载入内存，在临时目录写出并复核，通过后才发布文件。三份文件不是一个原子事务，但中断或错配时下一次复用校验会拒绝。不要同时对同一输出运行多个写入进程。

### 与ana_v4及本地CALVADOS代码的区别（仅比较居中）

本段依据本地`/Users/minato/Desktop/slab_idr/ana_v4/unified_analysis/geometry.py`及`/Users/minato/Desktop/slab_idr/calvados/analysis.py`。本地CALVADOS文件的版本尚未与服务器或上游发行版核对，因此下述差异只对应这些具体文件。

| 项目 | 本项目contact_partner | ana_v4 | 本地CALVADOS `SlabAnalysis.center` |
|---|---|---|---|
| 主要步骤 | 密度峰粗对齐→主要连续区域精对齐→沿骨架恢复完整链 | 两次密度对齐相同，末尾调用unwrap | 两次密度对齐相同，末尾调用unwrap |
| 用谁确定密度中心 | 全部珠 | 全部珠 | 可选参考组分`ref`或全部珠`all`；参考组分不一定是一条链 |
| 帧设置 | 固定start=750、step=3 | 默认start=总帧数的1/4、step=3，可指定stop | 接口默认start/end=None、step=1，可由调用者指定 |
| 密度薄片 | 1 Å | 1 Å | 本地`calc_z_Angstr`为1 Å |
| 连续区域选择 | 按非零箱的珠数总和选择；循环结束后再检查最后一段 | 同样补查最后一段 | 本地`calc_zpatch`未补查结束时仍连续非零的最后一段 |
| 没有可用区域时 | 报错，拒绝输出 | 不执行第二次平移 | 直接对区域求加权平均，无同样保护 |
| 骨架连接 | 按每segment相邻索引配对 | 按索引i连接i+1 | 按索引i连接i+1；当前连续索引拓扑下相同 |
| 输出与复用 | 独立目录，参数/来源/文件校验匹配才复用，写后逐帧复核 | 主流程调用时写出到配置的output目录 | 写到该SlabAnalysis设置的input_path/centered_dcd |

因此，三者采用的是同一类slab居中思路，不是三套独立的物理方法。当前单组分体系若CALVADOS参考组包含全部链，`ref`与`all`可使用相同的珠；多组分时则可能产生不同中心。contact_partner与ana_v4对末段的补查是相对于本地CALVADOS文件的实现差异，不代表已验证所有帧输出都相同或不同。

CALVADOS文件还有一个名字相近的`center_traj`函数：它只按全部珠的几何中心整体移动到盒中心，不进行上述密度峰和连续区域的两次对齐。这里对应比较的是`SlabAnalysis.center`，不要将两者混淆。

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
| `chain_partner_frame_means.csv` | 每帧稠相链的平均指标、稠相链数及有效比值链数 | 图a核对、逐帧波动检查 |
| `chain_partner_summary.csv` | 每帧等权均值、逐帧SD、时间块SD与SEM | 数值核对；notebook 01按同一函数重新汇总逐链CSV |
| `partner_count_distribution.csv` | 逐帧归一化后平均的伙伴数分布 | notebook 02 |
| `partner_count_distribution_by_block.csv` | 各时间块的伙伴数分布 | notebook 02 的波动带 |
| `partner_count_time_blocks.csv` | 稠相链数与平均伙伴数的时间块诊断数据 | 辅助检查；当前 PPT 图不读取 |
| `{变体}_edge_arrays.npz` | 每帧每对链的接触数，以及双方是否均在稠相的标记 | notebook 03 |

CSV 文件名的 `{变体}` 后缀为 `_T260`；NPZ 文件名没有该后缀。CSV 和 NPZ 是同级计算产物；不能只用逐链 CSV 重建 notebook 03 所需的链对数组。

`[compute]` 表示正在计算，`[reuse]` 表示对应文件已存在并被复用；最后的 `[ok] 750 frames per variant` 表示第二步完成。复用逐链 CSV 时，脚本仍会更新汇总与分布 CSV。`[reuse]` 也只根据文件存在判断。完整重算三组结果：

```bash
python 02_compute.py --force
```

第二步目前没有只重算单个变体的选项。可通过`--centered-dir`和`--results-dir`指定独立输入/输出；使用非默认轨迹目录时必须同时指定非默认结果目录。第二步旧结果的复用仍以存在为依据，改变轨迹/距离/选区后应使用新结果目录或`--force`，不能仅依赖第一步的校验。

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

## 结果1的平均与误差定义

每个变体、伙伴阈值分别汇总，顺序为：**每帧内对稠相链求平均 → 对750个帧平均值等权平均**。稠相链较多的帧不会获得更大权重。notebook 01调用 `dependencies/partner_analysis.py` 的同一绘图与汇总函数：图a每帧一个点；图b中心值为750帧平均，误差线为时间块SEM。

- `*_mean`：帧平均值的直接平均。
- `*_std_frames`、`*_count_frames`：帧平均值的样本标准差（ddof=1）和有效帧数。
- `*_std_blocks`：连续50个分析帧的组平均值之间的样本标准差。
- `*_sem_blocks`：组平均值的SD除以有效组数的平方根；通常15组。
- 原字段 `*_std`、`*_count` 保留，仍分别表示组平均值的SD、有效组数。
- `n_frames`、`n_blocks`、`block_frames_min/max`：实际帧数、组数、最短/最长组的帧数。

均值直接从帧级结果计算，不依赖分组长度。50帧是当前分组设置，尚不能据此认定各组独立；分组SEM是同一轨迹内的误差估计，使用于正式推断前需检查相关时间或分组长度稳定性。逐帧SD描述时间波动，不能与SEM互换。notebook 02、03原有的时间块SD波动带继续保留。

`contacts_per_partner`仍是各链接触数/伙伴数比值的平均，再进行帧平均；不是两个总体均值之比。零伙伴链的比值为NaN，仅此指标排除这些链。帧级表中 `n_valid_contacts_per_partner` 记录有效链数，汇总表的有效帧/组数可能少于750/15；完全没有有效值时，该指标及误差保留NaN。

更新代码后，运行 `python 02_compute.py`（不加 `--force`）即可复用现有逐链CSV与NPZ，更新帧级表、汇总表和既有分布表；然后从notebook目录重启内核，重新运行notebook 01。无需重新居中或计算接触。本次代码修改本身不会更新已有结果CSV或图片，重跑前它们仍是旧版产物。

## 文件夹位置与依赖文件

分析文件夹可以整体移动，只需保持 `centered/`、`results/`、`figures/` 与 `notebook/` 的相对关系。原始数据路径不会自动改变，因为 `DATA_ROOT` 是绝对路径。

`dependencies/center_local.py` 由第一步调用；`dependencies/partner_analysis.py` 和 `dependencies/contact_persistence_test.py` 由第二步调用。正常流程只运行顶层的两个脚本和所需 notebook，无需单独运行依赖文件。脚本使用临时英文路径处理 DCD，以兼容含日文字符的目录名。

### 接触距离的边界约定

残基接触采用严格 `d < 1.0 nm`。邻域搜索先取得候选点对，再用周期边界下的最小镜像距离排除恰好等于阈值的点对；伙伴链与网络边仍采用接触组数 `>= m`。本次对现有三组各750帧逐一核验，严格边界修正未改变任何链对接触计数，因此现有CSV、NPZ和结果数值可继续使用。以后若修改距离、稠相边界或轨迹，应使用 `--force` 重新计算。
