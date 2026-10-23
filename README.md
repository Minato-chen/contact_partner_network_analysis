# 分子内与分子间接触分析

按阶段处理模拟轨迹：**居中 → 密度与浓相边界 → 分子内/分子间接触 → Notebook 绘图**。各模块独立读写文件，原始数据只读；分析参数放在对应入口，不放进模拟清单。

## 1. 目录

```text
simulations.yaml                    原始模拟清单
01_center/
  run.py                            选帧与居中参数
  dependencies/                     居中实现
  output/<模拟编号>/
02_calculate/
  run.py                            依次运行边界和两类接触
  dependencies/analysis_io.py        读取轨迹、周期距离、统计
  phase_boundary/                    密度与浓相核心边界
  intramolecular/                    分子内接触与残基网络
  intermolecular/                    分子间接触、伙伴与链网络
```

第二步的三个计算目录各有 `run.py`、`dependencies/`、`output/` 和 `notebooks/`。`project.py` 只提供清单读取和入口工具；`result_files.py` 提供简单文件操作，没有统一 workflow。

## 2. 参数在哪里改

`simulations.yaml` **只保存原始信息**：模拟编号、条件、拓扑和轨迹路径。`data_root` 相对 YAML 解析，文件路径相对 `data_root` 解析；可用绝对路径或环境变量。

```yaml
data_root: ..
simulations:
  - id: A18_T260
    variant: A18
    temperature_K: 260
    topology: tau_k18_129aa_A18_T260_5/top.pdb
    trajectory: tau_k18_129aa_A18_T260_5/tau_k18_129aa_A18_T260_5.dcd
    # time_step_ps: 100.0
```

`time_step_ps` 是**原始相邻保存帧**的时间间隔，不是积分步长；可选 `time_origin_ps` 默认 0。缺少时间间隔时不生成实际时间。每个独立模拟使用不同 `id`。

修改各入口顶部的 `PARAMETERS`：

| 入口 | 默认参数 |
|---|---|
| `01_center/run.py` | `start=750`、`stop=None`、`step=3`、`expected_frames=750`、`max_bond_nm=1.0` |
| `02_calculate/phase_boundary/run.py` | `bin_width_nm=0.1`、`peak_halfwidth_bins=5`、`dense_coefficient=0.6`、`min_peak_concentration_mM=1.0` |
| `02_calculate/intramolecular/run.py` | `contact_cutoff_nm=0.7`、`exclude_sequence_separation=2` |
| `02_calculate/intermolecular/run.py` | `contact_cutoff_nm=0.7`、`partner_thresholds=[1,5,10]`、`network_thresholds=[1,5,10,20,30,35]` |

- 帧号从 0 开始，`stop` 不包含该帧；`None` 表示到结尾。
- `expected_frames=None` 可取消固定帧数要求；`max_bond_nm` 是整链重建检查，与接触距离无关。
- 两类接触默认 `block_frames=50`（统计块）；第二步各模块默认 `chunk_frames=25`（读取块，只影响内存）。
- 第一、二步总入口的 `SIMULATIONS=None` 表示全部；命令行 `--simulation` 可覆盖选择。单组取帧差异可在第一步的 `SIMULATION_PARAMETERS` 中设置。

## 3. 如何运行

在项目根目录、同一个 Python 环境中执行。需要 Python 3.10+：

```bash
python -m pip install -r requirements.txt

# 建议先完整跑一组；第一步成功后再执行第二步
python -u 01_center/run.py --simulation A18_T260
python -u 02_calculate/run.py --simulation A18_T260

# 全部模拟：去掉 --simulation
python -u 01_center/run.py
python -u 02_calculate/run.py
```

本机可用 `/opt/anaconda3/bin/python` 替换 `python`。项目位于 `/Users/minato/Desktop/slab_idr/contact_partner_network_analysis`。

选择多组可重复写 `--simulation`。第二步总入口总会先重算边界；只补算接触时使用独立入口：

```bash
# 已有第一步结果后，计算边界
python 02_calculate/phase_boundary/run.py --simulation A18_T260
# 已有有效边界后，独立计算两类接触
python 02_calculate/intramolecular/run.py --simulation A18_T260
python 02_calculate/intermolecular/run.py --simulation A18_T260
```

总入口也支持 `--analysis phase_boundary`、`--analysis intramolecular` 或 `--analysis intermolecular`。入口参数帮助见 `--help`。

**每次运行都会覆盖所选模拟在本模块的结果，不自动复用。** 修改居中/取帧后重跑全部下游；修改边界后重跑两类接触；修改某类接触参数只重跑该类；修改绘图样式只重画。

### 画图

```bash
python -m jupyterlab
```

打开各模块的 `notebooks/`，从头运行。只算了一组时将 `SELECTED` 改为 `["A18_T260"]`，全部完成后可用 `None`。Notebook 只读取数值结果；图片写到对应 `output/figures/`，并记录来源。重算后重启内核再画图。

### 服务器

从项目根目录提交，两次任务按顺序运行：

```bash
qsub -v STAGE=01_center,CONTACT_PYTHON=/path/to/env/bin/python jobs/pbs.sh
# 确认第一步成功后提交
qsub -v STAGE=02_calculate,CONTACT_PYTHON=/path/to/env/bin/python jobs/pbs.sh
```

Grid Engine 使用 `jobs/grid_engine.sh`。单组可加 `SIMULATION=A18_T260`。按集群补充队列和资源；模板未在服务器实测，不会自动等待上游任务。多组默认串行；第一步临时盘需容纳一份居中 DCD。

## 4. 每一步算什么

### 第一步：居中

1. 按 `start/stop/step` 选原始帧。
2. 沿 z 以不超过 1 Å 的均匀箱统计珠密度，将最高箱移到盒中心。
3. 找珠数最多的连续非空区域（合并周期首尾），再将其加权中心移到盒中心。
4. 沿相邻珠的周期最短位移恢复完整链，再整链平移，使链几何中心位于主盒内。

只负责居中和整链重建，不判断轨迹是否平衡或是否相分离。第一步使用 Å，第二步使用 nm；1 nm = 10 Å。原始输入应是生产轨迹，不要对旧版已截帧的轨迹重复处理。

### 第二步上游：密度与浓相核心

- 把珠按周期位置放回盒内，沿 z 每 **0.1 nm** 分箱。
- 每箱珠数除以链长和箱体积，换算成链当量浓度 mM，再对所有选中帧等权平均。它不是质量密度，也不是链中心直方图。
- 中心箱及左右各5箱，共 **11箱、约1.1 nm** 的平均密度作为峰高 h。
- 从中心向右找首个密度严格低于 h/2 的箱，其中心位置定义半高半宽 w；浓相核心半宽为 **0.6 × w**。不拟合、不插值，也不是找峰高60%的位置。
- 每帧按整链几何中心分类：相对盒中心的周期 z 位移满足 **|Δz| < 0.6w**，该链属于浓相。左侧穿越位置仅供查看，不参与边界计算。

中心峰高不足1 mM或没有有效穿越时失败，**不回退到固定宽度**。当前不定义稀相，“非浓相”不等于稀相。

### 分子内：链内接触如何组织

同链残基珠的周期最小镜像距离严格 **<0.7 nm**，且序列间隔 **>2**，记为接触。按链内顺序索引，每对 i–j 只计一次。

| 指标 | 定义 |
|---|---|
| 接触数 E | 每条链的唯一接触对数 |
| 度分布 P(k) | 接触了 k 个其他残基的残基比例，包含 k=0；平均度为 2E/N |
| 最大接触簇占比 | 最大连通分量的残基数 / 链长 N |
| 全局效率（辅助） | 对不同残基对平均最短接触路径长度的倒数；不连通贡献0 |
| 接触概率图 | 同一序列位置对接触的帧内链平均，再对有效帧平均 |

每帧每条链单独建无向无权网络，包含全部残基，**不加入主链边**。无接触时 E=0、P(0)=1、效率=0、最大簇占比=1/N。分别统计全体系和浓相。

这些指标描述接触组织，不代表稳定折叠或物理传递效率。序列弧线图显示瞬时接触位置；时间平均概率图不能当作某个真实瞬时网络。

### 分子间：接触量、伙伴广度与网络

不同链 a、b 之间，距离严格 **<0.7 nm** 的残基对数记为 Wab。底层保存全部链对；下列统计只使用**双方均属浓相**的链对。

| 指标 | 定义 |
|---|---|
| 每链总接触数 Ca | 与其他浓相链的 Wab 之和 |
| 伙伴数 Ka(m) | 满足 Wab ≥ m 的不同链数，m默认1、5、10 |
| 每伙伴接触比值 | Ca/Ka(m)，零伙伴时为 NaN |
| 伙伴数分布 | 每帧各伙伴数的链比例，再对有效帧平均，包含零伙伴 |
| 最大连通分量比例 | 最大接触簇链数 / 全部浓相链数，包含孤立链 |
| 全局效率 | 对不同链对平均最短网络路径长度的倒数；不连通贡献0 |

网络节点是链，Wab ≥ m 时连接无向无权边，网络阈值为 **1、5、10、20、30、35**。路径长度计边数，不是空间距离。无浓相节点时两项网络指标为 NaN；只有一条链时效率为 NaN、最大分量比例为1。

接触比值的分子包含全部浓相接触，未筛掉低于伙伴阈值的接触。例如与三条链有12、3、1对接触，阈值5时伙伴数为1，比值是16而非12。体系所有链的 Ca 求和会把链间接触计两次。

## 5. 平均方式与输出

**逐链计算 → 帧内所选链等权平均 → 有效帧等权平均。** 网络指标不从跨帧合并的边计算；浓相链多的帧不额外加权。度分布先在每条链内归一化，再按上述顺序平均。

- 0：有样本但无接触；NaN：无浓相链、无有效比值或样本不足。
- `*_mean_frames`、`*_std_frames`、`*_count_frames`：帧均值的平均、样本SD、有效帧数。
- 每50个分析帧为一块，保存块均值及块间SD、`SEM=块间SD/√有效块数`；少于两块时SD/SEM为NaN。总体均值直接从帧计算，不受末块长度影响。
- 时间块不是独立模拟重复，固定50帧不保证独立。分子内节点是残基，分子间节点是链，适合比较变化趋势，不直接比较网络指标绝对大小。

各模块输出到自己的 `output/<模拟编号>/`：

| 模块 | 主要文件 |
|---|---|
| 居中 | `centered.dcd`、`topology.pdb`、`chains.json`、`frame_index.csv`、居中平移记录 |
| 边界 | `density_profile.csv`、`boundary.json` |
| 分子内 | `chain_contacts.csv.gz`、`contact_map.npz`、`degree_distribution/`、`edges/` |
| 分子间 | `chain_partners.csv.gz`、伙伴分布表、`network_frames.csv`、`network/`、`edges/` |

两类接触另有 `frame_means.csv`、`summary.csv`、`block_means.csv`。各模块保存 `metadata.json`（参数及上游来源）和成功后的 `COMPLETE`。帧号和残基索引从0开始，原始帧号及时间来自第一步。

## 6. 覆盖与失败

开始计算时清空本模块该模拟的旧输出；成功后才生成 `COMPLETE`。失败/中断的目录不可作为完成结果，重新运行对应模块即可。上游变化后，下游和图片需手动按顺序重算。

批量结束打印 `SUCCESS`、`FAILED`、`SKIPPED` 汇总；失败日志包含模拟编号、模块和原因。边界失败跳过两类接触；分子内外彼此失败不互相阻断。任一失败返回非零退出码。qsub 提交成功不等于计算成功，需查看作业日志。

不提供锁、自动复用、历史代次或自动依赖调度。不要同时读写同模块的同一模拟目录。旧 `centered/`、`results/`、`figures/` 保留，不能直接作为新版产物。

## 7. 逐步替代 ana_v4

### 已覆盖和暂缺内容

| 内容 | 当前情况 |
|---|---|
| 居中、平均密度、浓相核心边界 | 已有；无固定宽度回退 |
| 分子内接触、分子间接触 | 已有；扩展残基网络、伙伴数和链网络 |
| 稀相边界及浓/稀/过渡分层 | 尚缺，目前为全体系/浓相 |
| 浓相浓度、稀相浓度/csat、contrast、CV、浓度比对数dG及误差 | 尚缺，密度剖面不等于这些标量已计算 |
| Rg、序列κ/SCD | 尚缺 |
| 软接触、电荷分类、分子间残基位置矩阵、稀相代表链接触 | 尚缺 |
| 构象PCA、密度时间演化、跨指标比较等旧图 | 按需要迁移计算与绘图 |
| 旧版完整性/收敛诊断 | 未整体迁移，不属于第一步居中的职责 |

### 替代策略

保留旧代码和旧结果作为参考。先跑通并固定当前接触分析，再按科研问题逐个增加同级模块：明确输入、定义和输出，验证后再替代旧版对应功能。暂不追求统一总表，也不整体搬入旧版总控或注册器。序列指标不依赖轨迹，无需强行依赖居中；正式计算放模块，Notebook只画图。

### 新旧结果对照注意点

以下为旧版默认/实现，实际运行可能被配置覆盖：

| 项目 | 差异 |
|---|---|
| 取帧 | 旧版默认从轨迹四分之一开始，新版固定从750开始；核对真实原始帧号和抽帧间隔 |
| 接触规则 | 旧版默认1.0 nm、分子内排除间隔≤3；新版0.7 nm、排除≤2；soft不能与二元count直接比较 |
| 分子内计数 | 旧版对称矩阵求和包含i–j和j–i；新版每对一次。旧版展开坐标欧氏距离与新版周期最小镜像也可能不同 |
| 分子间对象 | 旧版常规计算每帧中心代表链与所有其他链；新版伙伴/网络统计全部浓相链对，不能要求总量相等 |
| 边界与密度 | 核对分箱、周期位置、峰高窗口和0.6系数；旧版固定宽度回退不等于测得边界 |
| 平均与误差 | 旧版部分量合并链—帧样本；新版先帧内链平均。核对零/NaN、样本权重和分块方法 |

先统一数据、帧、单位和定义，再核对结果。**同定义的差异需要排查；定义不同则记录方法变化。** 正式条件比较应统一采用同一套设置重算，不混用新旧标量或按文件排列顺序拼接。尚缺的旧版指标可暂时保留旧分析来源，不将其误标为新版结果。

当前支持固定正交盒、每segment一条有序线性链、每残基一珠、等长链。数据和输出不纳入Git。
