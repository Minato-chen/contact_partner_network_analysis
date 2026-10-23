# 分步模拟分析：接触、构象与内外关系

原始轨迹只读，计算产生数值文件，Notebook只读取结果画图。按上游和研究对象分层，各文件夹的定义见本层README。

## 结构

```text
simulations.yaml                       只含原始模拟信息
01_center/                             选帧、居中、整链重建
02_calculate/
  run.py                               自动发现直接子目录的run.py
  01_phase_boundary/                   密度与浓相核心边界
  02_conformation/
    01_rg/                             所有帧、所有链的Rg，只计算一次
  03_contacts/
    00_prepare/                        一次扫描生成公共接触与对象选择
    dependencies/                      两类对象共用计算方法
    01_central_chain/                  每帧最靠中心的浓相链
    02_dense_chains/                   每帧全部浓相链
      01_intramolecular/               接触数、序列概率、Rg及补充网络
      02_intermolecular/               接触、伙伴、网络和路径分布
      03_joint_analysis/               同帧同链的内部—外部配对
```

中心链对象下也有相同三个子目录。正式输出统一为各计算目录的 `output/<模拟编号>/`，配套绘图放在 `notebooks/`。父入口按文件夹编号运行直接子入口；新增带run.py的子文件夹自动纳入，不维护名单。需要遵守统一 `calculate(simulation, centered_root)` 接口；序号控制顺序，不自动猜测依赖。

## 参数与运行

原始文件、编号、条件在 `simulations.yaml`。路径相对 `data_root`，后者相对YAML。可选 `time_step_ps` 是原始保存帧间隔，缺少则时间留空。变异体和独立重复只新增清单记录，不复制代码。

| 参数 | 修改位置 |
|---|---|
| start/stop/step、预期帧数 | `01_center/run.py` |
| 分箱、峰高窗口、浓相系数 | `02_calculate/01_phase_boundary/run.py` |
| Rg端基质量修正、读取块 | `02_calculate/02_conformation/01_rg/run.py` |
| 接触距离、序列排除 | `02_calculate/03_contacts/00_prepare/run.py` |
| 伙伴/网络阈值、统计块 | 对象下 `02_intermolecular/run.py` |
| 分子内/联合分析统计块 | 各自run.py |
| 选择模拟、图形样式 | 命令行或Notebook开头 |

当前保留本地接触定义：距离严格 <1.0 nm，链内排除序列间隔≤2；0.7–1.0 nm最终主阈值尚未决定。伙伴阈值1、5、10；网络阈值1、5、10、20、30、35。不要把重构与科学阈值选择混为一项修改。

```bash
# 在项目根目录、同一个Python 3.10+环境执行
python -m pip install -r requirements.txt
python -u 01_center/run.py --simulation A18_T260
# 确认第一步成功后
python -u 02_calculate/run.py --simulation A18_T260
```

去掉 `--simulation` 默认处理全部清单记录，也可重复该选项选择多组。本机可将python替换为 `/Users/minato/venvs/mdtraj-contact/bin/python`。

```bash
# 只运行接触分支：要求相界和Rg已完成
python 02_calculate/03_contacts/run.py --simulation A18_T260
# 只运行浓相对象：要求公共接触和Rg已完成
python 02_calculate/03_contacts/02_dense_chains/run.py --simulation A18_T260
# 只重画：打开对应notebooks，无需运行计算
python -m jupyterlab
```

父入口可用 `--task 直接子目录名` 选择子任务；不会自动补算选中分支以外的上游。例如总入口 `--task 03_contacts` 要求边界和Rg已存在。叶入口支持 `--output-dir`，但自定义输出不会自动修改下游读取位置；正常分析使用默认路径。所有第二阶段入口支持 `--centered-dir`，须保持一致。帮助见 `--help`。

## 核心定义

- **居中**：按原始帧号选帧，密度峰先平移至盒中心，再按最大连续占据段加权中心调整，最后恢复完整链。每segment一条有序等长链，每残基一珠，固定正交盒。不承担原始轨迹平衡/质量审查。Å与nm换算为10:1。
- **相界**：z方向0.1 nm分箱，珠数除链长和箱体积得到链当量mM。对帧平均密度；中心及左右各5箱（11箱）均值为峰高h，右侧首个密度<h/2的箱中心为半宽w，浓相核心半宽=0.6w。峰高不足1 mM或无穿越时报错，不使用固定比例回退。
- **对象**：链几何中心距slab中心的周期z位移严格小于核心半宽则为浓相。中心链从这些链中选|z|最小者，并列选最小链号；无人则该帧对象缺失。中心链可换身份，不能当固定链时间轨迹。分子间伙伴均为其他浓相链。
- **Rg**：质量加权质心下的残基均方距离开根号，单位nm。质量来自CALVADOS2残基表，首尾分别+2/+16 Da，与旧版定义一致；质量及实际用值保存到Rg输出。Rg与接触阈值无关。
- **链内接触**：同链珠的周期最小镜像距离<阈值且序列间隔>2，每对只计一次。序列曲线是该残基至少参与一次有效接触的概率，不是伙伴数或矩阵行和。
- **链间接触**：Wab为两链的残基接触对数。每链总量Ca=ΣWab；伙伴数为Wab≥m的不同链数。Ca/伙伴数在零伙伴时为NaN，分子仍包含未达伙伴阈值的接触。
- **浓相网络**：所有浓相链为节点，包括孤立链，Wab≥m连无向无权边。最大簇比例=最大簇链数/浓相链数；全局效率为不同节点对最短步数倒数的平均，不连通贡献0。中心链分支不计算全局浓相网络。
- **路径分布**：每帧以全部不同浓相链对为分母，不连通单列；图显示1–10步和>10。长尾相加但不重新归一化，尾部SD由各块尾部概率重新计算。少于两节点的帧缺失。
- **联合分析**：按simulation_id、source_frame、chain_id匹配内部接触、Rg、外部接触和伙伴数。每残基口径内部为2E/N、外部为Ca/N。中心链可逐帧变化；相关不等于竞争或因果。

标量和序列概率先在每帧对所选链等权平均，再对有效帧等权平均；浓相链多的帧不额外加权。0是有样本无接触，NaN是无样本。默认50分析帧一块，保存帧均值、SD、有效数、块均值及块SEM；不足两块时误差缺失，时间块不是独立模拟重复。

## 输出、绘图与失败

公共Rg保存 `chain_rg.csv.gz` 和 `bead_masses.csv`；公共接触保存 `edges/*.npz` 与 `chain_selection.csv.gz`。对象分析保留逐链明细、帧均值、summary和块均值。序列参与概率在 `residue_probability/`；成对概率矩阵只作补充保存，不画热图。联合配对表为 `paired_chain_frames.csv.gz`。

分子内Notebook：01接触数与Rg、02序列概率曲线、03补充度分布/网络。分子间：01接触与伙伴、02伙伴分布；浓相分支另有03网络及路径分布。联合Notebook展示内部—外部—Rg的配对关系。`SELECTED=None`读取所有模拟，只算一组时设为对应编号列表。误差线为块SEM，分布阴影为块SD，不是置信区间。

各模块开始时覆盖本模块选中模拟的旧输出，全部写完后才生成 `COMPLETE`。父入口会继续独立任务；缺少有效上游或来源不一致时标为 `SKIPPED`，计算异常为 `FAILED`，结束打印汇总，任一失败或跳过返回非零退出码。子入口不偷偷跨层重算。上游改变后手动按顺序重算下游，不同时读写同模块同模拟。不使用锁、哈希复用或隐藏结果代次。

旧第二阶段目录（含用户运行的Notebook和结果）保存在本地忽略目录 `.legacy_stage2/`，不自动冒充新版结果。第一步结果保留可复用。先运行新版第二阶段，再打开新路径Notebook；旧图不会自动更新。

## 服务器和扩展

`jobs/pbs.sh` 与 `jobs/grid_engine.sh` 仍分别运行 `01_center/run.py` 和 `02_calculate/run.py`。从项目根目录先提交第一阶段，确认成功后再提交第二阶段：

```bash
qsub -v STAGE=01_center,CONTACT_PYTHON=/path/to/python jobs/pbs.sh
qsub -v STAGE=02_calculate,CONTACT_PYTHON=/path/to/python jobs/pbs.sh
```

可加 `SIMULATION=A18_T260`。队列资源按集群调整，模板未在服务器实测；qsub返回编号不代表分析成功。

相比ana_v4，尚未补齐稀相定义、csat/浓度、contrast/CV/dG、κ/SCD、电荷分类、软接触、分子间残基位置图和部分构象图。Rg现已纳入公共上游。以后csat依赖密度与稀相定义，作为同级上游模块添加，不绑定接触。

逐项迁移时先统一原始帧、单位、接触阈值、近邻排除、对象和平均方式。旧版常规分子间采用中心链对所有其他链；新版中心链仅对其他浓相链，不能直接要求总量相等。旧版链内矩阵和计双向、部分量合并链—帧样本；新版E只计一次且先帧平均。旧版固定宽度回退与soft值也不能混入新版同一组比较。
