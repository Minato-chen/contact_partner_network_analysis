# 分子内、分子间与联合接触分析

本模块系统化分析液-液相分离（LLPS）浓相环境中无序蛋白质（IDP）的构象特征与相互作用模式，分为三大核心支柱：
1. **分子内分析 (`01_intra`)**：单链内部紧实度、残基折叠度与内部网络。
2. **分子间分析 (`02_inter`)**：三层级联物理架构（接触数 $\to$ 伙伴数 $\to$ 网络），解析相互作用的物理化学本质、多价分配与浓相介观网络。
3. **联合分析 (`03_joint`)**：同帧、同链、逐分子对齐分子内折叠与分子间交联，揭示能量竞争与多价平衡。

---

## 统一物理定义与考察对象

- **接触距离截断**：严格定义为空间周期性最小镜像距离（Minimum Image Distance）$< 0.8\text{ nm}$（8 Å），所有子模块完全统一。
- **两套标准考察对象**：
  - `central_chain`：每帧在浓相中几何中心最靠近 slab 中心的单链探针（若有并列取最小链编号）。作为探针链，不受相界面不对称形变影响。
  - `dense_chains`：每帧质心位于上游相边界内部的全部浓相链系综平均，反映浓相本体的整体统计热力学行为。
- **分子间有效伙伴**：限定为同一时间帧内的其他浓相链（不包含稀相孤立链）。

---

## 一、 分子内分析 (`01_intra`)

### 1. 整体物理结构
分子内分析关注单链高分子构象的空间压缩程度与序列残基间的接触网络。为排除相邻肽键共价约束的假象，严格排除序列间隔 $|i - j| \le 2$ 的残基对。

### 2. 严谨指标命名与物理定义

| 指标命名 | 英文标量 Key | 单位 | 物理定义与科学含义 |
| :--- | :--- | :--- | :--- |
| **分子内接触数** | `intra.contact_count` | `pairs/chain` | 链内部满足 $|i - j| > 2$ 且空间距离 $< 0.8\text{ nm}$ 的残基接触对总数。量化单链内部非共价折叠程度。 |
| **回转半径** | `intra.rg_nm` | `nm` | 单链空间均方根分布半径 $R_g$。反映高分子链的整体尺寸与空间扩展度。 |
| **长程接触占比** | `intra.long_range_fraction` | 无量纲 (`dimensionless`) | 序列间隔 $|i - j| \ge 12$ 的分子内接触数占总分子内接触数的比例。量化大环形成与深层长程塌缩程度。 |
| **残基网络全局效率** | `intra.global_efficiency` | 无量纲 (`dimensionless`) | 将残基作为节点、有效接触作为边构建的网络中，所有节点对最短拓扑路径倒数的均值。反映链内相互作用的连通紧密度。 |
| **最大残基簇比例** | `intra.largest_component_fraction` | 无量纲 (`dimensionless`) | 最大连通残基子图所含残基数占链全长（129残基）的比例。反映是否存在核心折叠簇。 |
| **残基接触概率分布** | `contact_probability.csv` | 概率 (`fraction`) | 链上各个残基位点在系综中至少参与一次分子内接触的概率。表征沿序列分布的内部作用热点区。 |

---

## 二、 分子间分析 (`02_inter`)

### 1. 整体三层物理架构
分子间分析采用由微观到介观的清晰三层递进结构：
$$\textbf{接触数 (Contacts)} \;\longrightarrow\; \textbf{伙伴数 (Partners)} \;\longrightarrow\; \textbf{网络 (Network)}$$
摒弃口语化比喻，全面采用客观严谨的物理化学术语。

### 2. 严谨指标命名与物理定义

#### 第一层：接触数（Contacts）—— 相互作用化学本质
直接回答：**“链在外界总共发生了多少接触？这些接触分别属于什么化学性质？”**

| 指标命名 | 英文标量 Key | 单位 | 物理定义与科学含义 |
| :--- | :--- | :--- | :--- |
| **总接触数** | `inter.total_contacts` | `pairs/chain` | 单链与浓相中所有其他链形成的分子间残基接触对总数。量化对外界环境的总吸附粘附量。 |
| **异号接触数** | `inter.attract_contacts` | `pairs/chain` | 发生接触的残基对中，电荷为一正一负（`KR` 与 `DE`）的接触数。代表驱动结合的核心**静电吸引作用**。 |
| **同号接触数** | `inter.repulse_contacts` | `pairs/chain` | 发生接触的残基对中，电荷均为正（`KR-KR`）或均为负（`DE-DE`）的接触数。代表**静电排斥/挫折作用**。 |
| **中性接触数** | `inter.neutral_contacts` | `pairs/chain` | 发生接触的残基对中，至少含有一个不带电残基的接触数。代表**疏水、极性及空间排阻背景**。 |
| **异号接触占比** | `inter.attract_fraction` | 无量纲 (`dimensionless`) | $\frac{\text{异号接触数}}{\text{总接触数}}$。量化分子间吸引作用在所有接触中所占的主导程度。 |

> **数学守恒关系**：$\text{总接触数} \equiv \text{异号接触数} + \text{同号接触数} + \text{中性接触数}$。各分类完全穷尽且互斥。

#### 第二层：伙伴数（Partners）—— 多价配位数与异号作用分配
直接回答：**“分子间的接触是分摊给多个伙伴，还是集中在少数伙伴身上？与单一伙伴的异号吸引强度有多大？”**
（采用标准接触梯度阈值考察：$k = 1, 5, 10$ 对接触，量化瞬态与强结合伙伴）

| 指标命名 | 英文标量 Key | 单位 | 物理定义与科学含义 |
| :--- | :--- | :--- | :--- |
| **伙伴数** | `inter.partner_count` `[threshold=k]` | `chains` | 与目标链至少形成 $k$ 对接触的独立浓相伙伴链数量。表征高分子多价物理配位数。 |
| **单伙伴平均接触数** | `inter.contacts_per_partner` `[threshold=k]` | `pairs/partner` | $\frac{\text{总接触数}}{\text{伙伴数}}$。反映平摊到每个有效伙伴上的平均接触密集度。 |
| **单伙伴平均异号接触数** | `inter.attract_contacts_per_partner` `[threshold=k]` | `pairs/partner` | $\frac{\text{异号接触数}}{\text{伙伴数}}$。反映向各个伙伴分摊静电吸引力时的平均强度。 |
| **单伙伴最大异号接触数** | `inter.max_attract_contacts_per_partner` `[threshold=k]` | `pairs/partner` | 目标链在当期所有合格伙伴中，**与单一伙伴所形成的最大异号接触数**（峰值异号接触对数）。客观量化是否存在局部的强二体嵌合锁死，还是保持均匀多价流动。 |
| **伙伴数概率分布** | `partner_distribution.csv` | 概率 (`fraction`) | 具有 $m$ 个伙伴链的系综概率分布，反映多价配位数的起伏涨落。 |

#### 第三层：网络（Network）—— 浓相介观拓扑连通度
直接回答：**“整个浓相凝聚体是否完全连通？大分子拓扑网络传递相互作用的效率如何？”**
（采用网络连通阈值：$k = 1, 5, 10, 20, 30, 35$ 对接触，探查从渗流骨架到高强核心网络的逾渗转变）

| 指标命名 | 英文标量 Key | 单位 | 物理定义与科学含义 |
| :--- | :--- | :--- | :--- |
| **全局效率** | `network.global_efficiency` `[threshold=k]` | 无量纲 (`dimensionless`) | 以浓相链为节点、两链间接触数 $\ge k$ 为边构成的网络中，所有不同链对间最短拓扑路径倒数的系综平均。表征凝聚体中应力或信息的全局传递效能。 |
| **巨组分比例** | `network.largest_component_fraction` `[threshold=k]` | 无量纲 (`dimensionless`) | 最大连通分支包含的链数占浓相总链数的比例。接近 1 表示凝聚体在拓扑上完全逾渗为单一宏观水凝胶网络。 |
| **最短路径分布** | `path_distribution.csv` | 概率 (`fraction`) | 浓相中任意两链间最短跃迁步数（1 到 10 步及 $>10$ 步）的概率分布；不可达链对单独统计。 |

---

## 三、 联合分析 (`03_joint`)

### 1. 整体物理结构
联合分析将相同模拟、相同时间帧、相同单链的分子内指标与分子间指标逐一配对（Strict 1:1 Matching）。旨在揭示单链在“向内自塌缩”与“向外交联结合”之间的多价竞争与热力学分配平衡。

### 2. 严谨指标命名与物理定义

| 指标命名 | 英文标量 Key | 单位 | 物理定义与科学含义 |
| :--- | :--- | :--- | :--- |
| **单残基分子内接触数** | `joint.intra_contacts_per_residue` | `contacts/residue` | $\frac{2 \times \text{intra.contact\_count}}{N_{\mathrm{residues}}}$。每残基平均参与的内部接触残基端数（单接触贡献2端）。 |
| **单残基分子间接触数** | `joint.inter_contacts_per_residue` | `contacts/residue` | $\frac{\text{inter.total\_contacts}}{N_{\mathrm{residues}}}$。每残基平均参与的外部接触对数。 |
| **分子间接触占比** | `joint.inter_fraction` | 无量纲 (`dimensionless`) | $\frac{\text{inter.total\_contacts}}{2 \times \text{intra.contact\_count} + \text{inter.total\_contacts}}$。总接触中外向连接所占的比例，量化链向外拉手的热力学倾向。 |
| **分子间与分子内接触比值** | `joint.inter_intra_ratio` | 无量纲 (`dimensionless`) | $\frac{\text{inter.total\_contacts}}{2 \times \text{intra.contact\_count}}$。内外接触的直接比例关系。 |
| **多价景观态图** | Notebook/绘图输出 | 联合概率密度分布 | $R_g$ 与分子间接触数、单残基分子内接触与单残基分子间接触的双参数联合态密度图，描绘多价驱动力在相变中的相空间轨迹。 |

---

## 四、 统计不确定度与数据流转

1. **误差估计（Block-Averaging SEM）**：
   - 默认每 50 分析帧为一个独立的时间块（Block）。
   - 计算各时间块内部均值，再以块均值的标准差除以 $\sqrt{N_{\mathrm{blocks}}}$ 得到标准误（SEM），有效消除分子动力学轨迹的时间自相关。
2. **标量与汇总集成 (`04_summary`)**：
   - 每个对象的计算产物均生成标准化的 `scalars.json`。
   - `04_summary/run.py` 自动检验跨体系的一致性定义，汇总输出 `contacts_central_chain.csv` 与 `contacts_dense_chains.csv`，作为全课题对比分析与论文绘图的标准数据源。
