# 可视化与学术论文绘图模块 (`05_plot`)

本模块作为整个分析流水线的**成果展示系统（Stage-5）**，专为发表国际顶尖学术期刊（Nature 系列、PNAS）高标准插图设计。模块深度融合了“**故事线驱动（Storyline-Driven）**”与“**两阶段多层级物理展开**”架构，支持从全局宏观相图到微观残基接触的全景证据链可视化。

> **严格输出规范**：根据论文准备与课题规范要求，**全部图表严格仅导出 300 DPI 高分辨率 `.png` 格式**（无损保存，不输出 `.pdf`、`.svg` 等冗余格式）。

---

## 一、 课题驱动目录树与交互组件

```text
05_plot/
├── common/                                      # 公共样式与数据加载底层支持
│   ├── style.py                                 # Nature 标准排版规范、调色板与严格 PNG 导出工具
│   └── data_loader.py                           # 桥接 04_summary 汇总表与 03_results 明细时序
├── 00_cross_topic_overview/                     # 跨课题全 7 变体全景总览
│   ├── 01_cross_topic_overview.ipynb            # 全景总览交互式 Jupyter Notebook
│   ├── plot_cross_topic.py                      # 全景图生成脚本
│   └── figures/{dense_chains,central_chain}/    # 严格 300 DPI PNG 图表导出目录
│       ├── cross_topic_phase_behavior.png       # 平板密度剖面 ρ(z) 与稀相平衡浓度 Csat
│       ├── cross_topic_multivalent_tradeoff_map.png # 多价相分离权衡图 (有效伙伴数 vs 单伙伴最大异号接触数)
│       └── cross_topic_descriptors_matrix.png   # 结构、接触、伙伴与网络六联指标对比矩阵
├── 01_topic1_blockiness/                        # 课题一：电荷块大小效应 (A18, B01, B10)
│   ├── 01_intra.ipynb / plot_topic1_intra.py    # 分子内构象紧实度、残基接触谱与环尺寸分布
│   ├── 02_inter.ipynb / plot_topic1_inter.py    # 分子间电荷拆解、伙伴数均值与完整 PMF 分布、网络逾渗
│   ├── 03_joint.ipynb / plot_topic1_joint.py    # 分子内外接触平衡与 2D KDE 多价景观态
│   └── figures/{dense_chains,central_chain}/    # 课题一纯 PNG 图表输出
├── 02_topic2_polarization/                      # 课题二：电荷极化与中性尾巴屏蔽 (A20, C04, B09, C22)
│   ├── 01_intra.ipynb / plot_topic2_intra.py    # 分子内构象、中性尾巴暴露与局部电荷塌缩
│   ├── 02_inter.ipynb / plot_topic2_inter.py    # 偶极哑铃陷阱、尾巴位阻排斥、伙伴数 PMF 与网络瓦解
│   ├── 03_joint.ipynb / plot_topic2_joint.py    # 分子内外平衡与 2D KDE 多价景观态
│   └── figures/{dense_chains,central_chain}/    # 课题二纯 PNG 图表输出
├── run.py                                       # 顶层一键批处理导出全套图表入口
└── PLOT_README.md                               # 本技术文档
```

---

## 二、 课题一故事线与图表系统 (Topic 1: Charge Blockiness)

### 1. 核心科学问题与反直觉现象
- **研究体系**：
  - **B01**（细交替，$\kappa=0.076$，$C_{\mathrm{sat}}=0.145$ mM）：超强相分离能力；
  - **A18**（野生型基准，$\kappa=0.228$，$C_{\mathrm{sat}}=0.949$ mM）：中等相分离能力；
  - **B10**（大电荷块，$\kappa=0.518$，$C_{\mathrm{sat}}=1.305$ mM）：弱相分离（$C_{\mathrm{sat}}$ 减弱近 10 倍）。
- **核心矛盾**：常识直觉认为“大电荷块局部吸引力更强，理应更容易凝聚”，但模拟显示相分离能力随着电荷块增大急剧下降。微观物理机制是什么？

### 2. 逻辑递进链条与图表对应
```text
[起] 分子内预构象消耗
     ├── B10 链内异号块强折叠塌缩，Rg 骤减至 1.84 nm，长程接触受抑；B01 开放伸展 (Rg ≈ 2.45 nm)
     └── 图表：topic1_intra_metrics.png & topic1_residue_and_loop_profile.png
     ↓
[承] 分子间局部过黏形成“静电拉链锁定 (Zipper Lock)”
     ├── B10 异号吸引集中在单一伙伴链上 (单伙伴最大异号接触飙升至 20.8 对，是 B01 的 3.7 倍)
     └── 图表：topic1_inter_charge_decomposition.png & topic1_inter_partners_and_max_adhesion.png
     ↓
[转] 多价配位数与概率分布向左塌缩
     ├── 强拉链消耗了多价性，B10 有效伙伴数跌至 8.13 (B01 为 10.61)
     └── 图表：topic1_inter_partner_distribution_pmf.png (展示 threshold 1, 5, 10 时的完整概率密度分布)
     ↓
[合] 超分子网络解体与多价自由能景观
     ├── B01 保持 98% 巨组分逾渗与高全局效率；B10 在强相互作用阈值下网络全面瓦解成孤岛
     ├── 2D KDE 景观清晰印证：B01 表现为大 Rg、广域分布的流体态；B10 跌入小 Rg、高单体接触的锁定态
     └── 图表：topic1_network_percolation.png & topic1_multivalent_landscape.png
```

---

## 三、 课题二故事线与图表系统 (Topic 2: Polarization & Neutral Tail)

### 1. 核心科学问题与双重破坏机制
- **研究体系**：
  - **A20**（对称两端交替，$\kappa=0.088$，$C_{\mathrm{sat}}=0.125$ mM）：超强相分离基准；
  - **B09**（单端非对称团簇，$\kappa=0.354$，$C_{\mathrm{sat}}=1.212$ mM）：适度相分离；
  - **C04**（两极电荷哑铃，$\kappa=0.887$，$C_{\mathrm{sat}}=4.142$ mM）：极弱相分离；
  - **C22**（密集电荷核心 + 60 残基完全中性尾巴，$\kappa=0.407$，$C_{\mathrm{sat}}=4.584$ mM）：极端弱相分离。
- **核心矛盾**：C04 与 C22 同样表现出极端微弱的相分离能力（$C_{\mathrm{sat}} > 4.1$ mM），但它们背后的分子机制截然不同：一个是“**偶极哑铃双体陷阱**”，另一个是“**内部自中和折叠 + 中性尾巴空间位阻排斥**”。

### 2. 逻辑递进链条与图表对应
```text
[起] 分子内序列区域分化
     ├── C22 前 60 残基完全不参与接触 (概率 ≈ 0)，而后半段电荷核心极端致密 (分子内接触高达 172.9 对/链)
     ├── C04 呈现显著的头尾相接极性环化
     └── 图表：topic2_intra_metrics.png & topic2_residue_and_loop_profile.png
     ↓
[承] 分子间微观破坏的分野：双体陷阱 vs 空间屏蔽
     ├── C04 异号占比高达 35.2%，单伙伴最大异号接触达 18.3，深陷反平行二聚体/少体寡聚陷阱
     ├── C22 分子间总接触仅 34.3 对 (全变体最低)，中性尾巴在外围产生强烈的空间排斥与溶剂化阻力
     └── 图表：topic2_inter_charge_decomposition.png & topic2_inter_partners_and_max_adhesion.png
     ↓
[转] 多价配位数与 PMF 分布异化
     ├── C04 伙伴数跌至 4.85；C22 伙伴数跌至极低极限 2.47，PMF 曲线显著聚集在 0～2 伙伴甚至孤立链
     └── 图表：topic2_inter_partner_distribution_pmf.png (清晰印证 C22 与 A20 的鲜明对比)
     ↓
[合] 超分子网络完全不逾渗与景观凝固
     ├── C22 网络巨组分比例仅 0.05 (根本无法形成物理逾渗网络，拓扑学上宣告相分离死亡)
     ├── 联合接触占比证实 C22 的 83.5% 接触被链内自我消耗；2D 景观直观再现构象态空间的极端偏移
     └── 图表：topic2_network_percolation.png & topic2_multivalent_landscape.png
```

---

## 四、 跨课题全景机理总图 (`00_cross_topic_overview`)

1. **相行为全景 (`cross_topic_phase_behavior.png`)**：
   - 面板 A：全 7 变体在平板几何下的密度对称剖面 $\rho(z)$；
   - 面板 B：稀相饱和浓度 $C_{\mathrm{sat}}$ 柱状对比（带两课题视觉分区与块 SEM）。
2. **多价权衡与相分离机理总图 (`cross_topic_multivalent_tradeoff_map.png`)**：
   - **横轴**：单伙伴最大异号接触数（表征强相互作用局域黏附上限）；
   - **纵轴**：有效相互作用伙伴数（表征链在浓相中的多向配位度）；
   - **气泡注记**：标明各变体代号及其实测 $C_{\mathrm{sat}}$，并配备严格防重叠智能正交偏移；
   - **三个物理分区**：
     - *流体多价逾渗区 (左上)*：高伙伴数 + 适度黏附（B01, A20），形成低 $C_{\mathrm{sat}}$ 超强相分离；
     - *静电拉链锁定区 (右中)*：极端局部强黏，消耗多价性（B10, C04, B09），产生异常升高的 $C_{\mathrm{sat}}$；
     - *自中和与尾巴屏蔽区 (左下)*：内部自锁 + 空间位阻排斥（C22），形成最高 $C_{\mathrm{sat}}$。
3. **综合指标对比矩阵 (`cross_topic_descriptors_matrix.png`)**：
   - 2×3 六联面板对比 $R_g$、分子间总接触、分子间异号接触、有效伙伴数、最大异号接触与网络巨组分比例。

---

## 五、 排版防重叠与留白工程规范

针对此前排版中的截断与重叠问题，本模块已全面实施如下工程规范：
1. **防止标题削切截断**：
   - 所有多面板图均采用 `plt.tight_layout(rect=[0, 0, 1, 0.92])`，大标题统一置于 `y=0.98`，预留 8% 的顶部呼吸空间；
2. **消除柱顶数值与误差棒重叠**：
   - 柱顶数值标注点统一算法：`y = bar.get_height() + sem + offset`，严格置于误差棒帽子上方；
   - Y 轴上限显式设定为 `(max + sem) * 1.30`，彻底杜绝文字溢出外边框；
3. **散点图智能正交避让（Anti-Collision）**：
   - 针对坐标接近的点（B10 与 B09，B01 与 A20），赋予正负交替的独立偏移向量（如 B10 偏右上，B09 偏左下），杜绝文字与背景提示块相互重叠。

---

## 六、 运行指引与使用方法

### 1. 命令行一键生成全量纯 PNG 图表

在项目根目录或 `05_plot` 目录下运行：

```bash
# 自动生成全课题、全对象 (dense_chains 与 central_chain) 对应的全部 36 张 300 DPI PNG 图表
python 05_plot/run.py

# 仅针对浓相全系综 (dense_chains) 导出
python 05_plot/run.py --object dense_chains

# 仅针对中心单链探针 (central_chain) 导出
python 05_plot/run.py --object central_chain
```

### 2. Jupyter Notebook 交互式挑选与探索

在 JupyterLab 中打开任意 Notebook：
- `05_plot/01_topic1_blockiness/01_intra.ipynb`、`02_inter.ipynb`、`03_joint.ipynb`
- `05_plot/02_topic2_polarization/01_intra.ipynb`、`02_inter.ipynb`、`03_joint.ipynb`
- `05_plot/00_cross_topic_overview/01_cross_topic_overview.ipynb`

执行后图表将直接在页面内高清交互渲染，并自动以 300 DPI PNG 同步落盘至对应 `figures/` 目录中。
