# 分子内与分子间接触分析

研究链内部接触、与其他链的连接，以及同帧同链的内外关系。分子内、分子间各自扫描轨迹，两个对象共用算法，不再有独立prepare步骤。

```text
run.py                 按顺序运行全部子计算
01_intra/
  run.py               分子内入口
  intra.py             参数、接触数、序列概率及补充网络
02_inter/
  run.py               分子间入口
  inter.py             参数、伙伴及浓相网络统计
  path_distribution.py 最短路径分布
03_joint/
  run.py               联合入口
  joint.py             按原始帧号和链编号配对内外接触及Rg
common/                对象选择、网络和结果读取工具
```

## 输入、参数与运行

上游为01的居中轨迹和基础相界；分子内还需要基础Rg，联合分析需要上述分子内与分子间结果。参数在intra.py、inter.py、joint.py顶部，run.py只负责运行。共同距离改变时同时改intra.py与inter.py；联合分析拒绝距离不一致的结果。

```bash
python 02_calculate/02_contacts/run.py
python 02_calculate/02_contacts/01_intra/run.py --simulation B01_T260 --object dense_chains
python 02_calculate/02_contacts/02_inter/run.py --object central_chain
```

对象默认都计算。dense_chains是每帧几何中心位于上游浓相核心内的全部链；central_chain是其中最靠近slab中心的一条，并列取最小链号。中心链身份可以逐帧变化，不能当作固定链时间轨迹。分子间伙伴只包括其他浓相链。

## 核心定义

- 接触：周期最小镜像距离严格<0.8 nm。分子内排除序列间隔≤2，每对仅计一次。
- 分子内序列概率：残基至少参与一次有效内部接触的概率，不是接触矩阵行和。残基为接触网络节点，主链不计边；最大簇、效率与度分布仅为补充，不直接称为折叠度。
- 分子间：两链残基对数为Wab；每链总接触数Ca=ΣWab。Wab≥m为伙伴，m=1/5/10。Ca/伙伴数在零伙伴时缺失；m>1时分子仍包含未达m的接触，不能直接称为合格伙伴的平均接触数。
- 浓相网络：全部浓相链为节点，包括孤立链；Wab≥m连无向无权边，m=1/5/10/20/30/35。最大簇比例是最大簇链数/节点数；全局效率是不同节点对最短步数倒数均值，不连通贡献0。
- 路径分布：每帧以全部不同浓相链对为分母，不可达单列；不足两个节点时缺失。图显示1–10步与>10，尾部块SD重新计算，不对可达部分重新归一化。
- 联合：按simulation_id、source_frame、chain_id匹配；每残基内部接触为2E/N，外部为Ca/N。关联不等于竞争或因果。

标量、序列概率先对帧内所选链等权平均，再对有效帧等权平均。0是有样本无接触，NaN是无样本。默认每50分析帧一块，保存块SD和SEM；时间块不是独立模拟重复。

## 输出与绘图

03_results/<编号>/02_contacts/下按01_intra、02_inter、03_joint分类，再按central_chain/dense_chains区分对象。各自common/保存本计算的底层分块数据，不是独立任务。重扫底层数据会使未选对象的旧完成标记失效。

逐链表、帧表、分布等仍保存；每个对象另附scalars.json用于04自动汇总，数值来自原有summary.csv，不改算法。04生成contacts_central_chain.csv和contacts_dense_chains.csv；网络只出现在浓相表。详细图在05_plot/02_contacts/的Notebook中，切换OBJECT后从头运行。
