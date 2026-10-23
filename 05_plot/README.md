# 绘图

Notebook只读数值结果，不承担轨迹计算。顶部SELECTED选择模拟，OBJECT选择central_chain或dense_chains；切换后从头运行。

```text
01_basic/
  01_density_and_boundary.ipynb  密度、半高位置和浓相核心
  rg.ipynb                      全部链Rg
02_contacts/
  01_intra/                     接触数与Rg、序列概率、补充网络
  02_inter/                     接触与伙伴、伙伴分布、浓相网络与路径
  03_joint/                     同帧同链的内外接触与Rg
```

详细图直接读取03_results；跨指标对比也可读取04_summary/output中的宽表。现有Notebook保持直接读取03的方式，不强制先运行04。改用自定义数值根目录时修改RESULTS_ROOT。

图片保存到本分类figures/<对象>/；基础图在01_basic/figures/。PNG/PDF/SVG及来源记录可重新导出，不进入Git。误差线为块SEM，分布阴影为块SD，不是独立重复实验的置信区间。网络快照位置是拓扑布局，不是物理坐标。
