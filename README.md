# 分步模拟分析

本仓库研究分子内、分子间接触及两者关系。原始数据只读，计算保存数值，汇总整理指标，Notebook负责绘图。

```text
simulations.yaml   原始模拟编号、条件和文件路径
01_center/         选帧、轨迹居中及整链重建
02_calculate/
  01_basic/        基础指标：相界、Rg、κ、SCD；一个指标一个py
  02_contacts/     接触课题：分子内、分子间、联合分析
03_results/        第一层是模拟编号，内部分类与02一致
04_summary/        读取03，按课题生成每模拟一行的对比表
05_plot/           Notebook与图片，分类对应02
common/            模拟清单读取、结果文件和状态打印工具
```

## 运行

在项目根目录，使用同一个Python环境：

```bash
source ~/venvs/mdtraj-contact/bin/activate
python -u 01_center/run.py
python -u 02_calculate/run.py
python 04_summary/run.py
python -m jupyterlab
```

已有01结果时无需重跑居中。--simulation B01_T260选择单组，选项可重复；省略则处理全部清单。参数在各计算py顶部，模拟清单只放原始信息。父入口自动发现本层计算；选择--task时不会补算上游。命令详情用--help。

详细结果在03_results/<编号>/；汇总在04_summary/output/；图片在05_plot对应分类的figures/。03数值、04输出和导出图片不进入Git，需要另行备份。05可以读取03画详细图，也可以读取04比较汇总指标。

## 定义与使用说明

- [居中](01_center/README.md)：选帧、单位与输出。
- [基础指标](02_calculate/01_basic/README.md)：相界、Rg、κ、SCD，以及单值指标的输出约定。
- [接触分析](02_calculate/02_contacts/README.md)：子结构、对象、接触/伙伴/网络定义及独立运行方法。
- [汇总](04_summary/README.md)：只运行run.py即可更新各课题表，不维护指标清单。
- [绘图](05_plot/README.md)：Notebook与对象选择。

计算成功后写COMPLETE和metadata.json；失败或上游缺失会打印FAILED/SKIPPED，返回非零。上游改变后按依赖顺序重算下游；不要同时读写同一计算结果。

## 服务器与后续扩充

jobs/pbs.sh、grid_engine.sh分别支持STAGE=01_center或02_calculate；CONTACT_PYTHON指定服务器Python，SIMULATION可选。先确认01成功再运行02。04汇总等计算结束后单独执行，不在并行计算中共同写总表。

csat尚未实现；未来新增基础指标文件即可，不提供空算法。新课题与02_contacts并列，在03采用同名分类并输出汇总记录，04自动产生新表。cutoff_comparison是可独立删除的四距离对照，正式接触定义目前严格<0.8 nm。

相比ana_v4，还未补齐csat/稀相浓度、contrast/CV/dG、电荷分类、软接触和部分图形。逐项替代前统一原始帧、单位、距离、序列排除、对象与平均方式。旧版中心链常对所有链，新版仅对其他浓相链；旧版链内矩阵可能双向计数，新版每对一次，不能直接要求总量相等。
