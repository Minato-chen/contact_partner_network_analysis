# 基础计算

基础指标供不同研究课题复用，一个指标一个py；参数在各文件顶部。

```text
run.py                自动发现并运行本层指标，相界优先
phase_boundary.py     时间平均密度与浓相核心边界
kappa.py              使用localCIDER计算序列电荷分布参数κ
scd.py                按序列电荷及序列间距计算SCD
rg.py                 全部链、全部分析帧的质量加权Rg
residue_masses.csv    CALVADOS2残基质量表
```

## 方法与结果

相界：z方向0.1 nm分箱，珠数除以链长换算成链等效mM。时间平均后，中心及左右各5箱（11箱）的均值为峰高；右侧首个低于半峰高的箱中心距离为w，核心半宽为0.6w。左侧穿越仅作诊断。这是半高规则，不是拟合；无有效穿越时失败，不使用固定宽度回退。

Rg：整链坐标相对质量加权质心的均方距离开根号，单位nm。每残基一珠，首尾质量分别加2/+16 Da。Rg不依赖接触阈值。基础汇总的Rg是全部链先帧内平均、再对帧等权平均；这里不提供尚未计算的块SEM，浓相/中心链Rg及误差见接触结果。

详细结果在03_results/<编号>/01_basic/{phase_boundary,rg}/，成功后有COMPLETE、metadata.json和scalars.json。后者仅提取最终汇总值，04不重复计算。

```bash
python 02_calculate/01_basic/run.py
python 02_calculate/01_basic/rg.py --simulation B01_T260
```

## κ

直接调用官方`from localcider.sequenceParameters import SequenceParameters`和`SequenceParameters(sequence).get_kappa()`，不自行实现公式。κ描述正负电荷沿序列的分布，通常越大表示同号电荷越集中；它是序列属性，不是轨迹平均量，也没有时间SEM。定义见[localCIDER文档](https://pappulab.github.io/localCIDER/)。

序列读取01输出的topology.pdb和chains.json，检查所有链序列相同后计算一次，不拼接多条链。仅接受20种标准氨基酸；混合序列或官方返回未定义κ时明确失败。结果记录序列、上游元数据和库版本；不依赖相界、Rg或接触阈值。

```bash
python 02_calculate/01_basic/kappa.py
python 02_calculate/01_basic/kappa.py --simulation B01_T260
```

结果保存到03_results/<编号>/01_basic/kappa.scalars.json。运行`python 04_summary/run.py`自动加入basic.csv。

## SCD

采用`SCD = Σ(n<m) q[m]q[n]√(m−n) / N`，每对只计算一次。D/E为−1，K/R为+1，其他残基（包括H）为0；N为完整链长。保留符号，不取绝对值、不归一化到0–1，不对序列近邻作排除。它描述电荷与序列间距的组合，不是空间接触或轨迹平均量。

与κ共用序列读取：来自01输出，所有链序列必须一致，只计算一次。输出为03_results/<编号>/01_basic/scd.scalars.json，包含序列、公式、电荷规则及上游元数据，无时间SEM。

```bash
python 02_calculate/01_basic/scd.py
python 04_summary/run.py
```

基础run自动包含SCD，汇总自动加入basic.csv，无需另行注册。

## 以后新增单值指标

csat尚未实现。类似κ的单值指标可以直接写03_results/<编号>/01_basic/<指标名>.scalars.json，不再创建目录、单值CSV或额外完成标记。示例输出接口（不包含κ算法）：

```python
from common import results as files
path = results_root / simulation["id"] / "01_basic/kappa.scalars.json"
files.start_scalar(path, simulation, parameters=PARAMETERS, inputs={"sequence": sequence})
value = ...  # 此处执行指标本身的计算
files.finish_scalar(path, [files.scalar_record("kappa", value, "dimensionless")])
```

先start再计算，失败会保持未完成，避免旧数值冒充新结果。记录注明单位与定义；不同方法可用definition字典区分。基础指标object留空，合入同一张basic.csv。

csat这类轨迹统计量应保留逐帧、时间块和区域等详细结果，再附scalars.json；不能只留下一个均值。辅助函数放公共工具目录或以下划线命名，避免被当作新指标执行。
