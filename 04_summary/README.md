# 汇总

计算完成后，在项目根目录运行：

```bash
python 04_summary/run.py
```

自动读取03_results的最终指标，不扫描轨迹、不运行科学计算。输出在本目录output：

- basic.csv：基础指标，每个模拟一行。
- contacts_central_chain.csv、contacts_dense_chains.csv：接触指标，按对象分表。
- status.csv：来源是否成功、缺失、未完成或过期。
- sources.json：实际来源、定义与元数据，便于追溯。

列名注明指标、条件（例如threshold=5）和单位；有块SEM的指标额外生成SEM列。未完成的结果不进入数值表，缺失项留空而非0。完全没有任何记录时明确报错；同一指标重复或单位/定义冲突时停止，本次不更新旧表。不同计算参数要在definition中注明，不能把不同定义混进同列。若要并列不同方法，请在metric或conditions中明确方法名；definition用于检查同一列是否可比。

## 自动扩展规则

03_results/<模拟编号>/<分类>/.../scalars.json会自动收集。分类名去掉数字前缀作为表名，例如01_basic→basic.csv；记录的object非空时再分对象表。

新增已有分类的指标，会自动增加列；新增分类，会增加表。04不维护指标名单，也不猜测任意CSV的含义。尚未输出过记录的计算，04无法凭空知道；已知记录在部分模拟中缺失时会列入status.csv。

## 计算端输出约定

详细计算在自身结果目录生成scalars.json，同时保留COMPLETE和metadata.json；记录中的source必须与元数据一致。用common.results.finish(..., scalars=records)即可自动生成。

单值计算可直接保存<指标>.scalars.json，使用start_scalar/finish_scalar标记完成，不再建一个只有一个值的文件夹。每条记录由scalar_record生成：metric、value、unit、object、conditions、definition、sem。基础指标object留空；多个指标可以放在同一计算的records列表中。

04每次从记录重新生成表，清理上次生成而本次已不存在的分类表。不要在output里手改汇总数值。并行任务只写各自记录，结束后再运行一次本入口。

可选--simulation B01_T260、--results-dir /path/to/results、--output-dir /path/to/summary。原有结果已补充记录，无需重算轨迹。
