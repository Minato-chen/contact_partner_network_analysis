"""Test Self-Neutralization and Sticker-Frustration Model in new_measure/.
Investigating how long-range attraction transitions from sticker to self-neutralization trap.
"""
from pathlib import Path
import sys
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, pearsonr

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from importlib import import_module
read_sequence = import_module("02_basic_calculate._sequence").read_sequence
import yaml

with open(ROOT / "simulations.yaml") as f:
    cfg = yaml.safe_load(f)

CSAT_POSTER = {
    "A18": 0.949, "B01": 0.145, "B10": 1.305,
    "A20": 0.125, "C04": 4.142, "B09": 1.212, "C22": 4.584,
}

variants = []
for s in cfg["simulations"]:
    var = s["variant"]
    src = ROOT / "01_center/output" / s["id"]
    seq, _ = read_sequence(src)
    q = np.array([1 if x in "KR" else (-1 if x in "DE" else 0) for x in seq], dtype=float)
    variants.append({
        "variant": var,
        "charges": q,
        "csat": CSAT_POSTER[var],
        "ln_csat": np.log(CSAT_POSTER[var]),
    })

# 物理模型：
# 1. 短程/中程异号电荷相互作用 (k <= k_cutoff): 形成分子表面可用的多价斑块 (Stickers), 促进相分离 (驱动力为正)
# 2. 超长程异号电荷相互作用 (k > k_cutoff): 导致头尾自锁 (Self-neutralizing hairpin), 破坏多价性 (惩罚项)
# 3. 同号电荷排斥: 始终产生静电膨胀，防止过度折叠

def compute_effective_charge_valence(q, k_cutoff=35, lambda_penalty=1.5):
    N = len(q)
    f_crosslink = 0.0
    f_trap = 0.0
    
    for m in range(1, N):
        for n in range(m):
            k = m - n
            prod = q[m] * q[n]
            if prod < 0:  # 异号吸引
                if k <= k_cutoff:
                    # 适度间距：暴露在分子表面的有效结合斑块 (Sticker)
                    f_crosslink += -prod * np.sqrt(k)
                else:
                    # 超长间距：分子内头尾自闭锁定 (Frustration trap)
                    f_trap += -prod * np.sqrt(k)
            elif prod > 0:  # 同号排斥
                # 同号电荷排斥促使链伸展，防止自闭
                f_crosslink -= prod * np.sqrt(k)
                
    # 净有效交联驱动力 = 表面多价驱动力 - 自闭锁定惩罚
    # 驱动力越大，相分离越强 (Csat 越低)，所以我们定义预测得分 Score:
    score = (f_crosslink - lambda_penalty * f_trap) / N
    return score

print("Scanning k_cutoff and lambda_penalty...")
best_k, best_lam, best_rho, best_r = None, None, 1.0, 1.0
best_grid = []

# 我们期望 Score 越大，Csat 越低 (即负相关, rho < 0 且接近 -1)
for k_cut in range(15, 60, 5):
    for lam in np.linspace(0.2, 3.0, 29):
        scores = [compute_effective_charge_valence(v["charges"], k_cutoff=k_cut, lambda_penalty=lam) for v in variants]
        csats = [v["csat"] for v in variants]
        ln_csats = [v["ln_csat"] for v in variants]
        
        rho, _ = spearmanr(scores, csats)
        r, _ = pearsonr(scores, ln_csats)
        best_grid.append({"k_cut": k_cut, "lambda": lam, "rho": rho, "r": r, "r2": r**2})

df_res = pd.DataFrame(best_grid)
# 寻找与 Csat 负相关最强的参数 (理论上 Score 越高，相分离越容易，Csat 越低)
best_match = df_res.sort_values("rho").iloc[0]
print("Best Match:")
print(best_match)

# 打印在该模型下的各变异体得分
opt_k = int(best_match["k_cut"])
opt_lam = float(best_match["lambda"])

for v in variants:
    v["eff_score"] = compute_effective_charge_valence(v["charges"], k_cutoff=opt_k, lambda_penalty=opt_lam)

df_out = pd.DataFrame(variants)[["variant", "csat", "eff_score"]].sort_values("csat")
print(f"\n--- Ranked by True Csat (lowest Csat = strongest phase separation) ---")
print(df_out)
