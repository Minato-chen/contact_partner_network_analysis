"""Test Multivalent Sticker Dispersion and Self-Locking Theory in new_measure/.
"""
from pathlib import Path
import sys
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, pearsonr
import matplotlib.pyplot as plt

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

# 物理创新假说：
# 经典高分子相分离理论（Rubinstein-Semenov, Choi-Pappu）：
# 相分离能力取决于有效多聚价态 (Effective Sticker Valence, v_eff)。
# 如果一个链的正负电荷高度集中成各 1 个超大块 (如 B10, C04)，它的价态退化为二价 (v = 2), 无法逾渗!
# 只有当一种电荷集中(作为高亲和核心), 另一种电荷沿链多点分散(作为多价分支点), 才能形成 v >> 2 的三维网络!

def analyze_charge_architecture(q):
    """提取电荷拓扑特征"""
    pos_idx = np.where(q > 0)[0]
    neg_idx = np.where(q < 0)[0]
    
    # 1. 负电荷的空间分散度 (Negative charge spatial spread / dispersion)
    # std 越大说明负电荷分散在整条链上，提供了多个独立的交联锚点
    neg_spread = np.std(neg_idx) if len(neg_idx) > 1 else 0.0
    
    # 2. 负电荷之间的平均间距 (Spacing between negative stickers)
    neg_spacing = np.mean(np.diff(neg_idx)) if len(neg_idx) > 1 else 1.0
    
    # 3. 头尾极端长程对消指数 (Terminal head-to-tail lock intensity)
    # 计算正负电荷距离大于 70 的相互作用强度
    terminal_lock = 0.0
    for p in pos_idx:
        for n in neg_idx:
            if abs(p - n) > 60:
                terminal_lock += 1.0
    terminal_lock /= (len(pos_idx) * len(neg_idx))
    
    # 4. 链两端非对称跨度 (Effective Span / Coverage of the chain)
    # 粘性基团覆盖链长的范围
    all_sticky = np.concatenate([pos_idx, neg_idx])
    chain_coverage = (all_sticky.max() - all_sticky.min()) / len(q)
    
    return {
        "neg_spread": neg_spread,
        "neg_spacing": neg_spacing,
        "terminal_lock": terminal_lock,
        "chain_coverage": chain_coverage,
    }

for v in variants:
    v.update(analyze_charge_architecture(v["charges"]))

df_arch = pd.DataFrame(variants).sort_values("csat")
print(df_arch[["variant", "csat", "neg_spread", "terminal_lock", "chain_coverage"]])

# 建立综合物理描述符：
# Effective Percolation Valence Index (EPVI):
# 越分散、覆盖度越大、末端自锁越小，则 EPVI 越大，Csat 越低！
# EPVI = (neg_spread * chain_coverage) / (1 + beta * terminal_lock)

for beta in [0.0, 1.0, 2.0, 3.0, 5.0]:
    scores = [(r.neg_spread * r.chain_coverage) / (1.0 + beta * r.terminal_lock) for _, r in df_arch.iterrows()]
    rho, _ = spearmanr(scores, df_arch["csat"])
    r, _ = pearsonr(scores, df_arch["ln_csat"])
    print(f"Beta = {beta:.1f} -> Spearman rho: {rho:.3f}, Pearson r with ln(Csat): {r:.3f} (R2={r**2:.3f})")
