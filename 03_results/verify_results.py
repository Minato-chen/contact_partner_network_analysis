"""Verify completeness of 03_results across all simulated variants."""
from pathlib import Path
import sys
import yaml

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
RESULTS_ROOT = HERE

with open(ROOT / "simulations.yaml") as f:
    cfg = yaml.safe_load(f)

simulations = cfg["simulations"]

CHECKS = [
    # (Category, Task, Object/File, Subpath, Description)
    ("01_basic", "phase_boundary", "", "01_basic/phase_boundary/COMPLETE", "Phase boundary"),
    ("01_basic", "csat", "", "01_basic/csat/COMPLETE", "Dilute Csat"),
    ("01_basic", "rg", "", "01_basic/rg/COMPLETE", "All-chain Rg"),
    ("01_basic", "kappa", "", "01_basic/kappa.scalars.json", "Kappa parameter"),
    ("01_basic", "scd", "", "01_basic/scd.scalars.json", "SCD parameter"),
    ("02_contacts", "intra", "central", "02_contacts/01_intra/central_chain/COMPLETE", "Intra (central)"),
    ("02_contacts", "intra", "dense", "02_contacts/01_intra/dense_chains/COMPLETE", "Intra (dense)"),
    ("02_contacts", "inter", "central", "02_contacts/02_inter/central_chain/COMPLETE", "Inter (central)"),
    ("02_contacts", "inter", "dense", "02_contacts/02_inter/dense_chains/COMPLETE", "Inter (dense)"),
    ("02_contacts", "joint", "central", "02_contacts/03_joint/central_chain/COMPLETE", "Joint (central)"),
    ("02_contacts", "joint", "dense", "02_contacts/03_joint/dense_chains/COMPLETE", "Joint (dense)"),
]


def main():
    print(f"\n{'='*78}")
    print(f" 03_results Verification Report ({len(simulations)} simulations)")
    print(f"{'='*78}")
    
    header = f"{'Simulation':<12} | " + " | ".join(f"{c[4]:<14}" for c in CHECKS)
    print(f"{'Simulation':<12} | Basic (5 items) | Contacts (6 items) | Status")
    print(f"{'-'*12}-+-{'-'*16}-+-{'-'*18}-+-{'-'*10}")
    
    all_ok = True
    for sim in simulations:
        sid = sim["id"]
        sim_dir = RESULTS_ROOT / sid
        
        basic_ok = 0
        contacts_ok = 0
        
        for c in CHECKS[:5]:
            if (sim_dir / c[3]).is_file():
                basic_ok += 1
        for c in CHECKS[5:]:
            if (sim_dir / c[3]).is_file():
                contacts_ok += 1
                
        status = "COMPLETE" if (basic_ok == 5 and contacts_ok == 6) else "INCOMPLETE"
        if status != "COMPLETE":
            all_ok = False
        print(f"{sid:<12} | {basic_ok:>5}/5 OK      | {contacts_ok:>6}/6 OK       | {status}")
        
    print(f"{'='*78}")
    if all_ok:
        print(" [SUCCESS] All 7 simulations have 100% complete result artifacts.")
    else:
        print(" [WARNING] Some simulations have missing results.")
    print(f"{'='*78}\n")
    return 0 if all_ok else 1


if __name__ == "__main__":
    sys.exit(main())
