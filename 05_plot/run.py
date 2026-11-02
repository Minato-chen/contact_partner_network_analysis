"""Top-level batch script for 05_plot.
Executes all plotting pipelines across Topic 1, Topic 2, and Cross-Topic Overview.
Strictly exports only PNG figures (300 DPI) as requested.
"""
from pathlib import Path
import sys
import argparse
import importlib

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

# Dynamic import because directory names start with digits
topic1_intra = importlib.import_module("01_topic1_blockiness.plot_topic1_intra")
topic1_inter = importlib.import_module("01_topic1_blockiness.plot_topic1_inter")
topic1_joint = importlib.import_module("01_topic1_blockiness.plot_topic1_joint")

topic2_intra = importlib.import_module("02_topic2_polarization.plot_topic2_intra")
topic2_inter = importlib.import_module("02_topic2_polarization.plot_topic2_inter")
topic2_joint = importlib.import_module("02_topic2_polarization.plot_topic2_joint")

cross_topic = importlib.import_module("00_cross_topic_overview.plot_cross_topic")


def run_all(object_names=("dense_chains", "central_chain")):
    print("=" * 60)
    print("Executing 05_plot batch figure generation (PNG only, 300 DPI)")
    print("=" * 60)
    
    for obj in object_names:
        print(f"\n>>> Target Object: {obj} <<<")
        print("\n--- Topic 1: Charge Blockiness (A18, B01, B10) ---")
        topic1_intra.main(obj)
        topic1_inter.main(obj)
        topic1_joint.main(obj)
        
        print("\n--- Topic 2: Charge Polarization & Neutral Tail (A20, C04, B09, C22) ---")
        topic2_intra.main(obj)
        topic2_inter.main(obj)
        topic2_joint.main(obj)
        
        print("\n--- Cross-Topic 7-Variant Overview ---")
        cross_topic.main(obj)

    print("\n" + "=" * 60)
    print("[SUCCESS] All 05_plot figures generated successfully.")
    print("=" * 60)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Batch generate 05_plot figures")
    parser.add_argument("--object", choices=["dense_chains", "central_chain", "both"], default="both")
    args = parser.parse_args()
    
    objs = ["dense_chains", "central_chain"] if args.object == "both" else [args.object]
    run_all(objs)
