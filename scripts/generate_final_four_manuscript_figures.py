#!/usr/bin/env python3
"""Regenerate the four final manuscript data figures from archived inputs.

Run: python scripts/generate_final_four_manuscript_figures.py
Writes to figures/manuscript/reproduced/; never overwrites submitted images.
"""
from pathlib import Path
from collections import defaultdict
import csv
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "figures/manuscript/reproduced"
LENGTHS = (30, 40, 50, 70)
THRESHOLDS = (0, 10, 30)

def rows(path, delimiter=","):
    with (ROOT / path).open(newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f, delimiter=delimiter))

def save(fig, name):
    OUT.mkdir(parents=True, exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(OUT / (name + "." + ext), dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("Saved", name)

def main():
    summary = rows("results/circular_pmd_replicate_summary.csv")
    all_rows = rows("results/circular_pmd_replicates_all.csv")
    assert len(summary) == 144 and len(all_rows) == 2976, "Unexpected result CSV row count"
    summary_idx = {(r["length"],r["damage"],r["aligner"],r["reference"],r["group"]):r for r in summary}
    def value(L, aligner, ref, group, col):
        r = summary_idx[(str(L),"high",aligner,ref,group)]
        assert int(r["replicates"]) == 10
        return 100*float(r[col])

    groups = ("original_boundary","shifted_boundary","interior")
    x=np.arange(3)
    fig,axes=plt.subplots(1,2,figsize=(12,4.6),sharey=True)
    for ax,aligner in zip(axes,("bwa","bowtie2")):
        for offset,ref,label in ((-0.17,"original","Original reference"),(0.17,"shifted","Shifted reference")):
            ax.bar(x+offset,[value(30,aligner,ref,g,"mean_recall") for g in groups],0.34,label=label)
        ax.set_xticks(x,["Original\nboundary","Shifted\nboundary","Interior"])
        ax.set_ylim(0,110)
        ax.set_title("BWA-aln" if aligner=="bwa" else "Bowtie2")
    axes[0].set_ylabel("Strict recall (%)")
    h,l=axes[0].get_legend_handles_labels()
    fig.legend(h,l,loc="lower center",ncol=2)
    fig.tight_layout(rect=(0,0.09,1,1))
    save(fig,"figure1_boundary_groups_30bp")

    fig,ax=plt.subplots(figsize=(12,4.6))
    for aligner,marker,label in (("bwa","o","BWA-aln"),("bowtie2","s","Bowtie2")):
        means=[value(L,aligner,"shifted","original_boundary","mean_recall") for L in LENGTHS]
        sds=[value(L,aligner,"shifted","original_boundary","sd_recall") for L in LENGTHS]
        ax.errorbar(LENGTHS,means,yerr=sds,marker=marker,capsize=4,label=label)
    ax.set(xlabel="Fragment length (bp)",ylabel="Strict recall (%)",ylim=(89,102))
    ax.set_xticks(LENGTHS)
    ax.legend()
    fig.tight_layout()
    save(fig,"figure2_length_recall_high_pmd")

    errors=defaultdict(int)
    seen=set()
    for r in all_rows:
        if r["damage"]!="high" or r["aligner"]!="bowtie2": continue
        ref=r["reference"]
        if r["group"]!=("original_boundary" if ref=="original" else "shifted_boundary"):continue
        L=int(r["length"])
        errors[(ref,L)]+=int(r["incorrect"])
        seen.add((ref,L,r["replicate"],r["strand"]))
    assert len(seen)==160, "Missing high-PMD replicate records"
    assert [errors[("original",L)] for L in LENGTHS]==[0,30,112,266]
    assert [errors[("shifted",L)] for L in LENGTHS]==[0,22,73,218]
    fig,ax=plt.subplots(figsize=(12,4.2))
    x=np.arange(4)
    ax.bar(x-0.175,[errors[("original",L)] for L in LENGTHS],0.35,label="Original reference")
    ax.bar(x+0.175,[errors[("shifted",L)] for L in LENGTHS],0.35,label="Shifted reference")
    ax.set_xticks(x,[str(L) for L in LENGTHS])
    ax.set(xlabel="Fragment length (bp)",ylabel="Incorrect active-boundary alignments")
    ax.legend()
    fig.tight_layout()
    save(fig,"figure3_bowtie2_boundary_errors")

    counts=defaultdict(int)
    records=0
    for L in LENGTHS:
        truth_rows=rows(f"ground_truth/circular_length_baseline/reads_{L}bp_none.tsv","\t")
        truth={r["read_id"]:r for r in truth_rows}
        assert len(truth)==len(truth_rows)
        for rep in range(1,11):
            path=(ROOT/"alignments"/("circular_pmd_main" if rep==1 else "circular_pmd_replicates")/
                  f"bowtie2_original_{L}bp_high{'' if rep==1 else '_rep'+str(rep)}.sam")
            seen_ids=set()
            with path.open(encoding="utf-8") as f:
                for line in f:
                    if line.startswith("@"):continue
                    a=line.rstrip("\n").split("\t")
                    if len(a)<11:raise ValueError("Malformed SAM record: "+str(path))
                    rid,flag_s,_,pos_s,mapq_s,cigar=a[:6]
                    flag=int(flag_s)
                    if flag & 0x900:continue
                    if rid not in truth or rid in seen_ids:raise ValueError("Unexpected or duplicate primary record")
                    seen_ids.add(rid)
                    tr=truth[rid]
                    unmapped=bool(flag & 4)
                    correct=(not unmapped and int(pos_s)==int(tr["start_1based"])
                             and ("-" if flag&16 else "+")==tr["strand"] and cigar==f"{L}M")
                    for th in THRESHOLDS:
                        status=("filtered" if unmapped or int(mapq_s)<th else
                                "correct" if correct else "incorrect")
                        counts[(L,tr["group"],th,status)]+=1
            if len(seen_ids)!=len(truth):raise ValueError(f"Missing primary reads in {path}")
            records+=len(seen_ids)
    boundary=[sum(counts[(L,"original_boundary",t,"incorrect")] for L in LENGTHS) for t in THRESHOLDS]
    assert boundary==[408,1,0],f"Unexpected MAPQ retained error counts: {boundary}"
    fig,axes=plt.subplots(1,2,figsize=(12,4.2))
    for L in LENGTHS:
        denominator=10*2*(L-1)
        recalls=[100*counts[(L,"interior",t,"correct")]/denominator for t in THRESHOLDS]
        axes[0].plot(THRESHOLDS,recalls,marker="o",label=f"{L} bp")
    axes[0].set(xlabel="Minimum MAPQ",ylabel="Interior strict recall (%)",ylim=(60,102))
    axes[0].set_xticks(THRESHOLDS)
    axes[0].legend(ncol=2)
    axes[1].bar([str(t) for t in THRESHOLDS],boundary)
    axes[1].set(xlabel="Minimum MAPQ",ylabel="Incorrect boundary mappings retained")
    fig.tight_layout()
    save(fig,"figure4_bowtie2_mapq_filtering")
    print(f"Validated {records} primary SAM records; retained boundary errors: {boundary}")
    print("Four manuscript data figures generated; original figures untouched.")

if __name__=="__main__":
    main()
