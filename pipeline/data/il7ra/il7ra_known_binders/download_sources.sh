#!/bin/bash
# Re-download the public sources used to build known_il7ra_binders.fasta
set -e
# 1) Proteinbase full snapshot (filter target == il7r)
curl -sL -o proteinbase.parquet "https://huggingface.co/api/datasets/adaptyvbio/proteinbase/parquet/default/proteins/0.parquet"
# 2) Anthropic binder-design release (filter target == IL-7Ra)
B=https://huggingface.co/datasets/Anthropic/claude-protein-binder-design/resolve/main/data
curl -sL -o design_summary.csv $B/tables/design_summary.csv
curl -sL -o designs.fasta $B/designs/designs.fasta
# 3) 7OPB minibinder
curl -sL https://www.rcsb.org/fasta/entry/7OPB/display
sha256sum proteinbase.parquet design_summary.csv designs.fasta
