#!/usr/bin/env bash
# provenance for one SRA accession
# usage: write_manifest.sh <sra_id> <max_reads>
set -euo pipefail
sra_id="$1"; max_reads="$2"
files=""; sep=""
for f in "${sra_id}"*.fastq.gz; do
  files="${files}${sep}\"${f}\": \"$(sha256sum "$f" | cut -d' ' -f1)\""
  sep=", "
done
printf '{"sra_id": "%s", "fetched_at": "%s", "subsample_max_reads": %s, "files": {%s}}\n' \
  "$sra_id" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$max_reads" "$files" > "${sra_id}.download.json"
