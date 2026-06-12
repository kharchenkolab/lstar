#!/usr/bin/env bash
# Robust scRNAseq sweep: run EACH dataset in its own Rscript subprocess (timeout-guarded), so a
# hard crash (segfault) or hang in one dataset only fails that one -- not the whole sweep.
RLIB="$(cd "$(dirname "$0")/../.." && pwd)/.Rlib"
OUT=/tmp/sweep_scrnaseq.tsv; echo -e "dataset\tstatus\tfields\taxes\tnote" > "$OUT"
CALLS=$(Rscript -e ".libPaths(c('$RLIB',.libPaths())); cat(scRNAseq::listDatasets()\$Call, sep='\n')" 2>/dev/null)
i=0
for call in $CALLS; do
  i=$((i+1)); [ $i -gt 61 ] && break
  timeout 180 Rscript -e ".libPaths(c('$RLIB',.libPaths())); suppressWarnings(suppressMessages({library(scRNAseq);library(SingleCellExperiment);library(lstar)}))
    r<-tryCatch({sce<-$call; if(is.list(sce)&&!is(sce,'SummarizedExperiment'))sce<-sce[[1]]
      if(ncol(sce)>100000)stop('too-big'); ds<-read_sce(sce); write_sce(ds); p<-tempfile(fileext='.zarr'); lstar_write(ds,p)
      cat('PASS',length(ds\$fields),length(ds\$axes),sep='\t')},
      error=function(e){m<-conditionMessage(e); s<-if(grepl('there is no package|could not find|too-big|namespace|openpyxl',m))'LOADERR' else 'FAIL'; cat(s,'','',substr(m,1,80),sep='\t')})" \
    > /tmp/sw_one.txt 2>/dev/null
  rc=$?; nm=$(echo "$call" | sed "s/(.*//")
  if [ $rc -eq 124 ]; then echo -e "$nm\tTIMEOUT\t\t\t(>180s)" >> "$OUT"
  elif ! grep -qE "PASS|FAIL|LOADERR" /tmp/sw_one.txt; then echo -e "$nm\tCRASH\t\t\t(segfault/OOM)" >> "$OUT"
  else echo -e "$nm\t$(cat /tmp/sw_one.txt)" >> "$OUT"; fi
done
echo "done: $(grep -cP '\tPASS\t' $OUT) PASS / $(grep -cP '\tFAIL\t' $OUT) FAIL / $(grep -cP '\tLOADERR\t' $OUT) load-skip / $(grep -cP '\tCRASH\t' $OUT) crash / $(grep -cP '\tTIMEOUT\t' $OUT) timeout"
