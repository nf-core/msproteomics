#!/usr/bin/env bash
# Export FragPipe dry-run plans (--dry-run --export-plan) for every workflow in <workflows dir>.
#
#   capture.sh <workflows dir> <plans dir>
#
# FRAGPIPE_IMAGE: the licensed FragPipe 24.0 image (tools and licence in /fragpipe_bin/fragpipe-24.0/fragpipe-24.0)
# FRAGPIPE_LIB:   lib/ of a FragPipe build with --export-plan (branch export-plan-upstream-min); it replaces
#                 the image's lib in the container only
set -euo pipefail
WF_DIR=$(realpath "$1")
OUT=$(realpath -m "$2")
LIB=$(realpath "${FRAGPIPE_LIB:?set FRAGPIPE_LIB}")
mkdir -p "$OUT"
docker run --rm --cpus 8 --memory 64g \
    -v "$LIB:$LIB:ro" -v "$WF_DIR:$WF_DIR:ro" -v "$OUT:$OUT" --entrypoint bash "${FRAGPIPE_IMAGE:?set FRAGPIPE_IMAGE}" -c "
set -euo pipefail
I=/fragpipe_bin/fragpipe-24.0/fragpipe-24.0
rm -rf \$I/lib && cp -r $LIB \$I/lib
IN=/tmp/in; mkdir -p \$IN \$IN/exp1 \$IN/exp2
printf '>sp|P1|A_HUMAN a\nPEPTIDEK\n>rev_sp|P1|A_HUMAN a\nKEDITPEP\n' > \$IN/db.fasta
: > \$IN/a.mzML; : > \$IN/b.mzML
printf '%s\texp1\t1\tDDA\n%s\texp2\t1\tDDA\n' \$IN/a.mzML \$IN/b.mzML > \$IN/m.fp-manifest
# TMT runs: one directory per experiment with its <exp>_annotation.txt, which FragPipe finds next to the
# LC-MS files (TmtiPanel.java:868-879); the pool sample is the reference channel (ref_tag pool)
for e in exp1 exp2; do
  : > \$IN/\$e/\$e.mzML
  for c in 126 127N 127C 128N 128C 129N 129C 130N 130C 131N 131C 132N 132C 133N 133C 134N 134C 135N; do
    [ \$c = 126 ] && echo \"\$c pool_\$e\" || echo \"\$c \${e}_\$c\"
  done > \$IN/\$e/channels.txt
done
printf '%s\texp1\t1\tDDA\n%s\texp2\t1\tDDA\n' \$IN/exp1/exp1.mzML \$IN/exp2/exp2.mzML > \$IN/tmt.fp-manifest
for wf in $WF_DIR/*.workflow; do
  name=\$(basename \$wf .workflow)
  manifest=\$IN/m.fp-manifest
  if grep -q '^tmtintegrator.run-tmtintegrator=true' \$wf; then
    channels=\$(sed -n 's/^tmtintegrator.channel_num=TMT-//p' \$wf)
    manifest=\$IN/tmt.fp-manifest
    for e in exp1 exp2; do head -n \$channels \$IN/\$e/channels.txt > \$IN/\$e/\${e}_annotation.txt; done
  fi
  { grep -v '^database.db-path=' \$wf; echo database.db-path=\$IN/db.fasta; } > \$IN/\$name.workflow
  # a fresh user.home per run: FragPipe keeps the last run's panel values in fragpipe-ui.cache there
  home=/tmp/home/\$name; mkdir -p \$home /tmp/wd/\$name
  java -Duser.home=\$home -Djava.io.tmpdir=/tmp -cp \"\$I/lib/*\" org.nesvilab.fragpipe.FragPipeMain --headless \
      --workflow \$IN/\$name.workflow --manifest \$manifest --workdir /tmp/wd/\$name --threads 4 --ram 8 \
      --config-tools-folder \$I/lib/../tools --dry-run --export-plan $OUT/\$name.plan.xml > $OUT/\$name.log 2>&1 \
    && echo \"OK \$name\" || echo \"REFUSED \$name\"
done
chown -R $(id -u):$(id -g) $OUT
"
