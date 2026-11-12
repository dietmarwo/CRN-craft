#!/bin/bash
# Alternating 0.6.0 / 0.7.0 timings of the mtDNA workload, each pair under the CPU lock.
M=/home/claude/m07/wt-int/models/rules/mtdna-three-class-cme.rules.json
OUT=/home/claude/m07/mtdna-speed/timings.jsonl
OLD=/home/claude/m07/venv060/bin/python
NEW=/home/claude/m07/venv-m6/bin/python
for threads in 1 2; do
  for spec in "mtdna_direct 1000" "static_direct 400" "tree 100"; do
    set -- $spec
    for rep in 1 2 3 4 5; do
      key="$1-$threads-$rep"
      grep -q "\"$key\"" /home/claude/m07/mtdna-speed/done.txt 2>/dev/null && continue
      if [ $((rep % 2)) -eq 1 ]; then first=$OLD; second=$NEW; else first=$NEW; second=$OLD; fi
      flock /home/claude/craft5/cpu.lock bash -c "cd /tmp && $first -I /home/claude/m07/mtdna-speed/bench.py $M $1 $threads $2 $OUT >/dev/null && $second -I /home/claude/m07/mtdna-speed/bench.py $M $1 $threads $2 $OUT >/dev/null" && echo "\"$key\"" >> /home/claude/m07/mtdna-speed/done.txt
    done
  done
done
echo finished >> /home/claude/m07/mtdna-speed/done.txt
