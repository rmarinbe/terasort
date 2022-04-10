#!/bin/sh
qsub -I -q coc-ice -l nodes=1:ppn=16:core24,walltime=00:20:00 -N lab3_interactive
exit 0
