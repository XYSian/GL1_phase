#!/bin/sh
set -vex
export PATH=$PATH:/bin
cd /home/data6/fenglijun/lxy/falcon/5-phase/output
/bin/bash task.sh
touch /home/data6/fenglijun/lxy/falcon/5-phase/output/run.sh.done
