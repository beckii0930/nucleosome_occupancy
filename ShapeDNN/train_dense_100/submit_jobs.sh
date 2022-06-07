#! /bin/bash

## first job - no dependencies
jid3=$(sbatch  train_ShapeDNN_run1.sl)
jid4=$(sbatch  train_ShapeDNN_run2.sl)
jid5=$(sbatch  train_ShapeDNN_run3.sl)
jid6=$(sbatch  train_ShapeDNN_run4.sl)
jid7=$(sbatch  train_ShapeDNN_run5.sl)
jid8=$(sbatch  train_ShapeDNN_run6.sl)

echo $jid3
echo $jid4
echo $jid5
echo $jid6
echo $jid7
echo $jid8

# start test after trained models
jid3a=$(sbatch  --dependency=afterany:${jid3##* } test_ShapeDNN_run1.sl)
jid4a=$(sbatch  --dependency=afterany:${jid4##* } test_ShapeDNN_run2.sl)
jid5a=$(sbatch  --dependency=afterany:${jid5##* } test_ShapeDNN_run3.sl)
jid6a=$(sbatch  --dependency=afterany:${jid6##* } test_ShapeDNN_run4.sl)
jid7a=$(sbatch  --dependency=afterany:${jid7##* } test_ShapeDNN_run5.sl)
jid8a=$(sbatch  --dependency=afterany:${jid8##* } test_ShapeDNN_run6.sl)

