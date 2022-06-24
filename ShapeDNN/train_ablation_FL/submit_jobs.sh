#! /bin/bash

## first job - no dependencies
## generate train.test split from depelete/enrich shape files
#jid1a=$(sbatch run_yeast_deplete_step1.sh)
#jid1b=$(sbatch run_yeast_enrich_step1.sh)
#
## convert train test to .mat inputs
#echo ${jid1a##* }
#echo ${jid1b##* }
#jid2a=$(sbatch  --dependency=afterany:${jid1a##* }:${jid1b##* } run_yeast_test_step2.sh)
#jid2b=$(sbatch  --dependency=afterany:${jid1a##* }:${jid1b##* } run_yeast_train_step2.sh)

## start training with input file
#echo ${jid1a##* }
#echo ${jid2b##* } 
#jid3=$(sbatch  --dependency=afterany:${jid2a##* }:${jid2b##* } train_ShapeDNN_run1.sl)
#jid4=$(sbatch  --dependency=afterany:${jid2a##* }:${jid2b##* } train_ShapeDNN_run2.sl)
#jid5=$(sbatch  --dependency=afterany:${jid2a##* }:${jid2b##* } train_ShapeDNN_run3.sl)
#jid6=$(sbatch  --dependency=afterany:${jid2a##* }:${jid2b##* } train_ShapeDNN_run4.sl)
#jid7=$(sbatch  --dependency=afterany:${jid2a##* }:${jid2b##* } train_ShapeDNN_run5.sl)
#jid8=$(sbatch  --dependency=afterany:${jid2a##* }:${jid2b##* } train_ShapeDNN_run6.sl)

#echo ${jid2a##* }
#echo ${jid2b##* } 
jid3=$(sbatch  train_ShapeDNN_run1.sl)
jid4=$(sbatch  train_ShapeDNN_run2.sl)
jid5=$(sbatch  train_ShapeDNN_run3.sl)
jid6=$(sbatch  train_ShapeDNN_run4.sl)
jid7=$(sbatch  train_ShapeDNN_run5.sl)
jid8=$(sbatch  train_ShapeDNN_run6.sl)

echo ${jid3##* }
echo ${jid4##* }
echo ${jid5##* }
echo ${jid6##* }
echo ${jid7##* }
echo ${jid8##* }

# start test after trained models
jid3a=$(sbatch  --dependency=afterany:${jid3##* } test_ShapeDNN_run1.sl)
jid4a=$(sbatch  --dependency=afterany:${jid4##* } test_ShapeDNN_run2.sl)
jid5a=$(sbatch  --dependency=afterany:${jid5##* } test_ShapeDNN_run3.sl)
jid6a=$(sbatch  --dependency=afterany:${jid6##* } test_ShapeDNN_run4.sl)
jid7a=$(sbatch  --dependency=afterany:${jid7##* } test_ShapeDNN_run5.sl)
jid8a=$(sbatch  --dependency=afterany:${jid8##* } test_ShapeDNN_run6.sl)


echo ${jid3a##* }
echo ${jid4a##* }
echo ${jid5a##* }
echo ${jid6a##* }
echo ${jid7a##* }
echo ${jid8a##* }
