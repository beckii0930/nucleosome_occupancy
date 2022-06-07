#! /bin/bash

dir_name=train_ablation_FL_run2/ 
mkdir $dir_name
cp test_ShapeDNN_run*.sl $dir_name
cp train_ShapeDNN_run*.sl $dir_name
cp trainShapeDNN.py $dir_name
cp testShapeDNN.py $dir_name
cp PreprocessShape.py $dir_name        
cp cluster2shape_makeTest.py $dir_name
cp cluster2shape_makeTrain.py $dir_name
cp submit_jobs.sh $dir_name
cp change_dir.sh $dir_name

new_sub="sed -i 's+ShapeDNN/+ShapeDNN/$dir_name+g' *_ShapeDNN_run*.sl"
cd $dir_name && echo $new_sub > change_dir.sh && ./change_dir.sh

#chmod 777 $d | ir_name'change_dir.sh'
#change_dir.sh'
