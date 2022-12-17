cd $3structural_effects/pipeline_strucomparison
mkdir $3structural_effects/$2
cp $3structural_effects/data/Spike.pdb $3structural_effects/$2
python3 pipeline.py -fo $3structural_effects/$2/ -rp Spike.pdb -rf $1 -lm None -rt 1 -pr $3

