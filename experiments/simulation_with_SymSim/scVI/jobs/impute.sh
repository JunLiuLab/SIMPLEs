#!/bin/bash

impute_script="scvi_imputation.py"

for i in {1..20};
do
    for p in umi nonumi;
    do
        python ../${impute_script} --seed ${i} --platform ${p}
    done
done
