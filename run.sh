#!/bin/bash

for wind in {1..20..1}
do
    rm -f wind.dat
    echo "$wind" >> wind.dat
    ./ps_1d.out
    wait
    dir=results/exp_v"$wind"_a1_goddard_savi_4m
    echo $dir
    cp -r results/exp $dir
    wait
done