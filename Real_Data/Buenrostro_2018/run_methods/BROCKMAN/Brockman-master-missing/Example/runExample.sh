#!/bin/bash
echo downloading example data
./downloadExample.sh 

#running brockman_pipeline for all samples, one at a time using config.sh
echo Running brockman_pipeline for sample 1
../brockman_pipeline config.sh 1
echo Running brockman_pipeline for sample 2
../brockman_pipeline config.sh 2
echo Running brockman_pipeline for sample 3
../brockman_pipeline config.sh 3
