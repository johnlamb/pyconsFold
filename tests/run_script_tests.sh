#!/bin/bash
# Test so all example scripts work
test_folder=`dirname "$0"`"/data/"
base_test_folder=`dirname "$0"`"/../examples/"

# Test normal modelling
${base_test_folder}run_model.py ${test_folder}short.fasta ${test_folder}short.rr script_test_out/model -ss ${test_folder}short.ss -stage2
echo "Passed test 1 of 3"

# Test distance modelling
${base_test_folder}run_model_dist.py ${test_folder}2HIUA.fasta ${test_folder}2HIUA.rr script_test_out/model_dist -stage2
echo "Passed test 2 of 3"

# Test docking modelling
${base_test_folder}run_model_dock.py ${test_folder}2HIUA.fasta ${test_folder}2HIUB.fasta ${test_folder}2HIU_dock.rr script_test_out/model_dock -stage2
echo "Passed test 3 of 3"
