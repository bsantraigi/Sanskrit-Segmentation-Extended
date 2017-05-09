#!/bin/bash
# A simple script

echo ''
echo '##############   DeepR multilayer perceptron (800, 800)   ######################'
echo ''

tail -20 logs/bigram_1L_mir_h800.log

echo ''
echo '###############  Learning Rate (3e-6, 5e-6) + dropout(k0.8)  #####################'
echo ''

tail -20 logs/bigram_1L_mir_drop_80_slow.log

echo ''
echo '####################################'
echo ''
