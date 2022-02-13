#!/bin/bash
asteroid=$1;
min_amp=$2;
./aspect ${asteroid} ${min_amp};
python stats.py ${asteroid};
# rm "${asteroid}-S_type_models.txt";
# rm "${asteroid}-C_type_models.txt";