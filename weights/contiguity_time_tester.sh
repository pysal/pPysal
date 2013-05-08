#!/bin/bash
#echo "Bash version ${BASH_VERSION}..."

python contiguity_apply_async.py
python contiguity_joinable_queue.py
python contiguity_managed_dict.py
python contiguity_map_async.py
python contiguity_serial.py

