#!/bin/bash
cd run
for FILE in *.png; do
	mv "${FILE}" /data/Project/Codes/NFLBB/plotsdata/Graphs/${FILE}
    
done
