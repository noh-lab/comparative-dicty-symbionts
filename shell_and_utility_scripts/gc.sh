#!/bin/bash

for FILE in GCF*.fna; do
	echo "$FILE"
	var1=$(grep -v ">" `echo "$FILE"` | tr -d -c GCgc | wc -c)
	var2=$(grep -v ">" `echo "$FILE"` | tr -d -c ATGCatgc | wc -c)
	#gc=`echo "$var1" / "$var2" | bc`
	echo "$var1, $var2"	
	awk -v v1=$var1 -v v2=$var2 'BEGIN {print v1/v2}'
done

