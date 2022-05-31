#!/bin/bash

# categorize cog annotations
for FILE in *redo.table; do
	KEEP="${FILE%.c*}" 
	echo "$KEEP"
	awk -F"\t" 'NR==FNR {a[$1]=$2;next} {if ($2 in a){print $1, $2, a[$2]} else {print $0}}' OFS="\t" cognames2003-2014.tab "$FILE" > temp
	awk -F"\t" '{if ( length($3)>1 ) { $3 = substr($3, 0, 1) } else { $3 = $3 }; print}' OFS="\t" temp > temp2
	awk -F"\t" 'NR==FNR {a[$1]=$2;next} {if ($3 in a){print $0, a[$3]} else {print $0}}' OFS="\t" fun2003-2014.tab temp2 > "$KEEP".cog.redo.categorized
done
rm temp temp2

