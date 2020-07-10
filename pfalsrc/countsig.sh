#! /bin/sh

sigfiles=`ls -d *sig*`

for entry in $sigfiles
do
	echo `wc -l $entry` > sig_cis.txt
 
 
done