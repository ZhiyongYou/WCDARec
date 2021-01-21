#!/bin/bash

export EOS_MGM_URL=root://eos01.ihep.ac.cn/

joboutdir=/lhaasofs/user/youzhiyong/WCDA/WCDARec/job
shelldir=/workfs/ybj/youzhiyong/WCDA/WCDARec/shell

for i in `seq 249 249`
do

	Year1=`date +%Y -d "-${i} day"`
	Date1=`date +%m%d -d "-${i} day"`

	echo $Year1 $Date1
	mkdir -p $joboutdir/$Year1/$Date1

	filelistdir=/workfs/ybj/youzhiyong/WCDA/WCDARec/iptlist/${Year1}

	for ifile in `cat $filelistdir/${Year1}${Date1}.txt`
	do
		b_file=`basename $ifile`
		echo $Year1 $Date1 $b_file
		hep_sub -o $joboutdir/$Year1/$Date1/ -e $joboutdir/$Year1/$Date1/ -g lhaaso $shelldir/Caliandrec.sh -argu $Year1 $Date1 $b_file
	done


done
