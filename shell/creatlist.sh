#!/bin/bash

export EOS_MGM_URL=root://eos01.ihep.ac.cn/
#Year1=2020
#Date1=$1
for i in `seq 249 249`
do
Year1=`date +%Y -d "-${i} day"`
Date1=`date +%m%d -d "-${i} day"`

echo $Year1 $Date1

infiledir=/eos/lhaaso/raw/wcda
filelistdir=/workfs/ybj/youzhiyong/WCDA/WCDARec/iptlist/${Year1}

mkdir -p $filelistdir

eos newfind -f $infiledir/$Year1/$Date1 | grep "ES" | grep "reduced.root" | grep -v ".sys." >$filelistdir/${Year1}${Date1}.txt

sed -i 's/^path=//g' $filelistdir/${Year1}${Date1}.txt

done
