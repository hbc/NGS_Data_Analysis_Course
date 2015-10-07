#!/bin/bash

for fq in ~/unix_oct2015/raw_fastq/*
do
  	sh ~/rnaseq_analysis_on_allfiles-echo.sh $fq > ~/temp
	bsub < ~/temp
	rm ~/temp
  	sleep 1
done



