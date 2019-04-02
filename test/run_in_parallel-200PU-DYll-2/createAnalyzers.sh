jobName="2019_Mar23-dyll-0p5-ext4"
#
j=0
for i in {0..75}
do
    cat test-Analyzer.py > SUB-Analyzer-${i}.py
    #echo "process.source.skipEvents = cms.untracked.uint32(${j})" >> SUB-Analyzer-${i}.py
    j=$(( $j + 100))
    cat submit-$i.py >> SUB-Analyzer-${i}.py
    cat submit.py >> SUB-Analyzer-${i}.py

    mkdir -p /nfs_scratch/ojalvo/${jobName}-DYLL-200PU-SUBPhase-$i/dags/daginputs
    
    farmoutAnalysisJobs --vsize-limit=7000 --site-requirements='OpSysAndVer == "SL6"' --assume-input-files-exist  --input-file-list=inputFileList-mini.txt --submit-dir=/nfs_scratch/ojalvo/${jobName}-DYLL-200PU-SUBPhase-$i/submit --output-dag-file=/nfs_scratch/ojalvo/${jobName}-DYLL-200PU-SUBPhase-$i/dags/dag  ${jobName}-DYLL-200PU  $CMSSW_BASE  $CMSSW_BASE/src/L1Trigger/L1MTD/test/run_in_parallel-200PU-DYll-2/SUB-Analyzer-$i.py     &


    if [ "$i" -eq "15" ]; then
	wait;
    fi

    if [ "$i" -eq "30" ]; then
	wait;
    fi

    if [ "$i" -eq "45" ]; then
	wait;
    fi

    if [ "$i" -eq "60" ]; then
	wait;
    fi


done



