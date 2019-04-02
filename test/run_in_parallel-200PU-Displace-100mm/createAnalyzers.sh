jobName="2019_Mar13-100mm"
#
j=0
for i in {0..100}
do
    cat test-Analyzer.py > SUB-Analyzer-${i}.py
    #echo "process.source.skipEvents = cms.untracked.uint32(${j})" >> SUB-Analyzer-${i}.py
    j=$(( $j + 100))
    cat submit-$i.py >> SUB-Analyzer-${i}.py
    cat submit.py >> SUB-Analyzer-${i}.py

    mkdir -p /nfs_scratch/ojalvo/${jobName}-stopChi-100mm-200PU-SUBPhase-$i/dags/daginputs
    
    farmoutAnalysisJobs --vsize-limit=7000 --site-requirements='OpSysAndVer == "SL6"' --assume-input-files-exist  --input-file-list=inputFileList-mini.txt --submit-dir=/nfs_scratch/ojalvo/${jobName}-stopChi-100mm-200PU-SUBPhase-$i/submit --output-dag-file=/nfs_scratch/ojalvo/${jobName}-stopChi-100mm-200PU-SUBPhase-$i/dags/dag  ${jobName}-stopChi-100mm-200PU  $CMSSW_BASE  $CMSSW_BASE/src/L1Trigger/L1MTD/test/run_in_parallel-200PU-Displace-100mm/SUB-Analyzer-$i.py     

done



