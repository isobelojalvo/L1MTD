jobName="2018_Nov29-mu"
#
j=0
for i in {0..89}
do
    cat test-Analyzer-grow-l1t.py > SUB-Analyzer-${i}.py
    #echo "process.source.skipEvents = cms.untracked.uint32(${j})" >> SUB-Analyzer-${i}.py
    j=$(( $j + 100))
    cat submit-$i.py >> SUB-Analyzer-${i}.py
    cat submit.py >> SUB-Analyzer-${i}.py

    mkdir -p /nfs_scratch/ojalvo/${jobName}-ZMM-200PU-SUBPhase-$i/dags/daginputs
    
    farmoutAnalysisJobs --site-requirements='OpSysAndVer == "SL6"' --assume-input-files-exist  --input-file-list=inputFileList.txt --submit-dir=/nfs_scratch/ojalvo/${jobName}-ZMM-200PU-SUBPhase-$i/submit --output-dag-file=/nfs_scratch/ojalvo/${jobName}-ZMM-200PU-SUBPhase-$i/dags/dag  ${jobName}-ZMM-200PU  $CMSSW_BASE  $CMSSW_BASE/src/L1Trigger/L1MTD/test/run_in_parallel-200PU-mu/SUB-Analyzer-$i.py     

done



