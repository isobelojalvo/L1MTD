jobName="2019_Mar27-qcd-LLJets"
#
j=0

for j in {0..40}
do
    for i in {0..24}
    do
	k=$(( $j * 25 ))
	k=$(( $k + $i ))
	cat test-Analyzer.py > SUB-Analyzer-${k}.py
    #echo "process.source.skipEvents = cms.untracked.uint32(${j})" >> SUB-Analyzer-${i}.py
	cat submit.py >> SUB-Analyzer-${k}.py
	cat submit-$k.py >> SUB-Analyzer-${k}.py
	
	mkdir -p /nfs_scratch/ojalvo/${jobName}-QCD-200PU-SUBPhase-$k/dags/daginputs
	
	farmoutAnalysisJobs --vsize-limit=7000 --site-requirements='OpSysAndVer == "SL6"' --assume-input-files-exist  --input-file-list=inputFileList-mini.txt --submit-dir=/nfs_scratch/ojalvo/${jobName}-QCD-200PU-SUBPhase-$k/submit --output-dag-file=/nfs_scratch/ojalvo/${jobName}-QCD-200PU-SUBPhase-$k/dags/dag  ${jobName}-QCD-200PU  $CMSSW_BASE  $CMSSW_BASE/src/L1Trigger/L1MTD/test/run_in_parallel-200PU-QCD/SUB-Analyzer-$k.py     &
    done
    wait;
done



