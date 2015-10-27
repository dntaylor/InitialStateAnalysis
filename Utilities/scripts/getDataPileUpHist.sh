# setup environment for farmout summary
export PATH=/afs/hep.wisc.edu/cms/cmsprod/farmoutCmsJobs/:$PATH
# get processed json
jobReportSummary.py /nfs_scratch/dntaylor/2015-10-24_WZNtuples_13TeV_25ns_v1/data_*_Run2015*_*_25ns/submit/*/*.xml  --json-out data.json
# calculate pileup
pileupCalc.py -i data.json --inputLumiJSON /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions15/13TeV/PileUp/pileup_latest.txt --calcMode true  --minBiasXsec 80000 --maxPileupBin 80 --numPileupBins 80 PileUpData.root
pileupCalc.py -i data.json --inputLumiJSON /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions15/13TeV/PileUp/pileup_latest.txt --calcMode true  --minBiasXsec 88000 --maxPileupBin 80 --numPileupBins 80 PileUpData_up.root
pileupCalc.py -i data.json --inputLumiJSON /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions15/13TeV/PileUp/pileup_latest.txt --calcMode true  --minBiasXsec 72000 --maxPileupBin 80 --numPileupBins 80 PileUpData_down.root
