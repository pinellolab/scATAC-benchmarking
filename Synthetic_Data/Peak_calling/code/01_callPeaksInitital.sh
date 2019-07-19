
readFile="/data/aryee/caleb/hemeATAC/Bcell_reads.bed"
bsub -q big macs2 callpeak -t $readFile --nomodel --shift -100 --extsize 200 --keep-dup all -n ../peak_calls/Bcell_allData -q 0.01

readFile="/data/aryee/caleb/hemeATAC/CD4_reads.bed"
bsub -q big macs2 callpeak -t $readFile --nomodel --shift -100 --extsize 200 --keep-dup all -n ../peak_calls/CD4_allData -q 0.01


readFile="/data/aryee/caleb/hemeATAC/Mono_reads.bed"
bsub -q big macs2 callpeak -t $readFile --nomodel --shift -100 --extsize 200 --keep-dup all -n ../peak_calls/Mono_allData -q 0.01


