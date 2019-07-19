
for i in *bed; do echo $i;  bsub -q big -o /dev/null macs2 callpeak -t $i --nomodel --shift -100 --extsize 200 --keep-dup all -n "../peak_calls/${i}" -q 0.01; sleep 1; done


