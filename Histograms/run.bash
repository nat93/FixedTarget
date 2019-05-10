#!/bin/bash

make clean; make histograms;
./histograms 1 /media/andrii/F492773C92770302/SPS_SimulationData/hardccbar2.root



#echo "Bash version ${BASH_VERSION}..."

#make clean; make -j2;

#for i in {1..20..1}
#do
##	screen -X -S natochii_$i.thr kill;
#	echo "#!/bin/bash" > run_$i.sh;
#	echo "source /home/gred/root_34_36/bin/thisroot.sh;" >> run_$i.sh;
#	echo "time ./histograms 1 ../../home2/SPS_Sim_Data/hardccbar_$i.root;" >> run_$i.sh;
#	echo "time ./histograms 2 ../../home2/SPS_Sim_Data/hardccbar_$i.root;" >> run_$i.sh;
#	echo "exec sh" >> run_$i.sh;
#	chmod a+x run_$i.sh;
#	screen -S natochii_$i.thr -L -d -m ./run_$i.sh;
#done

#	echo "#!/bin/bash" > runAll_21.sh;
#	echo "source /home/gred/root_34_36/bin/thisroot.sh;" >> runAll_21.sh;
#	echo "time ./histograms 1 ../../home2/SPS_Sim_Data/hardccbar_.root;" >> runAll_21.sh;
#	echo "exec sh" >> runAll_21.sh;
#	chmod a+x runAll_21.sh;
#	screen -S natochii_21.thr -L -d -m ./runAll_21.sh;

#echo "All done";

#screen -ls;
