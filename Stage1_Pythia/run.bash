make clean; make -j2;

# threebody decay
#./hardccbar1 270 1000 /media/andrii/F492773C92770302/SPS_SimulationData/hardccbar1.root | tee hardccbar1.log

# only stable
./hardccbar2 270 1000000 /media/andrii/F492773C92770302/SPS_SimulationData/hardccbar2.root | tee hardccbar2.log


#useroot; time ./hardccbar2 270 20000000 /home/gred/home2/SPS_Sim_Data/hardccbar_1.root | tee hardccbar_1.log
