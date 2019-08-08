Has been tested with Ubuntu 17.04 64Bit + Matlab 2015 

1. compile mhs_distmapC.cc via 
mex mhs_distmapC.cc -lgomp CXXFLAGS=" -O3   -Wfatal-errors  -std=c++11 -fopenmp-simd -fopenmp  -fPIC -march=native" 
2. type help mhs_simualte_axions
to get furhter infos about parameters to play with
3. run 
data=mhs_simualte_axions('shape',[512,512,512],'mode','chaos','density',1,'max_b',1,'axon_thickness',[2,3],'thick_min',1.5,'flatty_axions',true);
to get a dense sample of a "chaotic" group of axons, or
data=mhs_simualte_axions('shape',[256,256,512],'mode','bundle','density',5,'dir_noise_bias',0.005,'bundle_rad',0.75);
to get a dense bundle


 license:
  -------
    GPL 3.0
    see http://www.gnu.org/licenses/