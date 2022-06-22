#!/bin/tcsh

#module load /apps/modulefiles/singularity/3.4.0

#gcc and root
#gcc 4.9.1
#source /cvmfs/sft.cern.ch/lcg/external/gcc/4.9.1/x86_64-slc6/setup.csh
#root 5.34.36
#source /cvmfs/sft.cern.ch/lcg/app/releases/ROOT/5.34.36/x86_64-slc6-gcc49-opt/root/bin/thisroot.csh
setenv ROOTSYS /cvmfs/sft.cern.ch/lcg/app/releases/ROOT/5.34.36/x86_64-slc6-gcc48-opt/root
setenv LD_LIBRARY_PATH $ROOTSYS/lib
setenv PATH .:$ROOTSYS/bin:$PATH



#CERNlib steve wood's installation
setenv CERN /apps/cernlib/x86_64_rhel7
setenv CERN_ROOT /apps/cernlib/x86_64_rhel7
setenv CERNLIB $CERN/2005/lib
setenv CERNBIN $CERN/2005/bin
setenv PATH ${CERNBIN}:$PATH
setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:${CERNLIB}

setenv LIBRARY_PATH ${CERNLIB}


#add in LHAPDF lib
#setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:/cvmfs/eic.opensciencegrid.org/x8664_sl7/MCEG/releases/env/EIC2020b/lib/LHAPDF5
#setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:/u/group/eic/users/yxzhao/EW/lhapdf5/install/lib
#setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:/work/eic/users/xiaochao/lhapdf630/lib
setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:/work/eic/users/xiaochao/lhapdf640/lib   

setenv PYTHON /apps/python/2.7.12/bin
setenv PATH ${PYTHON}:${PATH}

#PDFsets should be under
#/cvmfs/sft.cern.ch/lcg/external/lhapdfsets/current

