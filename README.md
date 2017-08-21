# guiPMC

This repository contains the source code for a Qt user interface which works as
a frontend for goPMC (THIS REPOSITORY DOES NOT CONTAIN goPMC ITSELF. Contact the authors of [1-3] for more information regarding goPMC).

It also provides DICOM readers for CT images and Ion RT plans based on grassroots dicom (GDCM) as included with ITK.
Thus ITK with the ImageIO module is a dependency. (ITK is also used for geometric transformations of the particle sources).

This code uses Qt for the user interface and a heavily modified Qt example for displaying particle sources in 3D.
(Licensing is still unclear to me, but I think re-use of that example requires LGPL/GPL compliance)


System requirements:
OpenCL 1.1 or higher.
CPU/GPU supporting OpenCL.

NEEDED (NOT PROVIDED) goPMC files:
Subdirectories:
include:
	goPMC.h: goPMC interface.
	cl.hpp: c++ wrapper for OpenCL 1.1.
bin:
	dcmtk.dll(.lib) and libDicomRT.dll(.lib): Libraries for reading and processing Dicom CT data.
	goPMC.dll(.lib): goPMC library.
input:
	Physics input data.

Related publications:

[1] GPU-based fast Monte Carlo dose calculation for proton therapy
Xun Jia, Jan Schuemann, Harald Paganetti and Steve B. Jiang
Physics in Medicine and Biology, Volume 57, Number 23

[2] Validation of a GPU-based Monte Carlo code (gPMC) for proton radiation therapy: clinical cases study
Drosoula Giantsoudi, Jan Schuemann, Xun Jia, Stephen Dowdell, Steve B. Jiang and Harald Paganetti
Physics in Medicine and Biology, Volume 60, Number 6

[3] Recent developments and comprehensive evaluations of a GPU-based Monte Carlo package for proton therapy
Nan Qin, Pablo Botas, Drosoula Giantsoudi, Jan Schuemann, Zhen Tian, Steve B. Jiang, Harald Paganetti and Xun Jia
Accepted by Physics in Medicine and Biology


For more questions please email to andreasg@phys.au.dk
