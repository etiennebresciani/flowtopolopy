import os
import flowtopolopy as ft

# Working folder and flow file
work_folder = os.path.join(".", "test_Toth_2D")
fname = 'q_vert'
fext = '.vti'

# Topology
flowFile = os.path.join(work_folder, fname+fext)
ft.topology(flowFile=flowFile, separatrixDist=0.1,
            integrationStepSize=0.1, maxNumSteps=10000)

# Segmentation
linesFile = os.path.join(work_folder, fname+'_separatrices.vtp')
ft.segmentation(flowFile=flowFile, linesFile=linesFile)
# ft.segmentation_simpler(flowFile=flowFile, linesFile=linesFile)

# Transects
segmentationFile = os.path.join(work_folder, fname+'_segmentation.vtp')
separatricesCleanFile = os.path.join(work_folder, fname+'_separatricesClean.vtp')
ft.transects(segmentationFile=segmentationFile, linesFile=separatricesCleanFile, tol=0.01,
             integrationStepSize=0.6, maxNumSteps=1000)

# Flow-equally-spaced points along transects
transectsFile = os.path.join(work_folder, fname+'_segmentation_transects.vtp')
ft.flow_weighted_spacing(transectsFile=transectsFile, Npts=100)
