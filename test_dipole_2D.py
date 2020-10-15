import os
import flowtopolopy as ft

# Folder
work_folder = os.path.join(".", "test_dipole_2D")

# Topology
flowFile = os.path.join(work_folder, 'q_vert.vti')
# ft.topology(flowFile=flowFile, separatrixDist=0.1,
#             integrationStepSize=0.1, maxNumSteps=10000)

# Segmentation
linesFile = os.path.join(work_folder, 'q_vert_separatrices.vtp')
ft.segmentation(flowFile=flowFile, linesFile=linesFile)
# ft.segmentation_simpler(flowFile=flowFile, linesFile=linesFile)

# Transects
segmentationFile = os.path.join(work_folder, 'q_vert_segmentation.vtp')
separatricesCleanFile = os.path.join(work_folder, 'q_vert_separatricesClean.vtp')
ft.transects(segmentationFile=segmentationFile, linesFile=separatricesCleanFile, tol=0.04,
             integrationStepSize=0.5)
