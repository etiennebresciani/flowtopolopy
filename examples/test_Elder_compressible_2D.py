import os
import flowtopolopy as ft

def get_script_path():
    return os.path.dirname(os.path.realpath(__file__))

# Working folder and flow file
work_folder = os.path.join(get_script_path(), "Elder_compressible_2D_files")
fname = 'q_vert'
fext = '.vti'

# Topology
flowFile = os.path.join(work_folder, fname+fext)
ft.topology(flowFile=flowFile, integrationStepSize=0.1, maxNumSteps=10000,
            separatrixDist=0.1)

# Segmentation
linesFile = os.path.join(work_folder, fname+'_separatrices.vtp')
ft.segmentation(flowFile=flowFile, linesFile=linesFile)
# ft.segmentation_simpler(flowFile=flowFile, linesFile=linesFile)

# Transects
segmentationFile = os.path.join(work_folder, fname+'_segmentation.vtp')
separatricesCleanFile = os.path.join(work_folder,
                                     fname+'_separatricesClean.vtp')
ft.transects(segmentationFile=segmentationFile,
             linesFile=separatricesCleanFile, resolution=0.01,
             integrationStepSize=0.5, maxNumSteps=1000, terminalSpeed=1e-3,
             output='longest', colorRegions=True)

# Flow-equally-spaced points along transects
transectsFile = os.path.join(work_folder, fname+'_segmentation_transects.vtp')
ft.flow_weighted_spacing(transectsFile=transectsFile, Npts=20)
