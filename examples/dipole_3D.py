import os
import flowtopolopy as ft

def get_script_path():
    return os.path.dirname(os.path.realpath(__file__))

# Working folder and flow file
work_folder = os.path.join(get_script_path(), "dipole_3D_files")
fname = 'q_vert'
fext = '.vti'

# Topology
flowFile = os.path.join(work_folder, fname+fext)
ft.topology(flowFile=flowFile, integrationStepSize=0.1, maxNumSteps=200,
            separatrixDist=0.1)
