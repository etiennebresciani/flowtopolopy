import vtk # this is the VTK Python wrapping from Kitware
from vtk.util.numpy_support import vtk_to_numpy
from vtk.util.numpy_support import numpy_to_vtk
import math
import os
import numpy as np
# import timeit

def topology(flowFile, integrationStepUnit=2, integrationStepSize=0.1,
             maxNumSteps=100, separatrixDist=0.1, useIterativeSeeding=True,
             computeSurfaces=True, excludeBoundary=False, vectorName=None):
    """
    Computes critical points and separatrices of a 2D or 3D flow field. Results
    are written to the following files:\n
    *_criticalPoints.vtp contains the critical points\n        
    *_separatrices.vtp contains the 1D separatrices (i.e., lines)\n        
    *_surfaces.vtp contains the 2D separatrices (i.e., surfaces; only in a 3D
    flow field)\n
        
    Parameters
    ----------
    flowFile: string
        Path to file containing the flow field to be processed. The file can
        be .vti, .vtr, .vtu or .vtk.
    integrationStepUnit: integer
        Unit for integrationStepSize in vtkStreamTracer and for separatrixDist:\n
        1 = LENGTH_UNIT, i.e. all sizes are expresed in coordinate scale\n
        2 = CELL_LENGTH_UNIT, i.e. all sizes are expresed in cell scale\n
    integrationStepSize: float
        Initial, minimum, and maximum step size in vtkStreamTracer
        (expressed in IntegrationStepUnit).
    maxNumSteps: integer
        Maximum number of iterations in vtkStreamTracer.
    separatrixDist: float
        Distance by which the seedpoints of the separatrices are placed away
        from the saddle (expressed in IntegrationStepUnit).
    useIterativeSeeding: bool
        Specify if the simple (fast) or iterative (correct) version is called.
    computeSurfaces: bool
        Specify if the separating surfaces (separatrices in 3D) are computed
        or not.
    excludeBoundary: bool
        Specify if the boundary cells are treated or not (it may be necesarry
        to avoid cells along no-flow boundaries).
    vectorName: string
        Specify the name of the velocity vector field (if None, takes the first
        one found).

    """
    
    # Create reader depending on file extension
    fname, fext = os.path.splitext(flowFile)
    if fext == ".vti":
        reader = vtk.vtkXMLImageDataReader()
    elif fext == ".vtr":
        reader = vtk.vtkXMLRectilinearGridReader()
    elif fext == ".vtu":
        reader = vtk.vtkXMLUnstructuredGridReader()
    elif fext == ".vtk":
        reader = vtk.vtkUnstructuredGridReader()
    else:
        raise ValueError("File extension not recognized.")

    # Read the flow file
    reader.SetFileName(flowFile)
    reader.Update()
    
    # time1 = timeit.default_timer()

    # Compute flow topology
    topology = vtk.vtkVectorFieldTopology()
    topology.SetInputData(reader.GetOutput())
    topology.SetIntegrationStepUnit(integrationStepUnit)
    topology.SetSeparatrixDistance(separatrixDist)
    topology.SetIntegrationStepSize(integrationStepSize)
    topology.SetMaxNumSteps(maxNumSteps)
    topology.SetComputeSurfaces(computeSurfaces)
    topology.SetUseIterativeSeeding(useIterativeSeeding)
    topology.SetExcludeBoundary(excludeBoundary)
    if vectorName is not None:
        topology.SetInputArrayToProcess(0, 0, 0, vtk.vtkDataObject.FIELD_ASSOCIATION_POINTS, vectorName)
    topology.Update()
    # print(topology)
    
    # time2 = timeit.default_timer()
    # timeAB = time2 - time1
    # print('timeAB = {}'.format(timeAB))

    # Output the critical points in a file
    writer = vtk.vtkXMLPolyDataWriter()
    writer.SetInputData(topology.GetOutput(0))
    writer.SetFileName(fname + "_criticalPoints" + ".vtp")
    writer.Write()

    # Output the 1D separatrices in a file
    writer = vtk.vtkXMLPolyDataWriter()
    writer.SetInputData(topology.GetOutput(1))
    writer.SetFileName(fname + "_separatrices" + ".vtp")
    # writer.SetDataModeToAscii()
    writer.Write()

    # Output the 2D separatrices in a file
    isReal3D = False
    length1 = reader.GetOutput().GetBounds()[1]-reader.GetOutput().GetBounds()[0]
    length3 = reader.GetOutput().GetBounds()[5]-reader.GetOutput().GetBounds()[4]
    if(length3/length1 > 1e-10):
        isReal3D = True
    if(computeSurfaces and isReal3D):
        writer = vtk.vtkXMLPolyDataWriter()
        writer.SetInputData(topology.GetOutput(2))
        writer.SetFileName(fname + "_surfaces" + ".vtp")
        writer.Write()

def segmentation(flowFile, linesFile):
    """
    Segments a 2D flow field into different regions based on lines taken as
    boundaries. Results are written to file *_segmentation.vtp where cell data
    contain a region identifier (RegionId) and the area of the region to which
    the cell belongs (CellRegionArea). In addition, the generated file contains
    point data interpolated from the original flowField. Some intermediary
    construction files are also created.
        
    Parameters
    ----------
    flowFile: string
        Path to file containing the flow field to be processed. The file can
        only be .vti.
    linesFile: string
        Path to file containing the lines taken as boundary for the
        segmentation. The file needs to be .vtp.

    """
    
    # Check the file extension
    fname, fext = os.path.splitext(flowFile)
    if (fext != ".vti"):
        raise ValueError("The file extension must be .vti")

    # Read the segmentation line boundaries
    linesReader = vtk.vtkXMLPolyDataReader()
    linesReader.SetFileName(linesFile)
    linesReader.Update()
    separatrices = linesReader.GetOutput()
    
    # time1 = timeit.default_timer()

    # Get the four corner points from the input image
    imageReader = vtk.vtkXMLImageDataReader()
    imageReader.SetFileName(flowFile)
    imageReader.Update()
    bounds = imageReader.GetOutput().GetBounds()

    # Find the largest distance between two consecutive points
    dist = 0
    for c in range(separatrices.GetNumberOfCells()):
      cell = separatrices.GetCell(c)
      for p in range(1, cell.GetNumberOfPoints()):
        dist = max(dist, vtk.vtkMath.Distance2BetweenPoints(separatrices.GetPoint(cell.GetPointId(p)), separatrices.GetPoint(cell.GetPointId(p-1))))
    dist = math.sqrt(dist)
    print('dist = {}'.format(dist))

    # Merge close points to avoid problems with Delaunay triangulation
    mergeClose = vtk.vtkCleanPolyData()
    mergeClose.SetInputData(separatrices)
    mergeClose.ToleranceIsAbsoluteOn()
    mergeClose.SetAbsoluteTolerance(dist)
    mergeClose.Update()
    separatrices = mergeClose.GetOutput()

    # Save the "clean" separatrices in a file (useful for transects)
    writer = vtk.vtkXMLPolyDataWriter()
    writer.SetInputData(separatrices)
    fnameLines, fextLines = os.path.splitext(linesFile)
    writer.SetFileName(fnameLines + "Clean" + ".vtp")
    writer.Write()

    # Construct a bounding array of points to improve the tessellation process.
    plane = vtk.vtkPlaneSource()
    Nx = math.ceil((bounds[1]-bounds[0]) / dist)
    Ny = math.ceil((bounds[3]-bounds[2]) / dist)
    plane.SetResolution(Nx, Ny)
    plane.SetOrigin(bounds[0], bounds[2], 0.0)
    plane.SetPoint1(bounds[1], bounds[2], 0.0)
    plane.SetPoint2(bounds[0], bounds[3], 0.0)

    # Extract edges
    edges = vtk.vtkFeatureEdges()
    edges.SetInputConnection(plane.GetOutputPort())
    edges.ExtractAllEdgeTypesOff()
    edges.BoundaryEdgesOn()
    edges.Update()

    # Append separatrices and boundary points
    append = vtk.vtkAppendPolyData()
    append.AddInputData(separatrices)
    append.AddInputData(edges.GetOutput())
    append.Update()

    # Merge again close points to avoid problems with Delaunay triangulation
    mergeClose = vtk.vtkCleanPolyData()
    mergeClose.SetInputConnection(append.GetOutputPort())
    mergeClose.ToleranceIsAbsoluteOn()
    mergeClose.SetAbsoluteTolerance(0.5*dist)
    mergeClose.Update()
    
    # time2 = timeit.default_timer()
    # timeB = time2 - time1
    # print('timeB = {}'.format(timeB))

    # Output separatrices and boundary points in a file
    writer = vtk.vtkXMLPolyDataWriter()
    writer.SetInputData(mergeClose.GetOutput())
    writer.SetFileName(fnameLines + "PlusBoundary" + ".vtp")
    writer.Write()
    
    # time3 = timeit.default_timer()

    # Tessellate
    tess = vtk.vtkDelaunay2D()
    tess.SetInputConnection(mergeClose.GetOutputPort())
    tess.Update()
    
    # time4 = timeit.default_timer()
    # timeC = time4 - time3
    # print('timeC = {}'.format(timeC))

    # Color via connected regions
    conn = vtk.vtkPolyDataEdgeConnectivityFilter()
    conn.SetInputConnection(tess.GetOutputPort())
    #conn.SetSourceConnection(edges.GetOutputPort())
    conn.BarrierEdgesOn()
    maxBarrierEdgeLength = 2*dist
    conn.SetBarrierEdgeLength(0.0, maxBarrierEdgeLength)
    conn.SetExtractionModeToAllRegions()
    conn.GrowLargeRegionsOn()
    threshold = maxBarrierEdgeLength**2 / ((bounds[1]-bounds[0])*(bounds[3]-bounds[2]))
    print('threshold = {}'.format(threshold))
    conn.SetLargeRegionThreshold(threshold)
    conn.ColorRegionsOn()
    conn.CellRegionAreasOn()
    conn.Update()
    
    # time5 = timeit.default_timer()
    # timeD = time5 - time4
    # print('timeD = {}'.format(timeD))

    # Transform input data to polydata for append to work
    geometry = vtk.vtkGeometryFilter()
    geometry.SetInputData(imageReader.GetOutput())
    geometry.Update()

    # Produce combined dataset of original data and separatrices
    append = vtk.vtkAppendPolyData()
    append.AddInputData(geometry.GetOutput())
    append.AddInputData(conn.GetOutput())
    append.Update()

    # Tessellate
    tess = vtk.vtkDelaunay2D()
    tess.SetInputConnection(append.GetOutputPort())
    tess.Update()

    # Interpolate the arrays of the input data to the combined data
    probe = vtk.vtkProbeFilter()
    probe.SetInputData(tess.GetOutput())
    probe.SetSourceData(imageReader.GetOutput())
    probe.Update()

    # Transfer the RegionId to the combined data by looking up the cell centers
    cellLocator = vtk.vtkCellLocator()
    cellLocator.SetDataSet(conn.GetOutput())
    cellLocator.BuildLocator()
    cellLocator.Update()

    cellCentersFilter = vtk.vtkCellCenters()
    cellCentersFilter.SetInputData(tess.GetOutput())
    cellCentersFilter.Update()

    regionId = vtk.vtkDoubleArray()
    regionId.SetNumberOfTuples(probe.GetOutput().GetNumberOfCells())
    regionId.SetName('RegionId')
    probe.GetOutput().GetCellData().AddArray(regionId)
    if conn.GetOutput().GetCellData().HasArray('CellRegionArea'):
        regionArea = vtk.vtkDoubleArray()
        regionArea.SetNumberOfTuples(probe.GetOutput().GetNumberOfCells())
        regionArea.SetName('CellRegionArea')
        probe.GetOutput().GetCellData().AddArray(regionArea)
    for c in range(probe.GetOutput().GetNumberOfCells()):
      point = cellCentersFilter.GetOutput().GetPoint(c)
      cell = cellLocator.FindCell(point)
      if cell > -1:
          regId = conn.GetOutput().GetCellData().GetArray('RegionId').GetTuple1(cell)
          regionId.SetTuple1(c, regId)
          if conn.GetOutput().GetCellData().HasArray('CellRegionArea'):
              regArea = conn.GetOutput().GetCellData().GetArray('CellRegionArea').GetTuple1(cell)
              regionArea.SetTuple1(c, regArea)
    
    # time6 = timeit.default_timer()
    # timeE = time6 - time5
    # print('timeE = {}'.format(timeE))

    # Remove 'vtkValidPointMask' array, which otherwise causes problems later
    # (when probing the velocity along the streamlines in transects function)
    probe.GetOutput().GetPointData().RemoveArray('vtkValidPointMask')

    # Output the result in a file
    writer = vtk.vtkXMLPolyDataWriter()
    writer.SetInputData(probe.GetOutput())
    writer.SetFileName(fname + "_segmentation" + ".vtp")
    writer.Write()

def segmentation_simpler(flowFile, linesFile):
    """
    Segments a 2D flow field into different regions based on lines taken as
    boundaries. Results are written to file *_segmentation.vtp where cell data
    contain a region identifier (RegionId) and the area of the region to which
    the cell belongs (CellRegionArea). In addition, the generated file contains
    point data interpolated from the original flowField. Some intermediary
    construction files are also created.
    
    This algorithm is somewhat simpler than in the "segmentation" function but
    should give similar results.
        
    Parameters
    ----------
    flowFile: string
        Path to file containing the flow field to be processed. The file can
        only be .vti.
    linesFile: string
        Path to file containing the lines taken as boundary for the
        segmentation. The file needs to be .vtp.

    """
    
    # Check the file extension
    fname, fext = os.path.splitext(flowFile)
    if (fext != ".vti"):
        raise ValueError("The file extension must be .vti")

    # Read the segmentation line boundaries
    linesReader = vtk.vtkXMLPolyDataReader()
    linesReader.SetFileName(linesFile)
    linesReader.Update()
    separatrices = linesReader.GetOutput()

    # Get the four corner points from the input image
    imageReader = vtk.vtkXMLImageDataReader()
    imageReader.SetFileName(flowFile)
    imageReader.Update()
    bounds = imageReader.GetOutput().GetBounds()

    # Find largest distance between two consecutive points
    dist = 0
    for c in range(separatrices.GetNumberOfCells()):
      cell = separatrices.GetCell(c)
      for p in range(1, cell.GetNumberOfPoints()):
        dist = max(dist, vtk.vtkMath.Distance2BetweenPoints(separatrices.GetPoint(cell.GetPointId(p)), separatrices.GetPoint(cell.GetPointId(p-1))))
    dist = math.sqrt(dist)
    print('dist = {}'.format(dist))

    # Merge close points to avoid problems with Delaunay triangulation
    mergeClose = vtk.vtkCleanPolyData()
    mergeClose.SetInputData(separatrices)
    mergeClose.ToleranceIsAbsoluteOn()
    mergeClose.SetAbsoluteTolerance(dist)
    mergeClose.Update()
    separatrices = mergeClose.GetOutput()

    # Save the "clean" separatrices in a file (useful for transects)
    writer = vtk.vtkXMLPolyDataWriter()
    writer.SetInputData(separatrices)
    fnameLines, fextLines = os.path.splitext(linesFile)
    writer.SetFileName(fnameLines + "Clean" + ".vtp")
    writer.Write()

    # Construct a bounding array of points to improve the tessellation process.
    plane = vtk.vtkPlaneSource()
    Nx = math.ceil((bounds[1]-bounds[0]) / dist)
    Ny = math.ceil((bounds[3]-bounds[2]) / dist)
    plane.SetResolution(Nx, Ny)
    plane.SetOrigin(bounds[0], bounds[2], 0.0)
    plane.SetPoint1(bounds[1], bounds[2], 0.0)
    plane.SetPoint2(bounds[0], bounds[3], 0.0)

    # Extract edges
    edges = vtk.vtkFeatureEdges()
    edges.SetInputConnection(plane.GetOutputPort())
    edges.ExtractAllEdgeTypesOff()
    edges.BoundaryEdgesOn()
    edges.Update()

    # Append separatrices and boundary points
    append = vtk.vtkAppendPolyData()
    append.AddInputData(separatrices)
    append.AddInputData(edges.GetOutput())
    append.Update()

    # Merge again close points to avoid problems with Delaunay triangulation
    mergeClose = vtk.vtkCleanPolyData()
    mergeClose.SetInputConnection(append.GetOutputPort())
    mergeClose.ToleranceIsAbsoluteOn()
    mergeClose.SetAbsoluteTolerance(0.5*dist)
    mergeClose.Update()

    # Output separatrices and boundary points in a file
    writer = vtk.vtkXMLPolyDataWriter()
    writer.SetInputData(mergeClose.GetOutput())
    writer.SetFileName(fnameLines + "PlusBoundary" + ".vtp")
    writer.Write()

    # Extract mesh
    mesh = vtk.vtkGeometryFilter()
    mesh.SetInputData(imageReader.GetOutput())
    mesh.Update()

    # Append mesh
    append = vtk.vtkAppendPolyData()
    append.AddInputData(mergeClose.GetOutput())
    append.AddInputData(mesh.GetOutput())
    append.Update()

    # Tessellate
    tess = vtk.vtkDelaunay2D()
    # tess.SetInputConnection(mergeClose.GetOutputPort())
    tess.SetInputConnection(append.GetOutputPort())
    tess.Update()

    # Color via connected regions
    conn = vtk.vtkPolyDataEdgeConnectivityFilter()
    conn.SetInputConnection(tess.GetOutputPort())
    #conn.SetSourceConnection(edges.GetOutputPort())
    conn.BarrierEdgesOn()
    maxBarrierEdgeLength = 2*dist
    conn.SetBarrierEdgeLength(0.0, maxBarrierEdgeLength)
    conn.SetExtractionModeToAllRegions()
    conn.GrowLargeRegionsOn()
    threshold = maxBarrierEdgeLength**2 / ((bounds[1]-bounds[0])*(bounds[3]-bounds[2]))
    print('threshold = {}'.format(threshold))
    conn.SetLargeRegionThreshold(threshold)
    conn.ColorRegionsOn()
    conn.Update()

    # Interpolate the arrays of the input data to the combined data
    probe = vtk.vtkProbeFilter()
    probe.SetInputData(conn.GetOutput())
    probe.SetSourceData(imageReader.GetOutput())
    probe.PassCellArraysOn()
    probe.Update()
    
    # Remove 'vtkValidPointMask' array, which otherwise causes problems later
    # (when probing the velocity along the streamlines in transects function)
    probe.GetOutput().GetPointData().RemoveArray('vtkValidPointMask')

    # Output the result in a file
    writer = vtk.vtkXMLPolyDataWriter()
    writer.SetInputData(probe.GetOutput())
    writer.SetFileName(fname + "_segmentation" + ".vtp")
    writer.Write()

def transects(segmentationFile, linesFile, resolution=0.01,
              integrationStepUnit=2, integrationStepSize=0.1, maxNumSteps=1000,
              terminalSpeed=1e-12, output='largestQ', colorRegions=False):
    """
    Computes transects orthogonal to flow in the regions of a segmented 2D flow
    field and the flow rate accross them. Results are written to file
    *_transects.vtp.
    
    A "thickness" (i.e., a length in the 3rd dimension) can be specified
    as a point data array having this name in the segmentationFile. If not
    found, "thickness" is assumed to be uniformly equal to 1.
        
    Parameters
    ----------
    segmentationFile: string
        Path to file containing the segmented flow field to be processed
        (result of the "segmentation" function).
    linesFile: string
        Path to file containing the lines taken as boundary for the
        segmentation after being "cleaned" in the segmentation function.
    resolution: float
        Coefficient (0<resolution<1) used to limit the number of generated
        transects. A large resolution implies less transects.
    integrationStepUnit: integer
        Unit for integrationStepSize in vtkStreamTracer:\n
        1 = LENGTH_UNIT, i.e. all sizes are expresed in coordinate scale\n
        2 = CELL_LENGTH_UNIT, i.e. all sizes are expresed in cell scale\n
    integrationStepSize: float
        Initial, minimum, and maximum step size in vtkStreamTracer
        (expressed in IntegrationStepUnit).
    maxNumSteps: integer
        Maximum number of iterations in vtkStreamTracer.
    terminalSpeed: float
        Terminal speep in vtkStreamTracer.
    output: string
        Controls the output:\n
        'all' = all generated transects are output.\n
        'longest' = only the longest transect of each region is output.\n
        'largestQ' = only the transect bearing the largest flow rate in each
        region is output.\n
    colorRegions: bool
        If true, RegionFlowRate and RegionMeanResidenceTime are added to the
        cell data of the segmentation file.

    """

    fname, fext = os.path.splitext(segmentationFile)

    # Read the segmentation data
    segmentationReader = vtk.vtkXMLPolyDataReader()
    segmentationReader.SetFileName(segmentationFile)
    segmentationReader.Update()

    # Read the separatrices
    linesReader = vtk.vtkXMLPolyDataReader()
    linesReader.SetFileName(linesFile)
    linesReader.Update()

    # Merge close points on lines to accelerate the process
    diagonal = segmentationReader.GetOutput().GetLength() # length of bounding box diagonal
    dist = resolution * diagonal
    print('dist = {}'.format(dist))
    mergeClose = vtk.vtkCleanPolyData()
    mergeClose.SetInputConnection(linesReader.GetOutputPort())
    mergeClose.ToleranceIsAbsoluteOn()
    mergeClose.SetAbsoluteTolerance(dist)
    mergeClose.Update()
    # mergeClose = linesReader

    # Compute orthogonal flow
    orthogonalFlow = vtk.vtkDoubleArray()
    orthogonalFlow.SetNumberOfComponents(3)
    orthogonalFlow.SetNumberOfTuples(segmentationReader.GetOutput().GetNumberOfPoints())
    orthogonalFlow.SetName('orthogonalFlow')
    segmentationReader.GetOutput().GetPointData().AddArray(orthogonalFlow)
    for p in range(segmentationReader.GetOutput().GetNumberOfPoints()):
        v = segmentationReader.GetOutput().GetPointData().GetVectors().GetTuple3(p)
        orthogonalFlow.SetTuple3(p, v[1], -v[0], 0)

    # Filter to move through the different regions
    threshold = vtk.vtkThreshold()
    threshold.SetInputData(segmentationReader.GetOutput())
    threshold.SetInputArrayToProcess(0, 0, 0, 1, 'RegionId') # (id=0 for first array, port=0, connection=0, pointData=0 and cellData=1, name)

    # Filter to compute streamlines
    tracer = vtk.vtkStreamTracer()
    tracer.SetSourceData(mergeClose.GetOutput())
    tracer.SetIntegratorTypeToRungeKutta4()
    tracer.SetIntegrationDirectionToBoth()
    tracer.SetIntegrationStepUnit(integrationStepUnit)
    tracer.SetInitialIntegrationStep(integrationStepSize)
    tracer.SetMaximumNumberOfSteps(maxNumSteps)
    tracer.SetMaximumPropagation(dist * maxNumSteps)
    tracer.SetTerminalSpeed(terminalSpeed)
    tracer.SetComputeVorticity(0)
    tracer.SetInputArrayToProcess(0, 0, 0, 0, 'orthogonalFlow')

    # Generate output data
    transects = vtk.vtkPolyData()
    points = vtk.vtkPoints()
    lines = vtk.vtkCellArray()
    transects.SetPoints(points)
    transects.SetLines(lines)
    transectsRegionId = vtk.vtkDoubleArray()
    transectsRegionId.SetName('RegionId')
    transects.GetCellData().AddArray(transectsRegionId)
    transectsRegionFlowRate = vtk.vtkDoubleArray()
    transectsRegionFlowRate.SetName('RegionFlowRate')
    transects.GetCellData().AddArray(transectsRegionFlowRate)
    transectsQcumul = vtk.vtkDoubleArray()
    transectsQcumul.SetName('Qcumul')
    transects.GetPointData().AddArray(transectsQcumul)
    if segmentationReader.GetOutput().GetCellData().HasArray('CellRegionArea'):
        transectsRegionArea = vtk.vtkDoubleArray()
        transectsRegionArea.SetName('CellRegionArea')
        transects.GetCellData().AddArray(transectsRegionArea)

    # Calculate transects and flow across them in each region
    for regId in range(int(segmentationReader.GetOutput().GetCellData().GetArray('RegionId').GetRange()[1])+1):
        threshold.ThresholdBetween(regId-0.5, regId+0.5)
        threshold.Update()
        if threshold.GetOutput().GetNumberOfPoints() > 0:
            # Orthogonal streamline tracing in the region
            print(regId, threshold.GetOutput().GetNumberOfPoints())
            tracer.SetInputData(threshold.GetOutput())
            tracer.Update()

            # Probe the velocity along the streamlines
            probe = vtk.vtkProbeFilter()
            probe.SetInputData(tracer.GetOutput())
            probe.SetSourceData(segmentationReader.GetOutput())
            probe.PassCellArraysOn()
            probe.Update()
            streamlines = probe.GetOutput()

            # Variables for finding longest line or line with the largest flow rate
            longestLineId = 0
            longestLength = 0
            longestLineQ = 0
            largestQLineId = 0
            largestQ = 0
            lengthall = np.empty(streamlines.GetNumberOfCells())
            lengthall[:] = np.nan
            Qall = np.empty(streamlines.GetNumberOfCells())
            Qall[:] = np.nan
            Qcumulall = ()
            
            # Loop through all transects
            outputID = -1
            for c in range(streamlines.GetNumberOfCells()):
                streamline = streamlines.GetCell(c)
                # ReasonForTermination must be 1 (OUT_OF_DOMAIN) or 6 (STAGNATION)
                termination = streamlines.GetCellData().GetArray('ReasonForTermination').GetTuple1(c)
                if termination not in (1, 6) or streamline.GetNumberOfPoints() < 2:
                    continue
                # Calculate transect length
                length = 0
                for p in range(streamline.GetNumberOfPoints()-1):
                    p0 = np.array(streamlines.GetPoint(streamline.GetPointId(p)))
                    p1 = np.array(streamlines.GetPoint(streamline.GetPointId(p+1)))
                    length = length + np.linalg.norm(p1-p0)
                lengthall[c] = length
                if longestLength < length:
                    longestLength = length
                    longestLineId = c
                
                # Calculate cumulative flow rate across the streamline
                Qcumul = np.zeros(streamline.GetNumberOfPoints())
                for p in range(streamline.GetNumberOfPoints() - 1):
                    p0 = np.array(streamlines.GetPoint(streamline.GetPointId(p)))
                    p1 = np.array(streamlines.GetPoint(streamline.GetPointId(p+1)))
                    d = np.linalg.norm(p1 - p0)
                    v0 = streamlines.GetPointData().GetVectors().GetTuple3(streamline.GetPointId(p))
                    v1 = streamlines.GetPointData().GetVectors().GetTuple3(streamline.GetPointId(p+1))
                    if streamlines.GetPointData().HasArray('thickness'):
                        thick0 = streamlines.GetPointData().GetArray('thickness').GetTuple1(streamline.GetPointId(p))
                        thick1 = streamlines.GetPointData().GetArray('thickness').GetTuple1(streamline.GetPointId(p+1))
                    else:
                        thick0 = 1.
                        thick1 = 1.
                    Qcumul[p+1] = Qcumul[p] + 0.5 * d * (np.linalg.norm(v0)*thick0 + np.linalg.norm(v1)*thick1)
                Q = Qcumul[-1]
                Qall[c] = Q
                Qcumulall = Qcumulall + (Qcumul,)
                if largestQ < Q:
                    largestQ = Q
                    largestQLineId = c
                    largestQLineQcumul = Qcumul
                if longestLineId == c:
                    longestLineQ = Q
                    longestLineQcumul = Qcumul
                
                # Assign the output
                if output == 'all':
                    outputID = outputID + 1
                    newLine = vtk.vtkPolyLine()
                    newLine.GetPointIds().SetNumberOfIds(streamline.GetNumberOfPoints())
                    for p in range(streamline.GetNumberOfPoints()):
                        points.InsertNextPoint(streamlines.GetPoint(streamline.GetPointId(p)))
                        newLine.GetPointIds().SetId(p, transects.GetNumberOfPoints() - 1)
                        transectsQcumul.InsertNextTuple1(Qcumul[p])
                    lines.InsertNextCell(newLine)
                    transectsRegionId.InsertNextTuple1(regId)
                    transectsRegionFlowRate.InsertNextTuple1(Q)
                    if segmentationReader.GetOutput().GetCellData().HasArray('CellRegionArea'):
                        regArea = threshold.GetOutput().GetCellData().GetArray('CellRegionArea').GetTuple1(0) # simply look at the first cell
                        transectsRegionArea.InsertNextTuple1(regArea)

            # Assign the output
            if output == 'longest' or output == 'largestQ':                
                if output == 'longest':
                    bestLineID = longestLineId
                else:
                    bestLineID = largestQLineId
                bestLine = streamlines.GetCell(bestLineID)
                Q = Qall [bestLineID]
                Qcumul = Qcumulall[bestLineID]                
                newLine = vtk.vtkPolyLine()
                newLine.GetPointIds().SetNumberOfIds(bestLine.GetNumberOfPoints())
                for p in range(bestLine.GetNumberOfPoints()):
                    points.InsertNextPoint(streamlines.GetPoint(bestLine.GetPointId(p)))
                    newLine.GetPointIds().SetId(p, transects.GetNumberOfPoints() - 1)
                    transectsQcumul.InsertNextTuple1(Qcumul[p])
                lines.InsertNextCell(newLine)
                transectsRegionId.InsertNextTuple1(regId)
                transectsRegionFlowRate.InsertNextTuple1(Q)
                if segmentationReader.GetOutput().GetCellData().HasArray('CellRegionArea'):
                    regArea = threshold.GetOutput().GetCellData().GetArray('CellRegionArea').GetTuple1(0) # simply look at the first cell
                    transectsRegionArea.InsertNextTuple1(regArea)

    # Add region flow rate and mean residence time info to the segmentation file
    if colorRegions:
        RegionFlowRate = np.zeros([segmentationReader.GetOutput().GetNumberOfCells()])
        RegionVolume = np.zeros([segmentationReader.GetOutput().GetNumberOfCells()])
        RegionId = segmentationReader.GetOutput().GetCellData().GetArray('RegionId')
        RegionId = vtk_to_numpy(RegionId)
        
        # if the input data do not have 'thickness', pretend thickness = 1
        if not segmentationReader.GetOutput().GetPointData().HasArray('thickness'):
            thickness = np.ones(segmentationReader.GetOutput().GetNumberOfPoints())
            thickness = numpy_to_vtk(thickness)
            thickness.SetName('thickness')
            segmentationReader.GetOutput().GetPointData().AddArray(thickness)
        
        # Calculate average Q across transects in each region (relevant when output == 'all')
        Ntransects = np.zeros([segmentationReader.GetOutput().GetNumberOfCells()])
        for c in range(transects.GetNumberOfCells()):
            regId = transects.GetCellData().GetArray('RegionId').GetTuple1(c)
            Q = transects.GetCellData().GetArray('RegionFlowRate').GetTuple1(c)
            RegionFlowRate[RegionId==regId] = RegionFlowRate[RegionId==regId] + Q
            Ntransects[RegionId==regId] = Ntransects[RegionId==regId] + 1
        RegionFlowRate = RegionFlowRate / Ntransects
        
        # Calculate residence time in each region
        for regId in range(int(segmentationReader.GetOutput().GetCellData().GetArray('RegionId').GetRange()[1])+1):
            threshold.ThresholdBetween(regId-0.5, regId+0.5)
            if threshold.GetOutput().GetNumberOfPoints() > 0:
                threshold.Update()
                integrate = vtk.vtkIntegrateAttributes()
                integrate.AddInputConnection(threshold.GetOutputPort())
                integrate.Update()
                RegionVolume[RegionId==regId] = integrate.GetOutput().GetPointData().GetArray('thickness').GetTuple1(0)
        RegionMeanResidenceTime = RegionVolume / RegionFlowRate
        
        # Assign output
        RegionFlowRate = numpy_to_vtk(RegionFlowRate)
        RegionFlowRate.SetName('RegionFlowRate')
        segmentationReader.GetOutput().GetCellData().AddArray(RegionFlowRate)
        RegionVolume = numpy_to_vtk(RegionVolume)
        RegionVolume.SetName('RegionVolume')
        segmentationReader.GetOutput().GetCellData().AddArray(RegionVolume)
        RegionMeanResidenceTime = numpy_to_vtk(RegionMeanResidenceTime)
        RegionMeanResidenceTime.SetName('RegionMeanResidenceTime')
        segmentationReader.GetOutput().GetCellData().AddArray(RegionMeanResidenceTime)
    
        # Overwrite the segmentation file
        writer = vtk.vtkXMLPolyDataWriter()
        writer.SetInputData(segmentationReader.GetOutput())
        writer.SetFileName(segmentationFile)
        writer.Write()

    # Probe information of the original data along the transects
    probe = vtk.vtkProbeFilter()
    probe.SetInputData(transects)
    probe.SetSourceData(segmentationReader.GetOutput())
    probe.PassCellArraysOn()
    probe.PassPointArraysOn()
    probe.Update()
    transects = probe.GetOutput()

    # Output the transects in a file
    writer = vtk.vtkXMLPolyDataWriter()
    writer.SetInputData(transects)
    writer.SetFileName(fname + "_transects" + ".vtp")
    writer.Write()

def flow_weighted_spacing(transectsFile, Npts=100):
    """
    Generates points separated by an equal among of flow along transects.
    Results are written to file *FlowWeighted.vtp.

    Parameters
    ----------
    transectsFile: string
        Path to file containing the orthogonal transects (result of the
        "transects" function).
    Npts: integer
        Desired total number of points to be generated. The actual number of
        points generated can be slightly different.

    """

    fname, fext = os.path.splitext(transectsFile)

    # Read the transects
    reader = vtk.vtkXMLPolyDataReader()
    reader.SetFileName(transectsFile)
    reader.Update()
    transects = reader.GetOutput()

    # Calculate the sum of flow rates of all transects
    RegionFlowRate = vtk_to_numpy(transects.GetCellData().GetArray('RegionFlowRate'))
    Qtot = np.sum(RegionFlowRate)

    # Calculate the flow rate between two consecutive points
    Qinter = Qtot / Npts

    # Prepare output
    transectsFlowWeighted = vtk.vtkPolyData()
    points = vtk.vtkPoints()
    lines = vtk.vtkCellArray()
    transectsFlowWeighted.SetPoints(points)
    transectsFlowWeighted.SetVerts(lines)

    # Find points separated by Qinter along each transect
    for c in range(transects.GetNumberOfCells()):
        cell = transects.GetCell(c)
        pointIds = cell.GetPointIds()
        Qcumul = vtk.vtkDoubleArray()
        Qcumul.SetNumberOfTuples(cell.GetNumberOfPoints())
        transects.GetPointData().GetArray('Qcumul').GetTuples(pointIds, Qcumul)
        Qcumul = vtk_to_numpy(Qcumul)
        Q = 0.5 * Qinter
        newLine = vtk.vtkPolyLine()
        while(Q < Qcumul[-1]):
            i1 = np.searchsorted(Qcumul, Q)
            if i1 == len(Qcumul):
                break
            i0 = i1 - 1
            p0 = np.array(transects.GetPoint(cell.GetPointId(i0)))
            p1 = np.array(transects.GetPoint(cell.GetPointId(i1)))
            ratio_from_p0 = (Q-Qcumul[i0]) / (Qcumul[i1]-Qcumul[i0])
            p = p0 + ratio_from_p0*(p1-p0)
            points.InsertNextPoint(p)
            newLine.GetPointIds().InsertNextId(transectsFlowWeighted.GetNumberOfPoints() - 1)
            Q = Q + Qinter
        if newLine.GetNumberOfPoints() > 0:
            lines.InsertNextCell(newLine)

    # Output the new transects in a file
    writer = vtk.vtkXMLPolyDataWriter()
    writer.SetInputData(transectsFlowWeighted)
    writer.SetFileName(fname + "FlowWeighted" + ".vtp")
    writer.Write()

def streamlines_all_cells(flowFile, integrationStepUnit=2,
                          integrationStepSize=0.1, maxNumSteps=100,
                          computeSurfaces=1, excludeBoundary=0,
                          vectorName=None):
    """
    Computes streamlines from all cell centers of a 2D or 3D flow field.
    
    RESULTS ARE NOT WRITTEN
        
    Parameters
    ----------
    flowFile: string
        Path to file containing the flow field to be processed. The file can
        be .vti, .vtr, .vtu or .vtk.
    integrationStepUnit: integer
        Unit for integrationStepSize in vtkStreamTracer and for separatrixDist:\n
        1 = LENGTH_UNIT, i.e. all sizes are expresed in coordinate scale\n
        2 = CELL_LENGTH_UNIT, i.e. all sizes are expresed in cell scale\n
    integrationStepSize: float
        Initial, minimum, and maximum step size in vtkStreamTracer
        (expressed in IntegrationStepUnit).
    maxNumSteps: integer
        Maximum number of iterations in vtkStreamTracer.

    """
    
    # Create reader depending on file extension
    fname, fext = os.path.splitext(flowFile)
    if fext == ".vti":
        reader = vtk.vtkXMLImageDataReader()
    elif fext == ".vtr":
        reader = vtk.vtkXMLRectilinearGridReader()
    elif fext == ".vtu":
        reader = vtk.vtkXMLUnstructuredGridReader()
    elif fext == ".vtk":
        reader = vtk.vtkUnstructuredGridReader()
    else:
        raise ValueError("File extension not recognized.")

    # Read the flow field
    reader.SetFileName(flowFile)
    reader.Update()
    
    # Get the cell centers
    cellCenters = vtk.vtkCellCenters()
    cellCenters.SetInputConnection(reader.GetOutputPort())
    cellCenters.Update()
    
    # Filter to compute streamlines
    tracer = vtk.vtkStreamTracer()
    tracer.SetInputData(reader.GetOutput())
    tracer.SetSourceData(cellCenters.GetOutput())
    tracer.SetIntegratorTypeToRungeKutta4()
    tracer.SetIntegrationDirectionToBoth()
    tracer.SetIntegrationStepUnit(integrationStepUnit)
    tracer.SetInitialIntegrationStep(integrationStepSize)
    tracer.SetMaximumNumberOfSteps(maxNumSteps)
    # tracer.SetMaximumPropagation(dist * maxNumSteps)
    # tracer.SetTerminalSpeed(terminalSpeed)
    tracer.SetComputeVorticity(0)
    if vectorName is not None:
        tracer.SetInputArrayToProcess(0, 0, 0, vtk.vtkDataObject.FIELD_ASSOCIATION_POINTS, vectorName)
    tracer.Update()

    # # Output the streamlines in a file
    # writer = vtk.vtkXMLPolyDataWriter()
    # writer.SetInputData(tracer.GetOutput(0))
    # writer.SetFileName(fname + "_streamlinesAllCells" + ".vtp")
    # writer.Write()