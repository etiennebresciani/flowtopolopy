import vtk # this is the VTK Python wrapping from Kitware
import math
import os
import numpy as np
import sys

def topology(flowFile, separatrixDist=0.1, integrationStepSize=0.1,
             maxNumSteps=100, computeSurfaces=1, excludeBoundary=0):

    filename, file_extension = os.path.splitext(flowFile)
    if file_extension == ".vti":
        reader = vtk.vtkXMLImageDataReader()
    elif file_extension == ".vtr":
        reader = vtk.vtkXMLRectilinearGridReader()
    elif file_extension == ".vtu":
        reader = vtk.vtkXMLUnstructuredGridReader()
    elif file_extension == ".vtk":
        reader = vtk.vtkUnstructuredGridReader()
    else:
        raise ValueError("File extension not recognized.")

    # Read the flow field
    reader.SetFileName(flowFile)
    reader.Update()

    # Extract the flow topology
    topology = vtk.vtkVectorFieldTopology()
    topology.SetInputData(reader.GetOutput())
    topology.SetIntegrationStepUnit(2)
    topology.SetSeparatrixDistance(separatrixDist)
    topology.SetIntegrationStepSize(integrationStepSize)
    topology.SetMaxNumSteps(maxNumSteps)
    topology.SetComputeSurfaces(computeSurfaces)
    topology.SetUseIterativeSeeding(1)
    topology.SetExcludeBoundary(excludeBoundary)
    topology.Update()
    # print(topology)

    # Output the critical points in a file
    writer = vtk.vtkXMLPolyDataWriter()
    writer.SetInputData(topology.GetOutput(0))
    writer.SetFileName(filename + "_criticalPoints" + ".vtp");
    writer.Write()

    # Output the 1D separatrices in a file
    writer = vtk.vtkXMLPolyDataWriter()
    writer.SetInputData(topology.GetOutput(1))
    writer.SetFileName(filename + "_separatrices" + ".vtp");
    # writer.SetDataModeToAscii()
    writer.Write()

    # Output the 2D separatrices in a file
    if(computeSurfaces and reader.GetOutput().GetBounds()[5]-reader.GetOutput().GetBounds()[4] > 1e-10):
        writer = vtk.vtkXMLPolyDataWriter()
        writer.SetInputData(topology.GetOutput(2))
        writer.SetFileName(filename + "_surfaces" + ".vtp");
        writer.Write()

def segmentation(flowFile, linesFile):

    filename, file_extension = os.path.splitext(flowFile)
    if (file_extension != ".vti"):
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

    # Merge exactly coincident points to avoid unnecessary shifts in merge close points
    mergeExact = vtk.vtkCleanPolyData()
    mergeExact.SetInputConnection(append.GetOutputPort())
    mergeExact.Update()

    # Merge close points to avoid problems with Delaunay triangulation
    mergeClose = vtk.vtkVoxelGrid()
    mergeClose.SetInputConnection(mergeExact.GetOutputPort())
    mergeClose.SetConfigurationStyleToManual()
    Nx = math.ceil((bounds[1]-bounds[0]) / (dist/math.sqrt(2)))
    Ny = math.ceil((bounds[3]-bounds[2]) / (dist/math.sqrt(2)))
    mergeClose.SetDivisions(Nx, Ny, 1)
    mergeClose.Update()

    # Output separatrices and boundary points in a file
    # (first convert vtkPoints into vtkVertex cells so that they can be displayed in Paraview)
    mergeCloseVertex = vtk.vtkVertexGlyphFilter()
    mergeCloseVertex.SetInputConnection(mergeClose.GetOutputPort())
    mergeCloseVertex.Update()
    writer = vtk.vtkXMLPolyDataWriter()
    writer.SetInputData(mergeCloseVertex.GetOutput())
    writer.SetFileName(filename + "_separatricesPlusBoundary" + ".vtp");
    writer.Write()

    # Tessellate
    tess = vtk.vtkDelaunay2D()
    tess.SetInputConnection(mergeClose.GetOutputPort())
    tess.Update()

    # Color via connected regions
    conn = vtk.vtkPolyDataEdgeConnectivityFilter()
    conn.SetInputConnection(tess.GetOutputPort());
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
    cellLocator = vtk.vtkCellLocator();
    cellLocator.SetDataSet(conn.GetOutput());
    cellLocator.BuildLocator();
    cellLocator.Update()

    cellCentersFilter = vtk.vtkCellCenters();
    cellCentersFilter.SetInputData(tess.GetOutput());
    cellCentersFilter.Update();

    regionId = vtk.vtkDoubleArray()
    regionId.SetNumberOfTuples(probe.GetOutput().GetNumberOfCells())
    regionId.SetName('RegionId')
    probe.GetOutput().GetCellData().AddArray(regionId)
    if conn.GetOutput().GetCellData().HasArray('RegionArea'):
        regionArea = vtk.vtkDoubleArray()
        regionArea.SetNumberOfTuples(probe.GetOutput().GetNumberOfCells())
        regionArea.SetName('RegionArea')
        probe.GetOutput().GetCellData().AddArray(regionArea)
    for c in range(probe.GetOutput().GetNumberOfCells()):
      point = cellCentersFilter.GetOutput().GetPoint(c)
      cell = cellLocator.FindCell(point)
      if cell > -1:
          regId = conn.GetOutput().GetCellData().GetArray('RegionId').GetTuple1(cell)
          regionId.SetTuple1(c, regId)
          if conn.GetOutput().GetCellData().HasArray('RegionArea'):
              regArea = conn.GetOutput().GetCellData().GetArray('RegionArea').GetTuple1(cell)
              regionArea.SetTuple1(c, regArea)

    # Output the result in a file
    writer = vtk.vtkXMLPolyDataWriter()
    writer.SetInputData(probe.GetOutput())
    writer.SetFileName(filename + "_segmentation" + ".vtp");
    writer.Write()


def segmentation_simpler(flowFile, linesFile):

    filename, file_extension = os.path.splitext(flowFile)
    if (file_extension != ".vti"):
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

    # Merge exactly coincident points to avoid unnecessary shifts in merge close points
    mergeExact = vtk.vtkCleanPolyData()
    mergeExact.SetInputConnection(append.GetOutputPort())
    mergeExact.Update()

    # Merge close points to avoid problems with Delaunay triangulation
    mergeClose = vtk.vtkVoxelGrid()
    mergeClose.SetInputConnection(mergeExact.GetOutputPort())
    mergeClose.SetConfigurationStyleToManual()
    Nx = math.ceil((bounds[1]-bounds[0]) / (dist/math.sqrt(2)))
    Ny = math.ceil((bounds[3]-bounds[2]) / (dist/math.sqrt(2)))
    mergeClose.SetDivisions(Nx, Ny, 1)
    mergeClose.Update()

    # Output separatrices and boundary points in a file
    # (first convert vtkPoints into vtkVertex cells so that they can be displayed in Paraview)
    mergeCloseVertex = vtk.vtkVertexGlyphFilter()
    mergeCloseVertex.SetInputConnection(mergeClose.GetOutputPort())
    mergeCloseVertex.Update()
    writer = vtk.vtkXMLPolyDataWriter()
    writer.SetInputData(mergeCloseVertex.GetOutput())
    writer.SetFileName(filename + "_separatricesPlusBoundary" + ".vtp");
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
    conn.SetInputConnection(tess.GetOutputPort());
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

    # Output the result in a file
    writer = vtk.vtkXMLPolyDataWriter()
    writer.SetInputData(probe.GetOutput())
    writer.SetFileName(filename + "_segmentation" + ".vtp");
    writer.Write()

def transects(segmentationFile, linesFile, tol=0.01, integrationStepSize=0.1,
              maxNumSteps=1000):

    filename, file_extension = os.path.splitext(segmentationFile)

    # Read the segmentation data
    segmentationReader = vtk.vtkXMLPolyDataReader()
    segmentationReader.SetFileName(segmentationFile)
    segmentationReader.Update()
    bounds = segmentationReader.GetOutput().GetBounds()

    # Read the lines
    linesReader = vtk.vtkXMLPolyDataReader()
    linesReader.SetFileName(linesFile)
    linesReader.Update()
    print(linesReader.GetNumberOfPoints())

    # Merge close points on lines to accelerate the process
    mergeClose = vtk.vtkVoxelGrid()
    mergeClose.SetInputConnection(linesReader.GetOutputPort())
    mergeClose.SetConfigurationStyleToManual()
    dist = tol * segmentationReader.GetOutput().GetLength() # length of bounding box diagonal
    print('dist = {}'.format(dist))
    Nx = math.ceil((bounds[1]-bounds[0])/dist)
    Ny = math.ceil((bounds[3]-bounds[2])/dist)
    mergeClose.SetDivisions(Nx, Ny, 1)
    mergeClose.Update()
    print(mergeClose.GetOutput().GetNumberOfPoints())

    # Compute orthogonal flow
    orthogonalFlow = vtk.vtkDoubleArray()
    orthogonalFlow.SetNumberOfComponents(3)
    orthogonalFlow.SetNumberOfTuples(segmentationReader.GetOutput().GetNumberOfPoints())
    orthogonalFlow.SetName('orthogonalFlow')
    segmentationReader.GetOutput().GetPointData().AddArray(orthogonalFlow)
    for p in range(segmentationReader.GetOutput().GetNumberOfPoints()):
      v = segmentationReader.GetOutput().GetPointData().GetVectors().GetTuple3(p)
      orthogonalFlow.SetTuple3(p, v[1], -v[0], 0)

    # Filter to move through the different segments
    threshold = vtk.vtkThreshold()
    threshold.SetInputData(segmentationReader.GetOutput())
    threshold.SetInputArrayToProcess(0, 0, 0, 1, 'RegionId') # (id=0 for first array, port=0, connection=0, pointData=0 and cellData=1, name)

    # Filter to compute streamlines
    tracer = vtk.vtkStreamTracer()
    tracer.SetInputData(segmentationReader.GetOutput())
    tracer.SetSourceData(mergeClose.GetOutput())
    tracer.SetIntegratorTypeToRungeKutta4()
    tracer.SetIntegrationDirectionToBoth()
    tracer.SetIntegrationStepUnit(2); #  2 --> CELL_LENGTH_UNIT
    tracer.SetInitialIntegrationStep(integrationStepSize)
    tracer.SetMaximumNumberOfSteps(maxNumSteps)
    tracer.SetMaximumPropagation(dist * maxNumSteps)
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
    if segmentationReader.GetOutput().GetCellData().HasArray('RegionArea'):
        transectsRegionArea = vtk.vtkDoubleArray()
        transectsRegionArea.SetName('RegionArea')
        transects.GetCellData().AddArray(transectsRegionArea)

    # Find longest line orthogonal to flow for each segment and store it in output
    for regId in range(int(segmentationReader.GetOutput().GetCellData().GetArray('RegionId').GetRange()[1])+1):
    #for regId in range(1):
      threshold.ThresholdBetween(regId-0.5, regId+0.5)
      threshold.Update()
      if threshold.GetOutput().GetNumberOfPoints() > 0:
        print(regId, threshold.GetOutput().GetNumberOfPoints())
        tracer.SetInputData(threshold.GetOutput())
        tracer.Update()
        longestLineId = 0
        longestLength = 0
        for c in range(tracer.GetOutput().GetNumberOfCells()):
          # ReasonForTermination must be 1 (OUT_OF_DOMAIN)
          if tracer.GetOutput().GetCellData().GetArray('ReasonForTermination').GetTuple1(c) != 1:
              continue
          polyline = tracer.GetOutput().GetCell(c)
          length = 0
          for p in range(polyline.GetPointIds().GetNumberOfIds()-1):
            p0 = np.array(tracer.GetOutput().GetPoint(polyline.GetPointId(p)))
            p1 = np.array(tracer.GetOutput().GetPoint(polyline.GetPointId(p+1)))
            length = length + np.linalg.norm(p1-p0)
          if longestLength < length:
            longestLength = length
            longestLineId = c
        if longestLength > 0:
          longestLine = tracer.GetOutput().GetCell(longestLineId)
          newLine = vtk.vtkPolyLine()
          newLine.GetPointIds().SetNumberOfIds(longestLine.GetPointIds().GetNumberOfIds() + 1)
          for p in range(longestLine.GetPointIds().GetNumberOfIds()):
            points.InsertNextPoint(tracer.GetOutput().GetPoint(longestLine.GetPointId(p)))
            newLine.GetPointIds().SetId(p, transects.GetNumberOfPoints() - 1)

          # Connect the end of the line to the boundary of the segment
          endPoint = np.array(transects.GetPoint(transects.GetNumberOfPoints() - 1))
          closestPointId = 0
          closestDist = sys.float_info.max
          for p in range(linesReader.GetOutput().GetNumberOfPoints()):
            point = linesReader.GetOutput().GetPoint(p)
            if closestDist > np.linalg.norm(endPoint - point):
              closestDist = np.linalg.norm(endPoint - point)
              closestPointId = p
          points.InsertNextPoint(linesReader.GetOutput().GetPoint(closestPointId))
          newLine.GetPointIds().SetId(newLine.GetPointIds().GetNumberOfIds() - 1, transects.GetNumberOfPoints() - 1)

          # Assign the output
          lines.InsertNextCell(newLine)
          transectsRegionId.InsertNextTuple1(regId)
          if segmentationReader.GetOutput().GetCellData().HasArray('RegionArea'):
              regArea = threshold.GetOutput().GetCellData().GetArray('RegionArea').GetTuple1(0) # simply look at the first cell
              transectsRegionArea.InsertNextTuple1(regArea)

    # Probe the velocity along the transect
    probe = vtk.vtkProbeFilter()
    probe.SetInputData(transects)
    probe.SetSourceData(segmentationReader.GetOutput())
    probe.PassCellArraysOn()
    probe.Update()

    # Output the result in a file
    writer = vtk.vtkXMLPolyDataWriter()
    writer.SetInputData(probe.GetOutput())
    writer.SetFileName(filename + "_transects" + ".vtp");
    writer.Write()
