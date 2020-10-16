import vtk # this is the VTK Python wrapping from Kitware
from vtk.util.numpy_support import vtk_to_numpy
from vtk.util.numpy_support import numpy_to_vtk
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
    writer.SetFileName(filename + "_separatricesClean" + ".vtp");
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
    writer.SetFileName(filename + "_separatricesClean" + ".vtp");
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

    # Read the separatrices
    linesReader = vtk.vtkXMLPolyDataReader()
    linesReader.SetFileName(linesFile)
    linesReader.Update()

    # Merge close points on lines to accelerate the process
    dist = tol * segmentationReader.GetOutput().GetLength() # length of bounding box diagonal
    print('dist = {}'.format(dist))
    mergeClose = vtk.vtkCleanPolyData()
    mergeClose.SetInputConnection(linesReader.GetOutputPort())
    mergeClose.ToleranceIsAbsoluteOn()
    mergeClose.SetAbsoluteTolerance(dist)
    mergeClose.Update()

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
    transectsRegionFlowRate = vtk.vtkDoubleArray()
    transectsRegionFlowRate.SetName('RegionFlowRate')
    transects.GetCellData().AddArray(transectsRegionFlowRate)
    transectsQcumul = vtk.vtkDoubleArray()
    transectsQcumul.SetName('Qcumul')
    transects.GetPointData().AddArray(transectsQcumul)
    if segmentationReader.GetOutput().GetCellData().HasArray('RegionArea'):
        transectsRegionArea = vtk.vtkDoubleArray()
        transectsRegionArea.SetName('RegionArea')
        transects.GetCellData().AddArray(transectsRegionArea)

    # Find longest line orthogonal to flow for each segment and store it in output
    for regId in range(int(segmentationReader.GetOutput().GetCellData().GetArray('RegionId').GetRange()[1])+1):
    # for regId in range(5, 6):
        threshold.ThresholdBetween(regId-0.5, regId+0.5)
        threshold.Update()
        if threshold.GetOutput().GetNumberOfPoints() > 0:
            # Orthogonal streamline tracing
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
            # writer = vtk.vtkXMLPolyDataWriter()
            # writer.SetInputData(streamlines)
            # writer.SetFileName(filename + "_regionAllTransects" + ".vtp");
            # writer.Write()
            # return

            # Find longest line or line with the largest flow rate
            longestLineId = 0
            longestLength = 0
            longestLineQ = 0
            largestQLineId = 0
            largestQ = 0
            Qall = np.empty(streamlines.GetNumberOfCells())
            Qall[:] = np.nan
            for c in range(streamlines.GetNumberOfCells()):
                # ReasonForTermination must be 1 (OUT_OF_DOMAIN)
                if streamlines.GetCellData().GetArray('ReasonForTermination').GetTuple1(c) != 1:
                    continue
                polyline = streamlines.GetCell(c)
                # Calculate transect length
                length = 0
                for p in range(polyline.GetNumberOfPoints()-1):
                    p0 = np.array(streamlines.GetPoint(polyline.GetPointId(p)))
                    p1 = np.array(streamlines.GetPoint(polyline.GetPointId(p+1)))
                    length = length + np.linalg.norm(p1-p0)
                # Discard short streamlines that may be due to seed points moved during merge
                # if length < 0.5*np.sqrt(2)*dist:
                #     continue
                if longestLength < length:
                    longestLength = length
                    longestLineId = c
                # Calculate cumulative flow rate across the streamline
                Qcumul = np.zeros(polyline.GetNumberOfPoints())
                for p in range(polyline.GetNumberOfPoints() - 1):
                    p0 = np.array(streamlines.GetPoint(polyline.GetPointId(p)))
                    p1 = np.array(streamlines.GetPoint(polyline.GetPointId(p+1)))
                    d = np.linalg.norm(p1 - p0)
                    v0 = streamlines.GetPointData().GetVectors().GetTuple3(polyline.GetPointId(p))
                    v1 = streamlines.GetPointData().GetVectors().GetTuple3(polyline.GetPointId(p+1))
                    if streamlines.GetPointData().HasArray('thickness'):
                        thick0 = streamlines.GetPointData().GetArray('thickness').GetTuple1(polyline.GetPointId(p))
                        thick1 = streamlines.GetPointData().GetArray('thickness').GetTuple1(polyline.GetPointId(p+1))
                    else:
                        thick0 = 1.
                        thick1 = 1.
                    Qcumul[p+1] = Qcumul[p] + 0.5 * d * (np.linalg.norm(v0)*thick0 + np.linalg.norm(v1)*thick1)
                Q = Qcumul[-1]
                Qall[c] = Q
                if largestQ < Q:
                    largestQ = Q
                    largestQLineId = c
                    largestQLineQcumul = Qcumul
                if longestLineId == c:
                    longestLineQ = Q
                    longestLineQcumul = Qcumul

            # Assign the output
            if longestLength > 0:
                Q_CV = np.nanstd(Qall) / np.nanmean(Qall)
                print('Q_CV = {}'.format(Q_CV))
                if Q_CV > 0.1:
                    bestLine = streamlines.GetCell(largestQLineId)
                    Q = largestQ
                    Qcumul = largestQLineQcumul
                else:
                    bestLine = streamlines.GetCell(longestLineId)
                    Q = longestLineQ
                    Qcumul = longestLineQcumul
                newLine = vtk.vtkPolyLine()
                newLine.GetPointIds().SetNumberOfIds(bestLine.GetNumberOfPoints())
                for p in range(bestLine.GetNumberOfPoints()):
                    points.InsertNextPoint(streamlines.GetPoint(bestLine.GetPointId(p)))
                    newLine.GetPointIds().SetId(p, transects.GetNumberOfPoints() - 1)
                    transectsQcumul.InsertNextTuple1(Qcumul[p])
                lines.InsertNextCell(newLine)
                transectsRegionId.InsertNextTuple1(regId)
                transectsRegionFlowRate.InsertNextTuple1(Q)
                if segmentationReader.GetOutput().GetCellData().HasArray('RegionArea'):
                    regArea = threshold.GetOutput().GetCellData().GetArray('RegionArea').GetTuple1(0) # simply look at the first cell
                    transectsRegionArea.InsertNextTuple1(regArea)

    # Add region flow rate and mean residence time info to the segmentation file
    RegionFlowRate = np.zeros([segmentationReader.GetOutput().GetNumberOfCells()])
    RegionVolume = np.zeros([segmentationReader.GetOutput().GetNumberOfCells()])
    RegionId = segmentationReader.GetOutput().GetCellData().GetArray('RegionId')
    RegionId = vtk_to_numpy(RegionId)
    if not segmentationReader.GetOutput().GetPointData().HasArray('thickness'):
        # if the input data do not have 'thickness', pretend thickness = 1
        thickness = np.ones(segmentationReader.GetOutput().GetNumberOfPoints())
        thickness = numpy_to_vtk(thickness)
        thickness.SetName('thickness')
        segmentationReader.GetOutput().GetPointData().AddArray(thickness)
    for c in range(transects.GetNumberOfCells()):
        regId = transects.GetCellData().GetArray('RegionId').GetTuple1(c)
        Q = transects.GetCellData().GetArray('RegionFlowRate').GetTuple1(c)
        RegionFlowRate[RegionId==regId] = Q
        threshold.ThresholdBetween(regId-0.5, regId+0.5)
        threshold.Update()
        integrate = vtk.vtkIntegrateAttributes()
        integrate.AddInputConnection(threshold.GetOutputPort())
        integrate.Update()
        RegionVolume[RegionId==regId] = integrate.GetOutput().GetPointData().GetArray('thickness').GetTuple1(0)
    RegionMeanResidenceTime = RegionVolume / RegionFlowRate
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
    writer.SetFileName(filename + "_transects" + ".vtp");
    writer.Write()

def flow_weighted_spacing(transectsFile, Npts=100):

    filename, file_extension = os.path.splitext(transectsFile)

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
    writer.SetFileName(filename + "FlowWeighted" + ".vtp");
    writer.Write()
