import gmsh
from itertools import chain
from pathlib import Path

DEFAULT_MESHING_OPTIONS = {
    
    "Mesh.MshFileVersion": 2.2,   # Mandatory for MFEM compatibility
    # "Mesh.MeshSizeFromCurvature": 50,
    "Mesh.ElementOrder": 1,
    "Mesh.ScalingFactor": 1e-3, # Assumes input in mm, output in m.
    "Mesh.SurfaceFaces": 1,
    "Mesh.MeshSizeMax": 250,

    "General.DrawBoundingBoxes": 1,
    "General.Axes": 1,

    "Geometry.SurfaceType": 2,    # Diplay surfaces as solids rather than dashed lines.
    # "Geometry.OCCBoundsUseStl": 1,
    # "Geometry.OCCSewFaces": 1,
    # "Geometry.Tolerance": 1e-3,
}

RUN_GUI=True

class ShapesClassification:
    def __init__(self, shapes):
        gmsh.model.occ.synchronize()

        self.allShapes = shapes

        self.pecs = self.get_surfaces_with_label(shapes, "PEC_")
        self.sma = self.get_surfaces_with_label(shapes, "SMA_")
        self.totalField = self.get_surfaces_with_label(shapes, "TF_")      

    @staticmethod
    def getNumberFromName(entity_name: str, label: str):
        ini = entity_name.rindex(label) + len(label)
        num = int(entity_name[ini:])
        return num

    @staticmethod
    def get_surfaces_with_label(entity_tags, label: str):
        surfaces = dict()
        for s in entity_tags:
            name = gmsh.model.get_entity_name(*s)
            if s[0] != 2 or label not in name:
                continue
            num = ShapesClassification.getNumberFromName(name, label)
            surfaces[num] = [s]

        return surfaces

    def buildVacuumDomain(self):       
        dom = gmsh.model.occ.cut(
            dom, surfsToRemove, removeObject=False, removeTool=False)[0]
        gmsh.model.occ.synchronize()

        return dom


def getPhysicalGrupWithName(name: str):
    pGs = gmsh.model.getPhysicalGroups()
    for pG in pGs:
        if gmsh.model.getPhysicalName(*pG) == name:
            return pG


def meshFromStep(
        inputFile: str,
        case_name: str,
        meshing_options=DEFAULT_MESHING_OPTIONS):
    gmsh.model.add(case_name)

    # Importing from FreeCAD generated steps.
    # STEP default units are mm.
    allShapes = ShapesClassification(
        gmsh.model.occ.importShapes(inputFile, highestDimOnly=False)
    )

    # --- Geometry manipulation ---
    # -- Domains
    

    # -- Boundaries
    pec_bdrs = extractBoundaries(allShapes.pecs)
    sma_bdrs = extractBoundaries(allShapes.sma)

   
    
    # Meshing.
    for [opt, val] in meshing_options.items():
        gmsh.option.setNumber(opt, val)

    gmsh.model.mesh.generate(2)
    

def runFromInput(inputFile):
    case_name = Path(inputFile).stem

    gmsh.initialize()

    meshFromStep(inputFile, case_name, DEFAULT_MESHING_OPTIONS)   

    gmsh.write(case_name + '.msh')
    gmsh.write(case_name + '.vtk')
    gmsh.finalize()

def runCase(
        folder: str,
        case_name: str,
        meshing_options=DEFAULT_MESHING_OPTIONS):

    gmsh.initialize()

    inputFile = folder + case_name + '/' + case_name + ".step"
    meshFromStep(inputFile, case_name, meshing_options)
    
    gmsh.write(case_name + '.msh')
    if RUN_GUI: # for debugging only.
        gmsh.fltk.run()

    gmsh.finalize()
