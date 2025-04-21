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

RUN_GUI=False

class ShapesClassification:
    def __init__(self, shapes):
        # Assumes that input shapes are volumes.
        
        gmsh.model.occ.synchronize()

        # Renames surfaces defining volumes.
        for s in shapes:
            if s[0] != 3:
                continue
            name = gmsh.model.get_entity_name(*s)
            boundaries = gmsh.model.get_boundary([s])
            for b in boundaries:
                gmsh.model.set_entity_name(*b, name)

        
        # Stores volumes.
        self.pec_volumes = self.get_entities_with_label_and_dim(shapes, "PEC_", 3)
        self.totalField_volumes = self.get_entities_with_label_and_dim(shapes, "TF_", 3)      
        self.sma_volumes = self.get_entities_with_label_and_dim(shapes, "SMA_", 3)
        
        # Stores surfaces
        self.pec_surfaces = self.get_entities_with_label_and_dim(shapes, "PEC_", 2)
        self.totalField_surfaces = self.get_entities_with_label_and_dim(shapes, "TF_", 2)
        self.sma_surfaces = self.get_entities_with_label_and_dim(shapes, "SMA_", 2)


    @staticmethod
    def getNumberFromName(entity_name: str, label: str):
        ini = entity_name.rindex(label) + len(label)
        num = int(entity_name[ini:])
        return num

    @staticmethod
    def get_entities_with_label_and_dim(entity_tags, label: str, dim: int):
        shapes = dict()
        for s in entity_tags:
            name = gmsh.model.get_entity_name(*s)
            if s[0] != dim or label not in name:
                continue
            num = ShapesClassification.getNumberFromName(name, label)
            shapes[num] = [s]

        return shapes


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
    classifiedShapes = ShapesClassification(
        gmsh.model.occ.importShapes(inputFile, highestDimOnly=False))
    
    # --- Geometry manipulation ---
    # Assumes SMA is the most external shape.
    tools = []
    tools.extend(classifiedShapes.pec_volumes[0])
    tools.extend(classifiedShapes.totalField_volumes[0])
    gmsh.model.occ.cut( classifiedShapes.sma_volumes[0], tools, removeObject=True, removeTool=False)
    gmsh.model.occ.synchronize()

    tools = []
    tools.extend(classifiedShapes.pec_volumes[0])
    gmsh.model.occ.cut( classifiedShapes.totalField_volumes[0], tools, removeObject=True, removeTool=True)
    gmsh.model.occ.synchronize()
    

    # --- Physical groups ---
    
    
    # Meshing.
    for [opt, val] in meshing_options.items():
        gmsh.option.setNumber(opt, val)
        
    gmsh.model.mesh.generate(3)

    gmsh.fltk.run() # For debugging only.
    

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
    if RUN_GUI: # For debugging only.
        gmsh.fltk.run()

    gmsh.write(case_name + '.vtk') # For debugging only

    gmsh.finalize()
