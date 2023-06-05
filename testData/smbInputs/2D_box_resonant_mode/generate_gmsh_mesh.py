import gmsh
import os 

dir_path = os.path.dirname(os.path.realpath(__file__))

case_name = "2D_box_resonant_mode"

material_tag = {
    "Vacuum": 1,
    "PEC": 2
}

gmsh.initialize()
gmsh.model.add(case_name)

# Importing from FreeCAD generated steps.
gmsh.option.setString('Geometry.OCCTargetUnit', 'M')
region = gmsh.model.occ.importShapes(dir_path + '/vacuum.step')
boundary = gmsh.model.occ.importShapes(dir_path + '/outer_boundary.step')
gmsh.model.occ.synchronize()

# Geometry manipulation.
embeddedFace, embFaceMap = gmsh.model.occ.fragment(region, boundary)
gmsh.model.occ.remove(gmsh.model.occ.getEntities(1), True)
gmsh.model.occ.synchronize()

# Creates physical groups
gmsh.model.addPhysicalGroup(2, [region[0][1]], material_tag['Vacuum'])

for line in boundary:
    gmsh.model.addPhysicalGroup(1, [line[1]], material_tag['PEC'])


# Meshing.
gmsh.option.setNumber("Mesh.MeshSizeMin", 0.1)
gmsh.option.setNumber("Mesh.MeshSizeMax", 0.1)
gmsh.option.setNumber("Mesh.ElementOrder", 1)
gmsh.model.mesh.generate(2)

# Exporting
gmsh.write(dir_path + '/' + case_name + ".vtk")

gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)

gmsh.write(dir_path + '/' + case_name + ".msh")

# if '-nopopup' not in sys.argv:
#     gmsh.fltk.run()

gmsh.finalize()
