import gmsh
import os
from resources.gmshScripts.mesher_3d import *
import sys

sys.path.insert(0, '.')


dir_path = os.path.dirname(os.path.realpath(__file__)) + '/'
testdata_path = dir_path + '../testData/'


def countEntitiesInPhysicalGroupWithName(name: str):
    return len(
        gmsh.model.getEntitiesForPhysicalGroup(
            *getPhysicalGrupWithName(name)
        )
    )

def inputFileFromCaseName(case_name):
    return testdata_path + case_name + '/' + case_name + ".step"

def test_getNumberFromEntityName():
    assert (ShapesClassification.getNumberFromName(
        'Shapes/PEC_1',
        'PEC_') == 1
    )

def test_sphere_rcs_r_1m():
    runCase(testdata_path, 'sphere_rcs_r_1m')
 
def test_meshFromStep_with_sphere_rcs_r_1m():
    gmsh.initialize()
    
    case_name = 'sphere_rcs_r_1m'
    meshFromStep(inputFileFromCaseName(case_name), case_name)

    # pGNames = [gmsh.model.getPhysicalName(*pG) for pG in pGs]
    # assert ('Conductor_0' in pGNames)
    # assert ('Conductor_1' in pGNames)
    # assert ('Dielectric_1' in pGNames)
    # assert ('Vacuum' in pGNames)

    # c0pG = getPhysicalGrupWithName('Conductor_0')
    # c0ents = gmsh.model.getEntitiesForPhysicalGroup(*c0pG)
    # assert (len(c0ents) == 1)

    # # gmsh.fltk.run()  # for debugging only.
    # gmsh.finalize()   

    assert False # WIP

