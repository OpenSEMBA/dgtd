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

def test_sphere_rcs_r_1m():
    runCase(testdata_path, '3D_RCS_Sphere_O1_r1m')
 