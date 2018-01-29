import KratosMultiphysics
# Loading the Applications is needed bcs that is where the elements are!
import KratosMultiphysics.StructuralMechanicsApplication as KratosStructural
import KratosMultiphysics.FluidDynamicsApplication as KratosCFD

from sys import argv
from os.path import splitext, basename
from os import remove

if len(argv) != 2:
    raise Exception("please provide an *.mdpa file")

mdpa_file_name = splitext(basename(argv[1]))[0]

model_part = KratosMultiphysics.ModelPart("MDPAToGID")

model_part_io = KratosMultiphysics.ModelPartIO(mdpa_file_name).ReadModelPart(model_part)

gid_mode = KratosMultiphysics.GiDPostMode.GiD_PostAscii
multifile = KratosMultiphysics.MultiFileFlag.MultipleFiles
deformed_mesh_flag = KratosMultiphysics.WriteDeformedMeshFlag.WriteUndeformed
write_conditions = KratosMultiphysics.WriteConditionsFlag.WriteConditions

gid_io = KratosMultiphysics.GidIO(mdpa_file_name, gid_mode, multifile,
                                   deformed_mesh_flag, write_conditions)

gid_io.InitializeMesh(0)
gid_io.WriteMesh(model_part.GetMesh())
gid_io.FinalizeMesh()

try:
    remove(mdpa_file_name + ".time")
    print("*.time file removed sucessfully")
except FileNotFoundError as e:
    print("*.time file could not be removed")