from sys import argv
from time import time
if len(argv) != 2:
    raise Exception("please provide an *.mdpa file")

def TryImportApplication(application_name):
    # Loading the Applications is needed bcs that is where the elements/conditions are!
    # reading the mdpa will fail if the mdpa contains elements/conditions that are defined
    # in applications that are not available
    try:
        __import__("KratosMultiphysics."+application_name)
        print("Sucessfully imported " + application_name)
    except:
        print(application_name + " is available but import failed")

def PrintTime(label, start_time):
    print(label+": {0:.{1}f} [s]".format(time()-start_time,2))

import KratosMultiphysics as KM
KM.Logger.GetDefaultOutput().SetSeverity(KM.Logger.Severity.WARNING)
import KratosMultiphysics.kratos_utilities as kratos_utils
for app_name in kratos_utils.GetListOfAvailableApplications():
    TryImportApplication(app_name)

mdpa_file_name = argv[1][:-5] # remove ".mdpa"-extension

current_model = KM.Model()
model_part = current_model.CreateModelPart(mdpa_file_name)

import_flags = KM.ModelPartIO.READ | KM.ModelPartIO.IGNORE_VARIABLES_ERROR | KM.ModelPartIO.SKIP_TIMER
start_time = time()
model_part_io = KM.ModelPartIO(mdpa_file_name, import_flags).ReadModelPart(model_part)
PrintTime("\nMdpa reading time", start_time)
print()

# ### GID ###
# start_time = time()
# gid_mode = KM.GiDPostMode.GiD_PostBinary
# multifile = KM.MultiFileFlag.MultipleFiles
# deformed_mesh_flag = KM.WriteDeformedMeshFlag.WriteUndeformed
# write_conditions = KM.WriteConditionsFlag.WriteConditions

# gid_io = KM.GidIO(mdpa_file_name, gid_mode, multifile,
#                                    deformed_mesh_flag, write_conditions)

# gid_io.InitializeMesh(0)
# gid_io.WriteMesh(model_part.GetMesh())
# gid_io.FinalizeMesh()
# PrintTime("GiD write time", start_time)
# print()

### VTK ###
start_time = time()
vtk_parameters = KM.Parameters("""{
    "output_control_type"                : "step",
    "output_sub_model_parts"             : true,
    "save_output_files_in_folder"        : true,
    "output_path"                        : ""
}""")
vtu_parameters = vtk_parameters.Clone()

vtk_parameters["output_path"].SetString("VTK_Output_"+mdpa_file_name)

vtk_io = KM.VtkOutput(model_part, vtk_parameters)
vtk_io.PrintOutput()
PrintTime("Vtk write time", start_time)
print()

### VTU ###
start_time = time()
vtu_parameters["output_path"].SetString("VTU_Output_"+mdpa_file_name)
vtu_io = KM.VtuOutput(model_part, vtu_parameters)
vtu_io.PrintOutput()
PrintTime("Vtu write time", start_time)

kratos_utils.DeleteFileIfExisting(mdpa_file_name + ".time")
