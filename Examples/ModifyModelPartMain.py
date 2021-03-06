'''
The file can take an mdpa file and translate and rotate the node and element positions.
This can also combine 2 more mdpa files into a single one without loosing the node,element,condition details.
Transaltion is to be given in the format (x,y,z)
Rotation is to be given ([x,y,z],angle) where [x,y,z] represent the roation vector and 'angle' is the angle of rotation.
Chair of Structural Analysis, Technical University of Munich
'''

import sys
sys.path.insert(0, '../')

import model_part_manipulator_utility as mdpa_util # This is on the PYTHONPATH

# Import the applications where the elements/conditions are defined
import KratosMultiphysics.StructuralMechanicsApplication

modelpart_1 = mdpa_util.ReadModelPart("Blade", "Blade1", "StructuralMaterials.json")
modelpart_2 = mdpa_util.ReadModelPart("Blade", "Blade2","StructuralMaterials.json")
modelpart_3 = mdpa_util.ReadModelPart("Blade", "Blade3","StructuralMaterials.json")

mdpa_util.TranslateModelPart(modelpart_1, [0,0,0.65])
mdpa_util.TranslateModelPart(modelpart_2, [0,0,0.65])
mdpa_util.TranslateModelPart(modelpart_3, [0,0,0.65])

mdpa_util.RotateModelPart(modelpart_2, [1,0,0], 120)
mdpa_util.RotateModelPart(modelpart_3, [1,0,0], 240)

model_part_0 = mdpa_util.GetDefaultModelPart("Structure") # using this for a "clean" start
mdpa_util.AddModelPart(model_part_0, modelpart_1, add_as_submodelpart=True)
mdpa_util.AddModelPart(model_part_0, modelpart_2, add_as_submodelpart=False)
mdpa_util.AddModelPart(model_part_0, modelpart_3, add_as_submodelpart=False)

# mdpa_util.AddModelPart(model_part_0, modelpart_4, add_as_submodelpart=True)

mdpa_util.WriteMdpaFile(model_part_0, "CompleteRotor")
