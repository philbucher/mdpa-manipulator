'''
The file can take an mdpa file and translate and rotate the node and element positions.
This can also combine 2 more mdpa files into a single one without loosing the node,element,condition details.
Transaltion is to be given in the format (x,y,z)
Rotation is to be given ([x,y,z],angle) where [x,y,z] represent the roation vector and 'angle' is the angle of rotation.
Chair of Structural Analysis, Technical University of Munich
'''

import model_part_manipulator_utility

modelpart_1 = model_part_manipulator_utility.ModelPartManipulator("Blade")

modelpart_1.TranslateModelPart([0,0,0.65])

modelpart_1.RotateModelPart([1,0,0], 0)


modelpart_2 = model_part_manipulator_utility.ModelPartManipulator("Blade.mdpa")

modelpart_2.TranslateModelPart([0,0,0.65])

modelpart_2.RotateModelPart([1,0,0], 120)

modelpart_3 = model_part_manipulator_utility.ModelPartManipulator("Blade.mdpa")

modelpart_3.TranslateModelPart([0,0,0.65])

modelpart_3.RotateModelPart([1,0,0], 240)

modelpart_1.AddModelPart(modelpart_2.GetModelPart())
modelpart_1.AddModelPart(modelpart_3.GetModelPart())

modelpart_1.WriteMdpaFile("CompleteRotor.mdpa")
