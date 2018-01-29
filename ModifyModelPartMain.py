'''
The file can take an mdpa file and translate and rotate the node and element positions.
This can also combine 2 more mdpa files into a single one without loosing the node,element,condition details.
Transaltion is to be given in the format (x,y,z)
Rotation is to be given ([x,y,z],angle) where [x,y,z] represent the roation vector and 'angle' is the angle of rotation.
Chair of Structural Analysis, Technical University of Munich
'''



import model_part_manipulator_utility

modelpart_1 = model_part_manipulator_utility.ModelPartManipulator("samplestructural")

modelpart_1.TranslateModelPart([0,0,0])

modelpart_1.RotateModelPart([0,0,4], 0)


modelpart_2 = model_part_manipulator_utility.ModelPartManipulator("samplestructural.mdpa")

modelpart_2.TranslateModelPart([0,0,0])

modelpart_2.RotateModelPart([0,0,5], 120)

modelpart_3 = model_part_manipulator_utility.ModelPartManipulator("samplestructural.mdpa")

modelpart_3.TranslateModelPart([0,0,0])

modelpart_3.RotateModelPart([0,0,5], 240)


modelpart_1.AddModelPart(modelpart_2.GetModelPart())
modelpart_1.AddModelPart(modelpart_3.GetModelPart())


modelpart_1.WriteMdpaFile("result.mdpa")

#modelpart_2.WriteMdpaFile("result2.mdpa")