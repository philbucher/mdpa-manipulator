'''
The file can take an mdpa file and translate and rotate the node and element positions.
This can also combine 2 or more mdpa files into a single one without loosing the node, element, condition details.
Transaltion is to be given in the format (x,y,z)
Rotation is to be given ([x,y,z],angle) where [x,y,z] represent the roation vector and 'angle' is the angle of rotation.
Chair of Structural Analysis, Technical University of Munich
'''

import KratosMultiphysics
import numpy as np
import math
import os
import time

def ReadModelPart(mdpa_file_name):
    if mdpa_file_name.endswith('.mdpa'):
        mdpa_file_name = mdpa_file_name[:-5]
    model_part = KratosMultiphysics.ModelPart(mdpa_file_name)
    # We reorder because otherwise the numbering might be screwed up when we combine the ModelParts later
    KratosMultiphysics.ReorderConsecutiveModelPartIO(mdpa_file_name).ReadModelPart(model_part)
    __RemoveTimeFiles()
    return model_part


def GetDefaultModelPart():
    model_part = KratosMultiphysics.ModelPart("Empty")
    prop_0 = model_part.GetProperties(0,0) # create dummy properties with id 0
    return model_part


def TranslateModelPart(model_part, translation_vector):
    '''
    Translate the ModelPart by the values in the translation_vector
    example:
    modelpart_1.TranslateModelPart([1,2,4])
    translation_vector has the format: [x, y, z]
    '''
    if (type(model_part) != KratosMultiphysics.ModelPart):
        raise Exception("input is expected to be provided as a Kratos ModelPart object")

    if len(translation_vector) != 3:
        raise Exception("Wrong length of input")

    translation_x = translation_vector[0]
    translation_y = translation_vector[1]
    translation_z = translation_vector[2]

    for node in model_part.Nodes:
        node.X0 += translation_x
        node.Y0 += translation_y
        node.Z0 += translation_z

        node.X += translation_x
        node.Y += translation_y
        node.Z += translation_z


def RotateModelPart(model_part, rotation_axis, rotation_angle, elemental_data_to_rotate=[]):
    '''
    Rotate the ModelPart around rotation_axis by rotation_angle
    example:
    modelpart_1.RotateModelPart([1,2,4], 45)
    rotation_axis has the format: [x, y, z]
    rotation_angle is in Degree
    '''
    if (type(model_part) != KratosMultiphysics.ModelPart):
            raise Exception("input is expected to be provided as a Kratos ModelPart object")

    for node in model_part.Nodes:
        RotX0Y0Z0=__RotateVector([node.X0,node.Y0,node.Z0],rotation_axis, rotation_angle)
        RotXYZ=__RotateVector([node.X,node.Y,node.Z],rotation_axis, rotation_angle)

        node.X0 = RotX0Y0Z0[0]
        node.Y0 = RotX0Y0Z0[1]
        node.Z0 = RotX0Y0Z0[2]
        node.X = RotXYZ[0]
        node.Y = RotXYZ[1]
        node.Z = RotXYZ[2]

    for elem_data_name in elemental_data_to_rotate:
        for elem in model_part.Elements:
            kratos_variable = KratosMultiphysics.KratosGlobals.GetVariable(elem_data_name)
            elem_data = elem.GetValue(kratos_variable)
            rotated_vec = __RotateVector(elem_data, rotation_axis, rotation_angle)
            elem.SetValue(kratos_variable, rotated_vec)


def AddModelPart(model_part_1, model_part_2):
    '''
    Adding the model_part_2 to model_part_1 (appending)
    '''
    if (type(model_part_1) != KratosMultiphysics.ModelPart):
            raise Exception("input is expected to be provided as a Kratos ModelPart object")
    if (type(model_part_2) != KratosMultiphysics.ModelPart):
            raise Exception("input is expected to be provided as a Kratos ModelPart object")

    num_nodes_self = model_part_1.NumberOfNodes()
    num_elements_self = model_part_1.NumberOfElements()
    num_conditions_self =model_part_1.NumberOfConditions()

    for node in model_part_2.Nodes:
        node.Id += num_nodes_self
    for element in model_part_2.Elements:
        element.Id += num_elements_self
    for condition in model_part_2.Conditions:
        condition.Id += num_conditions_self

    KratosMultiphysics.FastTransferBetweenModelPartsProcess(model_part_1, model_part_2, "All").Execute()

    # adding submodel parts of model_part_2 to model_part_1 (called recursively)
    __AddSubModelPart(model_part_1, model_part_2)


def WriteMdpaFile(model_part, new_mdpa_file_name, variables_to_write_to_gid=[]):
    if (type(model_part) != KratosMultiphysics.ModelPart):
            raise Exception("input is expected to be provided as a Kratos ModelPart object")

    if new_mdpa_file_name.endswith('.mdpa'):
        new_mdpa_file_name = new_mdpa_file_name[:-5]

    file = open(new_mdpa_file_name + ".mdpa","w")
    __WriteModelPartInfo(model_part, file)
    file.close()

    # using append bcs some info was written beforehand
    KratosMultiphysics.ModelPartIO(new_mdpa_file_name, KratosMultiphysics.IO.APPEND).WriteModelPart(model_part)
    print("#####\nWrote", model_part.Name, "to MDPA\n#####")

    ### Write the file for Visualizing in GiD
    model_part = KratosMultiphysics.ModelPart("MDPAToGID")
    model_part_io = KratosMultiphysics.ModelPartIO(new_mdpa_file_name).ReadModelPart(model_part)

    gid_mode = KratosMultiphysics.GiDPostMode.GiD_PostBinary
    multifile = KratosMultiphysics.MultiFileFlag.MultipleFiles
    deformed_mesh_flag = KratosMultiphysics.WriteDeformedMeshFlag.WriteUndeformed
    write_conditions = KratosMultiphysics.WriteConditionsFlag.WriteConditions

    gid_io = KratosMultiphysics.GidIO(new_mdpa_file_name, gid_mode, multifile,
                                deformed_mesh_flag, write_conditions)

    gid_io.InitializeMesh(0)
    gid_io.WriteMesh(model_part.GetMesh())
    gid_io.FinalizeMesh()

    gid_io.InitializeResults(0, model_part.GetMesh())

    kratos_variables_to_write = []
    for variable_name in variables_to_write_to_gid:
        kratos_variables_to_write.append(KratosMultiphysics.KratosGlobals.GetVariable(variable_name))

    for elem in model_part.Elements:
        for variable in kratos_variables_to_write:
            elem.GetNodes()[0].SetValue(variable, elem.GetValue(variable))

    for variable in kratos_variables_to_write:
        gid_io.WriteNodalResultsNonHistorical(variable, model_part.Nodes, 0)

    gid_io.FinalizeResults()

    print("#####\nWrote", model_part.Name, "to GID\n#####")

    __RemoveTimeFiles()


def __WriteModelPartInfo(model_part, open_file):
    localtime = time.asctime( time.localtime(time.time()) )
    open_file.write("// File created on " + localtime + "\n")
    open_file.write("// Mesh Information:\n")
    open_file.write("// Number of Nodes: " + str(model_part.NumberOfNodes()) + "\n")
    open_file.write("// Number of Elements: " + str(model_part.NumberOfElements()) + "\n")
    open_file.write("// Number of Conditions: " + str(model_part.NumberOfConditions()) + "\n")
    open_file.write("// Number of SubModelParts: " + str(model_part.NumberOfSubModelParts()) + "\n")
    for smp in model_part.SubModelParts:
        open_file.write("// SubModelPart " + smp.Name + "\n")
        open_file.write("//   Number of Nodes: " + str(smp.NumberOfNodes()) + "\n")
        open_file.write("//   Number of Elements: " + str(smp.NumberOfElements()) + "\n")
        open_file.write("//   Number of Conditions: " + str(smp.NumberOfConditions()) + "\n")
    open_file.write("\n")


def __AddSubModelPart(original_model_part, other_model_part):
    for smp_other in other_model_part.SubModelParts:
        if original_model_part.HasSubModelPart(smp_other.Name):
            smp_original = original_model_part.GetSubModelPart(smp_other.Name)
        else:
            smp_original = original_model_part.CreateSubModelPart(smp_other.Name)

        # making list containing node IDs of particular submodel part
        num_nodes_other = smp_other.NumberOfNodes()
        smp_node_id_array = np.zeros(num_nodes_other, dtype=np.int)
        for node, node_id in zip(smp_other.Nodes, range(num_nodes_other)):
            smp_node_id_array[node_id] = node.Id

        # making list containing element IDs of particular submodel part
        num_elements_other = smp_other.NumberOfElements()
        smp_element_id_array = np.zeros(num_elements_other, dtype=np.int)
        for element, element_id in zip(smp_other.Elements, range(num_elements_other)):
            smp_element_id_array[element_id] = element.Id

        # making list containing condition IDs of particular submodel part
        num_conditions_other = smp_other.NumberOfConditions()
        smp_condition_id_array = np.zeros(num_conditions_other, dtype=np.int)
        for condition, condition_id in zip(smp_other.Conditions, range(num_conditions_other)):
            smp_condition_id_array[condition_id] = condition.Id

        smp_original.AddNodes(smp_node_id_array.tolist())
        smp_original.AddElements(smp_element_id_array.tolist())
        smp_original.AddConditions(smp_condition_id_array.tolist())

        __AddSubModelPart(smp_original, smp_other) # call recursively to transfer nested SubModelParts


def __RotateVector(vec_to_rotate, rotation_axis, rotation_angle):
    '''
    Rotate the ModelPart around rotation_axis by rotation_angle
    example:
    modelpart_1.RotateModelPart([1,2,4], 45)
    rotation_axis has the format: [x, y, z]
    rotation_angle is in Degree
    The concept of Quaternion is used for the implementation
    Reference : https://en.wikipedia.org/wiki/Quaternions_and_spatial_rotation
    '''
    rotation_angle = rotation_angle*(math.pi/180)

    if len(rotation_axis) != 3:
        raise Exception("Wrong length of input")

    length_of_axis = math.sqrt(rotation_axis[0]**2+rotation_axis[1]**2+rotation_axis[2]**2)

    Qtnion = [0,0,0,0]
    Qtnion[0] = math.cos(rotation_angle/2)
    Qtnion[1] = (math.sin(rotation_angle/2)*rotation_axis[0])/length_of_axis
    Qtnion[2] = (math.sin(rotation_angle/2)*rotation_axis[1])/length_of_axis
    Qtnion[3] = (math.sin(rotation_angle/2)*rotation_axis[2])/length_of_axis

    x = vec_to_rotate[0]
    y = vec_to_rotate[1]
    z = vec_to_rotate[2]

    return_vector = KratosMultiphysics.Vector(3)

    return_vector[0] = 2*(Qtnion[1]*Qtnion[2]*y+Qtnion[2]*Qtnion[0]*z+Qtnion[1]*Qtnion[3]*z)+(Qtnion[1]*Qtnion[1]*x)+(Qtnion[0]*Qtnion[0]*x)-2*(Qtnion[0]*Qtnion[3]*y)-(Qtnion[3]*Qtnion[3]*x)-(Qtnion[2]*Qtnion[2]*x)

    return_vector[1] = 2*(Qtnion[1]*Qtnion[2]*x+Qtnion[0]*Qtnion[3]*x+Qtnion[2]*Qtnion[3]*z)+(Qtnion[2]*Qtnion[2]*y)+(Qtnion[0]*Qtnion[0]*y)-2*(Qtnion[0]*Qtnion[1]*z)-(Qtnion[1]*Qtnion[1]*y)-(Qtnion[3]*Qtnion[3]*y)

    return_vector[2] = 2*(Qtnion[3]*Qtnion[1]*x+Qtnion[3]*Qtnion[2]*y+Qtnion[0]*Qtnion[1]*y)+(Qtnion[0]*Qtnion[0]*z)+(Qtnion[3]*Qtnion[3]*z)-2*(Qtnion[2]*Qtnion[0]*x)-(Qtnion[2]*Qtnion[2]*z)-(Qtnion[1]*Qtnion[1]*z)

    return return_vector


def __RemoveTimeFiles():
    current_path = os.getcwd()
    files = os.listdir(current_path)
    for file in files:
        if file.endswith(".time") or file.endswith(".lst"):
            os.remove(os.path.join(current_path, file))