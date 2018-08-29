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


def ReadModelPart(mdpa_file_name, model_part_name):
    '''
    Read and return a ModelPart from a mdpa file
    '''
    if mdpa_file_name.endswith('.mdpa'):
        mdpa_file_name = mdpa_file_name[:-5]
    model_part = KratosMultiphysics.ModelPart(model_part_name)
    # We reorder because otherwise the numbering might be screwed up when we combine the ModelParts later
    KratosMultiphysics.ReorderConsecutiveModelPartIO(mdpa_file_name).ReadModelPart(model_part)
    __RemoveAuxFiles()
    return model_part


def GetDefaultModelPart(model_part_name):
    '''
    Create and return an empty "dummy" ModelPart to which the other ModelParts are to be added
    in order to have a clean "start"
    '''
    model_part = KratosMultiphysics.ModelPart(model_part_name)
    prop_0 = model_part.GetProperties(0,0) # create dummy properties with id
    return model_part


def TranslateModelPart(model_part,
                       translation_vector):
    '''
    Translate the ModelPart by the values in the translation_vector
    example:
    TranslateModelPart(model_part, [1,2,4])
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


def RotateModelPart(model_part,
                    rotation_axis,
                    rotation_angle,
                    elemental_data_to_rotate=[],
                    conditional_data_to_rotate=[]):
    '''
    Rotate the ModelPart around rotation_axis by rotation_angle
    It can also rotate vectorial elemental/conditional data
    example:
    RotateModelPart(model_part, [1,2,4], 45)
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
        kratos_variable = KratosMultiphysics.KratosGlobals.GetVariable(elem_data_name)
        for elem in model_part.Elements:
            if elem.Has(kratos_variable):
                elem_data = elem.GetValue(kratos_variable)
                rotated_vec = __RotateVector(elem_data, rotation_axis, rotation_angle)
                elem.SetValue(kratos_variable, rotated_vec)

    for cond_data_name in conditional_data_to_rotate:
        kratos_variable = KratosMultiphysics.KratosGlobals.GetVariable(cond_data_name)
        for cond in model_part.Conditions:
            if cond.Has(kratos_variable):
                cond_data = cond.GetValue(kratos_variable)
                rotated_vec = __RotateVector(cond_data, rotation_axis, rotation_angle)
                cond.SetValue(kratos_variable, rotated_vec)


def AddModelPart(model_part_1,
                 model_part_2,
                 add_as_submodelpart=False):
    '''
    Adding the model_part_2 to model_part_1 (appending)
    '''
    if (type(model_part_1) != KratosMultiphysics.ModelPart):
            raise Exception("input is expected to be provided as a Kratos ModelPart object")
    if (type(model_part_2) != KratosMultiphysics.ModelPart):
            raise Exception("input is expected to be provided as a Kratos ModelPart object")

    num_nodes_self = model_part_1.NumberOfNodes()
    num_elements_self = model_part_1.NumberOfElements()
    num_conditions_self = model_part_1.NumberOfConditions()

    for node in model_part_2.Nodes:
        node.Id += num_nodes_self
    for element in model_part_2.Elements:
        element.Id += num_elements_self
    for condition in model_part_2.Conditions:
        condition.Id += num_conditions_self

    KratosMultiphysics.FastTransferBetweenModelPartsProcess(model_part_1, model_part_2, KratosMultiphysics.FastTransferBetweenModelPartsProcess.EntityTransfered.ALL).Execute()

    if add_as_submodelpart:
        # adding model_part_2 as submodel part to model_part_1 (called recursively)
        __AddAsSubModelPart(model_part_1, model_part_2)
    else:
        # adding submodel parts of model_part_2 to model_part_1 (called recursively)
        __AddSubModelPart(model_part_1, model_part_2)


def WriteMdpaFile(model_part,
                  mdpa_file_name,
                  variables_to_write_to_gid=[],
                  assing_properties=False):
    '''
    Writes the mdpa file from a ModelPart with some additional information
    Also a GiD file is created for visual postprocessing
    '''
    if (type(model_part) != KratosMultiphysics.ModelPart):
            raise Exception("input is expected to be provided as a Kratos ModelPart object")

    if mdpa_file_name.endswith('.mdpa'):
        mdpa_file_name = mdpa_file_name[:-5]

    file = open(mdpa_file_name + ".mdpa","w")
    __WriteModelPartInfo(model_part, file)
    file.close()

    # using append bcs some info was written beforehand
    KratosMultiphysics.ModelPartIO(mdpa_file_name, KratosMultiphysics.IO.APPEND).WriteModelPart(model_part)
    print("#####\nWrote ModelPart to MDPA\n#####")

    ### Write the file for Visualizing in GiD
    gid_model_part = KratosMultiphysics.ModelPart("MDPAToGID")
    KratosMultiphysics.ModelPartIO(mdpa_file_name).ReadModelPart(gid_model_part)

    if assing_properties:
        # assign different Properties to the elems/conds to visualize the smps in GiD
        property_counter=1
        for smp in gid_model_part.SubModelParts:
            prop = model_part.GetProperties(property_counter,0) # mesh_id=0
            for elem in smp.Elements:
                elem.Properties = prop
            for cond in smp.Conditions:
                cond.Properties = prop
            property_counter += 1

    gid_mode = KratosMultiphysics.GiDPostMode.GiD_PostBinary
    multifile = KratosMultiphysics.MultiFileFlag.MultipleFiles
    deformed_mesh_flag = KratosMultiphysics.WriteDeformedMeshFlag.WriteUndeformed
    write_conditions = KratosMultiphysics.WriteConditionsFlag.WriteConditions

    gid_io = KratosMultiphysics.GidIO(mdpa_file_name, gid_mode, multifile,
                                deformed_mesh_flag, write_conditions)

    gid_io.InitializeMesh(0)
    gid_io.WriteMesh(gid_model_part.GetMesh())
    gid_io.FinalizeMesh()

    gid_io.InitializeResults(0, gid_model_part.GetMesh())

    kratos_variables_to_write = []
    for variable_name in variables_to_write_to_gid:
        kratos_variables_to_write.append(KratosMultiphysics.KratosGlobals.GetVariable(variable_name))

    for elem in gid_model_part.Elements:
        for variable in kratos_variables_to_write:
            elem.GetNodes()[0].SetValue(variable, elem.GetValue(variable))

    for variable in kratos_variables_to_write:
        gid_io.WriteNodalResultsNonHistorical(variable, gid_model_part.Nodes, 0)

    gid_io.FinalizeResults()

    print("#####\nWrote ModelPart to GID\n#####")

    __RemoveAuxFiles()


def __WriteModelPartInfo(model_part,
                         open_file):
    '''
    Writing some information about the ModelPart to the mdpa file
    '''
    localtime = time.asctime( time.localtime(time.time()) )
    open_file.write("// File created on " + localtime + "\n")
    open_file.write("// Mesh Information:\n")
    open_file.write("// ModelPartName: " + model_part.Name + "\n")
    open_file.write("// Number of Nodes: " + str(model_part.NumberOfNodes()) + "\n")
    open_file.write("// Number of Elements: " + str(model_part.NumberOfElements()) + "\n")
    open_file.write("// Number of Conditions: " + str(model_part.NumberOfConditions()) + "\n")
    open_file.write("// Number of SubModelParts: " + str(model_part.NumberOfSubModelParts()) + "\n")
    __WriteSubModelPartInfo(model_part,open_file, level=0)
    open_file.write("\n")

def __WriteSubModelPartInfo(model_part,
                            open_file,
                            level):
    SPACE = "    "
    for smp in model_part.SubModelParts:
        full_name = [smp.Name]
        psmp = smp
        while psmp.IsSubModelPart():
            psmp = psmp.GetParentModelPart()
            full_name.append(psmp.Name)
        open_file.write("// " + SPACE*level + "SubModelPart " + ".".join(full_name[::-1]) + "\n")
        open_file.write("// " + SPACE*level + "Number of Nodes: " + str(smp.NumberOfNodes()) + "\n")
        open_file.write("// " + SPACE*level + "Number of Elements: " + str(smp.NumberOfElements()) + "\n")
        open_file.write("// " + SPACE*level + "Number of Conditions: " + str(smp.NumberOfConditions()) + "\n")
        open_file.write("// " + SPACE*level + "Number of SubModelParts: " + str(smp.NumberOfSubModelParts()) + "\n")
        __WriteSubModelPartInfo(smp,open_file, level+1)

def __AddEntitiesToSubModelPart(original_sub_model_part,
                                other_sub_model_part):
    '''
    Adds the entities of (nodes, elements and conditions) from
    one SubModelPart to another
    '''
    # making list containing node IDs of particular submodel part
    num_nodes_other = other_sub_model_part.NumberOfNodes()
    smp_node_id_array = np.zeros(num_nodes_other, dtype=np.int)
    for node_i, node in enumerate(other_sub_model_part.Nodes):
        smp_node_id_array[node_i] = node.Id

    # making list containing element IDs of particular submodel part
    num_elements_other = other_sub_model_part.NumberOfElements()
    smp_element_id_array = np.zeros(num_elements_other, dtype=np.int)
    for element_i, element in enumerate(other_sub_model_part.Elements):
        smp_element_id_array[element_i] = element.Id

    # making list containing condition IDs of particular submodel part
    num_conditions_other = other_sub_model_part.NumberOfConditions()
    smp_condition_id_array = np.zeros(num_conditions_other, dtype=np.int)
    for condition_i, condition in enumerate(other_sub_model_part.Conditions):
        smp_condition_id_array[condition_i] = condition.Id

    original_sub_model_part.AddNodes(smp_node_id_array.tolist())
    original_sub_model_part.AddElements(smp_element_id_array.tolist())
    original_sub_model_part.AddConditions(smp_condition_id_array.tolist())

def __AddSubModelPart(original_model_part,
                      other_model_part):
    '''
    Adds the SubModelParts of one ModelPart to the other one
    If the original ModelPart already contains a SMP with the same name,
    the entities are added to it
    '''
    for smp_other in other_model_part.SubModelParts:
        if original_model_part.HasSubModelPart(smp_other.Name):
            smp_original = original_model_part.GetSubModelPart(smp_other.Name)
        else:
            smp_original = original_model_part.CreateSubModelPart(smp_other.Name)

        __AddEntitiesToSubModelPart(smp_original, smp_other)

        __AddSubModelPart(smp_original, smp_other) # call recursively to transfer nested SubModelParts

def __AddAsSubModelPart(original_model_part,
                        other_model_part):
    '''
    Adds the SubModelParts of one ModelPart to the other one
    If the original ModelPart already contains a SMP with the same name,
    the entities are added to it
    '''
    smp_original = original_model_part.CreateSubModelPart(other_model_part.Name)

    __AddEntitiesToSubModelPart(smp_original, other_model_part)

    for smp_other in other_model_part.SubModelParts:
        __AddAsSubModelPart(smp_original, smp_other) # call recursively to transfer nested SubModelParts

def __RotateVector(vec_to_rotate,
                   rotation_axis,
                   rotation_angle):
    '''
    Rotates a generic vector around rotation_axis by rotation_angle
    example:
    __RotateVector([0,12,-3], [1,2,4], 45)
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

def __RemoveAuxFiles():
    '''
    Removes auxiliary files from the directory
    '''
    current_path = os.getcwd()
    files = os.listdir(current_path)
    for file in files:
        if file.endswith(".time") or file.endswith(".lst"):
            os.remove(os.path.join(current_path, file))

def DictToPrettyString(o, level=0):
    '''
    This function returns a nicely formatted string of a dictionary
    This can e.g. be written to a file
    '''
    # source: https://stackoverflow.com/questions/10097477/python-json-array-newlines
    INDENT = 4
    SPACE = " "
    NEWLINE = "\n"
    ret = ""
    if isinstance(o, dict):
        if len(o) == 0:
            ret += "{}"
        elif len(o) == 1: # and list(o.keys())[0] == "@table"
            ret += "{"
            for k in o.keys():
                v = o[k]
                ret += '"' + str(k) + '" : '
                ret += WriteToJSON(v, level + 1)

            ret += "}"
        else:
            ret += "{" + NEWLINE
            comma = ""
            for k in reversed(sorted(o.keys())): # reversed to make i better readable
                v = o[k]
                ret += comma
                comma = ",\n"
                ret += SPACE * INDENT * (level+1)
                ret += '"' + str(k) + '" :' + SPACE
                ret += WriteToJSON(v, level + 1)

            ret += NEWLINE + SPACE * INDENT * level + "}"
    elif isinstance(o, str):
        ret += '"' + o + '"'
    elif isinstance(o, list):
        if len(o) > 0:
            if isinstance(o[0], (int,  float, str)):
                ret += NEWLINE + SPACE * INDENT * (level+1) + "["
                num_vals = len(o)
                for i in range(num_vals):
                    entry = o[i]
                    ret += WriteValue(entry)
                    if i <= num_vals-2:
                        ret += ", "
                    else:
                        ret += "]"
            else:
                ret += "[" + ",".join([WriteToJSON(e, level+1) for e in o]) + "]"
        else:
            ret += "[]"
    elif isinstance(o, bool):
        ret += "true" if o else "false"
    elif isinstance(o, int):
        ret += str(o)
    elif isinstance(o, float):
        ret += '%.7g' % o
    elif isinstance(o, np.ndarray) and np.issubdtype(o.dtype, np.integer):
        ret += "[" + ','.join(map(str, o.flatten().tolist())) + "]"
    elif isinstance(o, np.ndarray) and np.issubdtype(o.dtype, np.inexact):
        ret += "[" + ','.join(map(lambda x: '%.7g' % x, o.flatten().tolist())) + "]"
    elif o is None:
        ret += 'null'
    elif isinstance(o, TabularData):
        data = o.GetTabularData()
        ret += "[" + NEWLINE
        num_rows = len(data)
        for i in range(num_rows-1):
            row = data[i]
            ret += SPACE * INDENT * (level+1) + "[" + str(row[0]) + " , " + str(row[1]) + "]," + NEWLINE
        ret += SPACE * INDENT * (level+1) + "[" + str(row[0]) + " , " + str(row[1]) + "]" + NEWLINE
        ret += SPACE * INDENT * (level+1) + "]"
    else:
        raise TypeError("Unknown type '%s' for json serialization" % str(type(o)))
    return ret