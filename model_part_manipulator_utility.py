import KratosMultiphysics
import KratosMultiphysics.StructuralMechanicsApplication
import KratosMultiphysics.FluidDynamicsApplication
import numpy as np
import math
from math import cos,sin,pi
import os

'''
TODO:
- Copy Variable Data! => Check if it works internally when mdpa is written
- Apply Rotation to vectorial elemental data!
- Copy tables
- Copy Properties
(- Copy ModelpartData?)
'''


class ModelPartManipulator:
    def __init__(self, MdpaFileName):
        
        #To remove .mpda if the file name is given as .mdpa
        if MdpaFileName.endswith('.mdpa'):
            MdpaFileName = MdpaFileName[:-5]
        self.model_part = KratosMultiphysics.ModelPart(MdpaFileName)
        # We reorder because otherwise the numbering might be screwed up when we combine the ModelParts later
        KratosMultiphysics.ReorderConsecutiveModelPartIO(MdpaFileName).ReadModelPart(self.model_part)

        self.RemoveTimeFiles()


    def TranslateModelPart(self, TranslationVector):  
        '''
        Translate the ModelPart by the values in the TranslationVector
        example:
        modelpart_1.TranslateModelPart([1,2,4])
        TranslationVector has the format: [x, y, z]
        '''
        if len(TranslationVector) != 3:
            raise Exception("Wrong length of input")

        translation_x = TranslationVector[0]
        translation_y = TranslationVector[1]
        translation_z = TranslationVector[2]

        for node in self.model_part.Nodes:
            node.X0 += translation_x
            node.Y0 += translation_y
            node.Z0 += translation_z

            node.X += translation_x
            node.Y += translation_y
            node.Z += translation_z


    def RotateModelPart(self, RotationAxis, RotationAngle):
        '''
        Rotate the ModelPart around RotationAxis by RotationAngle
        example:
        modelpart_1.RotateModelPart([1,2,4], 45)
        RotationAxis has the format: [x, y, z]
        RotationAngle is in Degree
        The concept of Quaternion is used for the implementation
        Reference : https://en.wikipedia.org/wiki/Quaternions_and_spatial_rotation
        '''
        RotationAngle = RotationAngle*(pi/180)
        
        if len(RotationAxis) != 3:
            raise Exception("Wrong length of input")
        
        Length_of_axis = math.sqrt(RotationAxis[0]**2+RotationAxis[1]**2+RotationAxis[2]**2)
                
        Qtnion = [0,0,0,0]
        Qtnion[0] = cos(RotationAngle/2)
        Qtnion[1] = (sin(RotationAngle/2)*RotationAxis[0])/Length_of_axis
        Qtnion[2] = (sin(RotationAngle/2)*RotationAxis[1])/Length_of_axis
        Qtnion[3] = (sin(RotationAngle/2)*RotationAxis[2])/Length_of_axis
       
        for node in self.model_part.Nodes:
            x0 = node.X0
            y0 = node.Y0
            z0 = node.Z0
            x = node.X
            y = node.Y
            z = node.Z
            
            node.X0 = 2*(Qtnion[1]*Qtnion[2]*y0+Qtnion[2]*Qtnion[0]*z0+Qtnion[1]*Qtnion[3]*z0)+(Qtnion[1]*Qtnion[1]*x0)+(Qtnion[0]*Qtnion[0]*x0)-2*(Qtnion[0]*Qtnion[3]*y0)-(Qtnion[3]*Qtnion[3]*x0)-(Qtnion[2]*Qtnion[2]*x0)
            
            node.Y0 = 2*(Qtnion[1]*Qtnion[2]*x0+Qtnion[0]*Qtnion[3]*x0+Qtnion[2]*Qtnion[3]*z0)+(Qtnion[2]*Qtnion[2]*y0)+(Qtnion[0]*Qtnion[0]*y0)-2*(Qtnion[0]*Qtnion[1]*z0)-(Qtnion[1]*Qtnion[1]*y0)-(Qtnion[3]*Qtnion[3]*y0)
            
            node.Z0 = 2*(Qtnion[3]*Qtnion[1]*x0+Qtnion[3]*Qtnion[2]*y0+Qtnion[0]*Qtnion[1]*y0)+(Qtnion[0]*Qtnion[0]*z0)+(Qtnion[3]*Qtnion[3]*z0)-2*(Qtnion[2]*Qtnion[0]*x0)-(Qtnion[2]*Qtnion[2]*z0)-(Qtnion[1]*Qtnion[1]*z0)
            
            node.X = 2*(Qtnion[1]*Qtnion[2]*y+Qtnion[2]*Qtnion[0]*z+Qtnion[1]*Qtnion[3]*z)+(Qtnion[1]*Qtnion[1]*x)+(Qtnion[0]*Qtnion[0]*x)-2*(Qtnion[0]*Qtnion[3]*y)-(Qtnion[3]*Qtnion[3]*x)-(Qtnion[2]*Qtnion[2]*x)
            
            node.Y = 2*(Qtnion[1]*Qtnion[2]*x+Qtnion[0]*Qtnion[3]*x+Qtnion[2]*Qtnion[3]*z)+(Qtnion[2]*Qtnion[2]*y)+(Qtnion[0]*Qtnion[0]*y)-2*(Qtnion[0]*Qtnion[1]*z)-(Qtnion[1]*Qtnion[1]*y)-(Qtnion[3]*Qtnion[3]*y)
            
            node.Z = 2*(Qtnion[3]*Qtnion[1]*x+Qtnion[3]*Qtnion[2]*y+Qtnion[0]*Qtnion[1]*y)+(Qtnion[0]*Qtnion[0]*z)+(Qtnion[3]*Qtnion[3]*z)-2*(Qtnion[2]*Qtnion[0]*x)-(Qtnion[2]*Qtnion[2]*z)-(Qtnion[1]*Qtnion[1]*z)


    def AddModelPart(self, OtherModelPart):
        '''
        Adding the OtherModelPart to self.model_part (appending)
        '''
        num_nodes_self = self.model_part.NumberOfNodes()
        num_elements_self = self.model_part.NumberOfElements()
        num_conditions_self = self.model_part.NumberOfConditions()

        for node in OtherModelPart.Nodes:
            node.Id += num_nodes_self
        for element in OtherModelPart.Elements:
            element.Id += num_elements_self
        for condition in OtherModelPart.Conditions:
            condition.Id += num_conditions_self

        KratosMultiphysics.FastTransferBetweenModelPartsProcess(self.model_part, OtherModelPart, "All").Execute()
        
        #adding submodel parts of OtherModelPart to self.model_part (called recursively)
        AddSubModelPart(self, OtherModelPart)
    
                
    def GetModelPart(self):
        return self.model_part


    def WriteMdpaFile(self, NewMdpaFileName):
        if NewMdpaFileName.endswith('.mdpa'):
            NewMdpaFileName = NewMdpaFileName[:-5]
        file = open(NewMdpaFileName + ".mdpa","w")
        file.close()
       
        KratosMultiphysics.ModelPartIO(NewMdpaFileName, KratosMultiphysics.IO.WRITE).WriteModelPart(self.model_part)

        ### Write the file to msh for Visualizing in GiD
        model_part = KratosMultiphysics.ModelPart("MDPAToGID")

        model_part_io = KratosMultiphysics.ModelPartIO(NewMdpaFileName).ReadModelPart(model_part)

        gid_mode = KratosMultiphysics.GiDPostMode.GiD_PostAscii
        multifile = KratosMultiphysics.MultiFileFlag.MultipleFiles
        deformed_mesh_flag = KratosMultiphysics.WriteDeformedMeshFlag.WriteUndeformed
        write_conditions = KratosMultiphysics.WriteConditionsFlag.WriteConditions

        gid_io = KratosMultiphysics.GidIO(NewMdpaFileName, gid_mode, multifile,
                                   deformed_mesh_flag, write_conditions)

        gid_io.InitializeMesh(0)
        gid_io.WriteMesh(model_part.GetMesh())
        gid_io.FinalizeMesh()

        self.RemoveTimeFiles()


    def RemoveTimeFiles(self):
        current_path = os.getcwd()
        files = os.listdir(current_path)
        for file in files:
            if file.endswith(".time"):
                os.remove(os.path.join(current_path, file))
        

def AddSubModelPart(OriginalModelPart, OtherModelPart):
    for smp_other in OtherModelPart.SubModelParts:
        if OriginalModelPart.model_part.HasSubModelPart(smp_other.Name):
            smp_original = OriginalModelPart.model_part.GetSubModelPart(smp_other.Name)
        else:
            smp_original = OriginalModelPart.model_part.CreateSubModelPart(smp_other.Name)
        
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
    
        AddSubModelPart(smp_original, smp_other) # call recursively to transfer nested SubModelParts