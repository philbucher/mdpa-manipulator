import KratosMultiphysics
import KratosMultiphysics.StructuralMechanicsApplication as KratosStructural
import KratosMultiphysics.FluidDynamicsApplication as KratosCFD
import numpy as np
import math
from math import cos,sin,pi
from os import remove
import SubModelPartModifier as smpm


class ModelPartManipulator:
   
    def __init__(self, MdpaFileName):
        
        #To remove .mpda if the file name is given as .mdpa
        if MdpaFileName.endswith('.mdpa'):
            MdpaFileName = MdpaFileName[:-5]
        self.model_part = KratosMultiphysics.ModelPart(MdpaFileName)
        KratosMultiphysics.ModelPartIO(MdpaFileName).ReadModelPart(self.model_part)
        
        try:
            remove(MdpaFileName + ".time")
            print("*.time file removed sucessfully")
        except FileNotFoundError as e:
            print("*.time file could not be removed")

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
        Model_1_Total_Nodes = self.model_part.NumberOfNodes()
        Model_1_Total_Elements = self.model_part.NumberOfElements()
        Model_1_Total_Conditions = self.model_part.NumberOfConditions()
       
        for node in OtherModelPart.Nodes:
            node.Id += Model_1_Total_Nodes
        for element in OtherModelPart.Elements:
            element.Id += Model_1_Total_Elements
        for condition in OtherModelPart.Conditions:
            condition.Id += Model_1_Total_Conditions
            
        KratosMultiphysics.FastTransferBetweenModelPartsProcess(self.model_part, OtherModelPart, "All").Execute()
        
        #adding submodel parts of OtherModelPart to self.model_part
        smpm.AddSubModelPart(self,OtherModelPart)
    
                
    def GetModelPart(self):
        return self.model_part


    def WriteMdpaFile(self, NewMdpaFileName):
        if NewMdpaFileName.endswith('.mdpa'):
            NewMdpaFileName = NewMdpaFileName[:-5]
        file = open(NewMdpaFileName + ".mdpa","w")
        file.close()
       
        KratosMultiphysics.ModelPartIO(NewMdpaFileName, KratosMultiphysics.IO.WRITE).WriteModelPart(self.model_part)
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

        try:
            remove(NewMdpaFileName + ".time")
            print("*.time file removed sucessfully")
        except FileNotFoundError as e:
            print("*.time file could not be removed")
        