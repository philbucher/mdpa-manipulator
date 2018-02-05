'''
The file can take an mdpa file and translate and rotate the node and element positions.
This can also combine 2 or more mdpa files into a single one without loosing the node, element, condition details.
Transaltion is to be given in the format (x,y,z)
Rotation is to be given ([x,y,z],angle) where [x,y,z] represent the roation vector and 'angle' is the angle of rotation.
Chair of Structural Analysis, Technical University of Munich
'''

import KratosMultiphysics
import KratosMultiphysics.StructuralMechanicsApplication
import KratosMultiphysics.FluidDynamicsApplication
import numpy as np
import math
from math import cos,sin,pi
from collections import OrderedDict
import os

'''
TODO:
- Copy Variable Data! => Check if it works internally when mdpa is written (nope. not working, we have to do it manually)
- Apply Rotation to vectorial elemental data!  done 
- Copy Properties  
- Copy Tables   
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

        # Copy NodalData, ElementalData ConditionalData
        # Datastructure: {"variable_name" : { id : [value] } }
        self.model_part.nodal_data = OrderedDict()
        self.model_part.elemental_data = OrderedDict()
        self.model_part.conditional_data = OrderedDict()
        
        with open(MdpaFileName + ".mdpa", "r") as mdpa_file:
            self.ReadNodalData(mdpa_file)
            self.ReadElementalData(mdpa_file)
            self.ReadConditionalData(mdpa_file)


    def __str__(self):
        print(self.model_part)
        return ""
    
    __repr__ = __str__


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

        for node in self.model_part.Nodes:
            
            RotX0Y0Z0=RotateVector([node.X0,node.Y0,node.Z0],RotationAxis, RotationAngle)
            RotXYZ=RotateVector([node.X,node.Y,node.Z],RotationAxis, RotationAngle)

            node.X0 = RotX0Y0Z0[0]
            node.Y0 = RotX0Y0Z0[1]
            node.Z0 = RotX0Y0Z0[2]
            node.X = RotXYZ[0]
            node.Y = RotXYZ[1]
            node.Z = RotXYZ[2]
            
        # TODO check if vectorial elemental data is available and rotate as well
        if "LOCAL_AXIS_2" in self.model_part.elemental_data.keys():
            for keys,values in self.model_part.elemental_data["LOCAL_AXIS_2"].items():
                self.model_part.elemental_data["LOCAL_AXIS_2"][keys]=RotateVector([values[0],values[1],values[2]],RotationAxis, RotationAngle)
        
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
        
        #TODO add Data
        #add elemental data
        for Skeys,Svalues in self.model_part.elemental_data.items():
            for keys,values in OtherModelPart.elemental_data.items():
                if Skeys==keys:
                    temp_count = len(list(filter(bool, Svalues)))
                    for value in values:
                        self.model_part.elemental_data[Skeys][temp_count+value] = OtherModelPart.elemental_data[keys][value]
    def GetModelPart(self):
        return self.model_part


    def WriteMdpaFile(self, NewMdpaFileName):
        if NewMdpaFileName.endswith('.mdpa'):
            NewMdpaFileName = NewMdpaFileName[:-5]
        file = open(NewMdpaFileName + ".mdpa","w")
        file.close()
       
        KratosMultiphysics.ModelPartIO(NewMdpaFileName, KratosMultiphysics.IO.WRITE).WriteModelPart(self.model_part)
        print("#####\nWrote", self.model_part.Name, "to MDPA\n#####")
        
        myfile = open("sample"+ ".mdpa","a")
        
        for key in self.model_part.elemental_data:
            myfile.write('Begin ElementalData '+key+'\n')
            for key2,values in self.model_part.elemental_data[key].items():
                if key=="LOCAL_AXIS_2":
                    myfile.write('  '+str(key2)+' ('+str(values[0])+' , '+str(values[1])+' , '+str(values[2])+')\n')
                else:
                    myfile.write('  '+str(key2)+' '+str(self.model_part.elemental_data[key][key2])+'\n')
                    
            myfile.write('End ElementalData '+'\n'+'\n')
        
        myfile.close()
        
        print("##### wrote elemental datas to mdpa file #####")
        
        # TODO: Append ElementalData

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

        print("#####\nWrote", self.model_part.Name, "to GiD-Msh\n#####")

        self.RemoveTimeFiles()

    def ReadNodalData(self, MdpaFile):
        MdpaFile.seek(0)
        for line in MdpaFile:
            words = line.split()
            if len(words) > 2:
                if words[0] == "Begin" and words[1] == "NodalData":
                    print("Found NodalData for", words[2])
        
    def ReadElementalData(self, MdpaFile):
        MdpaFile.seek(0)
        
        flag='NONE'        
        for line in MdpaFile:
            words = line.split()
            if flag=="AUX_INDEX" and len(words)==2:
                    if words[0]!="End":
                        self.model_part.elemental_data["AUX_INDEX"][float(words[0])]=float(words[1])
            if flag=="LOCAL_AXIS_2" and len(words)>7:
                    self.model_part.elemental_data["LOCAL_AXIS_2"][float(words[0])]=(float(words[3]),float(words[5]),float(words[7]))

            if len(words) > 2:         
                if words[0] == "Begin" and words[1] == "ElementalData":
                    print("Found ElementalData for", words[2])
                    
                    # TODO read and append ElementalData

                    self.model_part.elemental_data[words[2]]=OrderedDict()
                    flag=words[2]
                
                if words[0] == "End" and words[1] == "ElementalData":
                    print("Ends ElementalData for", words[2])                
                    flag = 'NONE'        

    def ReadConditionalData(self, MdpaFile):
        MdpaFile.seek(0)
        for line in MdpaFile:
            words = line.split()
            if len(words) > 2:
                if words[0] == "Begin" and words[1] == "ConditionalData":
                    print("Found ConditionalData for", words[2])


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
        
        

def RotateVector(Vector, RotationAxis, RotationAngle):
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
        
        length_of_axis = math.sqrt(RotationAxis[0]**2+RotationAxis[1]**2+RotationAxis[2]**2)
                
        Qtnion = [0,0,0,0]
        Qtnion[0] = cos(RotationAngle/2)
        Qtnion[1] = (sin(RotationAngle/2)*RotationAxis[0])/length_of_axis
        Qtnion[2] = (sin(RotationAngle/2)*RotationAxis[1])/length_of_axis
        Qtnion[3] = (sin(RotationAngle/2)*RotationAxis[2])/length_of_axis
       
        x = Vector[0]
        y = Vector[1]
        z = Vector[2]          
        
        Vector[0] = 2*(Qtnion[1]*Qtnion[2]*y+Qtnion[2]*Qtnion[0]*z+Qtnion[1]*Qtnion[3]*z)+(Qtnion[1]*Qtnion[1]*x)+(Qtnion[0]*Qtnion[0]*x)-2*(Qtnion[0]*Qtnion[3]*y)-(Qtnion[3]*Qtnion[3]*x)-(Qtnion[2]*Qtnion[2]*x)
            
        Vector[1] = 2*(Qtnion[1]*Qtnion[2]*x+Qtnion[0]*Qtnion[3]*x+Qtnion[2]*Qtnion[3]*z)+(Qtnion[2]*Qtnion[2]*y)+(Qtnion[0]*Qtnion[0]*y)-2*(Qtnion[0]*Qtnion[1]*z)-(Qtnion[1]*Qtnion[1]*y)-(Qtnion[3]*Qtnion[3]*y)
            
        Vector[2] = 2*(Qtnion[3]*Qtnion[1]*x+Qtnion[3]*Qtnion[2]*y+Qtnion[0]*Qtnion[1]*y)+(Qtnion[0]*Qtnion[0]*z)+(Qtnion[3]*Qtnion[3]*z)-2*(Qtnion[2]*Qtnion[0]*x)-(Qtnion[2]*Qtnion[2]*z)-(Qtnion[1]*Qtnion[1]*z)
        
        return Vector
        
        