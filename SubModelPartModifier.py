import KratosMultiphysics
import KratosMultiphysics.StructuralMechanicsApplication as KratosStructural
import KratosMultiphysics.FluidDynamicsApplication as KratosCFD
import numpy as np



def AddSubModelPart(OriginalModelPart,OtherModelPart):
        
        for smp in OtherModelPart.SubModelParts:
            if OriginalModelPart.model_part.HasSubModelPart(smp.Name):
                sssmp = OriginalModelPart.model_part.GetSubModelPart(smp.Name)
            else:
                sssmp = OriginalModelPart.model_part.CreateSubModelPart(smp.Name)
            
            # making list containing node IDs of particular submodel part
            NodeIdArray = np.zeros(smp.NumberOfNodes(),dtype=np.int)
            counter =0
            for node in smp.Nodes:
                NodeIdArray[counter] = node.Id 
                counter = counter+1            
            NodeIdList = np.array(NodeIdArray).tolist()
            
            # making list containing element IDs of particular submodel part
            ElementIdArray = np.zeros(smp.NumberOfElements(),dtype=np.int)
            counter =0
            for element in smp.Elements:
                ElementIdArray[counter] = element.Id 
                counter = counter+1            
            ElementIdList = np.array(ElementIdArray).tolist()
            
            # making list containing condition IDs of particular submodel part
            ConditionIdArray = np.zeros(smp.NumberOfConditions(),dtype=np.int)
            counter =0
            for condition in smp.Conditions:
                ConditionIdArray[counter] =condition.Id 
                counter = counter+1            
            ConditionIdList = np.array(ConditionIdArray).tolist()
                      
            sssmp.AddNodes(NodeIdList)
            sssmp.AddElements(ElementIdList)
            sssmp.AddConditions(ConditionIdList)
        
            AddSubModelPart(sssmp,smp)
            