#!/usr/bin/env python2
## #!/usr/bin/env python3
# filename: reconstruct_aggregates.py
# author: Jon David
# date: Wednesday, March 21, 2018
# description:
#   Reads a vizmo path file are renders a subset of those
#   configurations in PyMOL.
#--------------------------------------------------------------------
# notes:
#   3.21.18: replaced all "  " to " " in p0.mb1n1c.path.noHeader
#--------------------------------------------------------------------


class PathData(object):
    
    def __init__(self, fname, numModels):
        self.stepDict = {}
        self.readFile(fname, numModels)


    def readFile(self, fname, numModels):
        """Reads and parses and input vizmo path file to populate
        a dictionary of type (StepID, [ModelConfiguration]),
        where ModelConfiguration is of type (ModelID,x,y,z,a,b,g).
        Note: fname is assumed to have no header lines."""
        stepID = 0
        with open(fname, 'r') as f:
            for line in f:
                modelConfigList = self.parse(line, numModels)
                #print(modelConfigList)
                #print("stepID: {}".format(stepID))
                self.stepDict[stepID] = modelConfigList
                stepID = stepID + 1
        return

    
    def parse(self, line, numModels):
        """Parses a line of the form AB, where A is of type [(x,y,z)],
        and B is of type [(a,g,b)].
        Return is of type [(ModelID,x,y,z,a,b,g)]."""
        line = line.strip()
        itemList = line.split(' ')      ## line.split(sep=' ') for python3
        #print("len(itemList): {}".format( len(itemList )))
        #print("itemList: {}".format( itemList ))
        
        positionList = []   ## has type: [(x,y,z)]
        receptorID = 0
        for i in range(0, numModels):
            xID = 3*receptorID  ## x_index, where to find this x in itemList
            yID = xID + 1
            zID = xID + 2
            #print("(xID,yID,zID) = ({}, {}, {})".format(xID, yID, zID))
            position = (itemList[xID], itemList[yID], itemList[zID])
            positionList.append(position)
            receptorID = receptorID + 1

        rotationList = []  ## has type: [(a,b,g)]
        baseOffset = numModels*3
        allergenID = 0
        for j in range(0, numModels):
            aID = baseOffset + 3*allergenID
            bID = aID + 1
            gID = aID + 2
            #print("(aID, bID, gID) = ({}, {}, {}) ".format(aID,bID,gID))
            rotation = (itemList[aID], itemList[bID], itemList[gID])
            rotationList.append(rotation)
            allergenID = allergenID + 1

        modelConfigList = []
        for k in range(0, numModels):
            modelConfig = positionList[k] + rotationList[k]  #one 6-tuple
            modelConfigList.append(modelConfig)
            
        return modelConfigList

    def getConfiguration(self, stepID, modelID):
        """Gets the configuration of model  modelID in step stepID."""
        return self.stepDict[stepID][modelID]



##---- test functions ----

def test_PathData():
    fname = "input/test_input"
    numModels = 2
    data = PathData(fname, numModels)
    print("(Step0, Model0) is: {}".format( data.getConfiguration(0,0) ))
    print("(Step0, Model1) is: {}".format( data.getConfiguration(0,1) ))
    print("(Step1, Model0) is: {}".format( data.getConfiguration(1,0) ))
    print("(Step1, Model1) is: {}".format( data.getConfiguration(1,1) ))

def test_outputAllConfigsInLastStep():
    """I used this function to display all configurations. This helps me
    find all molecules in a particular quadrant."""
    fname = "input/p0.mb1n1c.path.noHeader"
    numModels = 20
    data = PathData(fname, numModels)
    for modelID in range(0, numModels):
        print("(Step999,Model{}) is: {}".format(modelID, data.getConfiguration(999, modelID)))

    ## modelIDs in quadrant 3 (-x values, +z values):
    # (Step999,Model1) is: ('-43.3077', '-4.59132', '27.605', '0', '0.122539', '0')
    # (Step999,Model2) is: ('-56.2992', '-4.59132', '13.7883', '0', '0.392702', '0')
    # (Step999,Model4) is: ('-58.7896', '-4.59132', '27.6803', '0', '0.682062', '0')
    # (Step999,Model6) is: ('-41.8501', '-4.59132', '41.8614', '0', '0.964607', '0')
    # (Step999,Model13) is: ('-45.4026', '0.01561', '33.9709', '0', '0.372954', '0')
    # (Step999,Model14) is: ('-51.3592', '0.01561', '24.7466', '0', '0.130871', '0')
    # (Step999,Model16) is: ('-50.752', '0.01561', '16.224', '0', '0.528833', '0')
    # (Step999,Model18) is: ('-57.1006', '0.01561', '22.1472', '0', '0.775501', '0')



##---- main ----
#test_PathData()
test_outputAllConfigsInLastStep()
