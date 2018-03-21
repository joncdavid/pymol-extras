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


from pymol import cmd


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
                self.stepDict[stepID] = modelConfigList
                stepID = stepID + 1
        return

    
    def parse(self, line, numModels):
        """Parses a line of the form AB, where A is of type [(x,y,z)],
        and B is of type [(a,g,b)].
        Return is of type [(x,y,z,a,b,g)]."""
        line = line.strip()
        itemList = line.split(' ')      ## line.split(sep=' ') for python3
        
        positionList = []   ## has type: [(x,y,z)]
        receptorID = 0
        for i in range(0, numModels):
            xID = 3*receptorID  ## x_index, where to find this x in itemList
            yID = xID + 1
            zID = xID + 2
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


class ConfigurationRenderer(object):
    def __init__(self):
        ## initialize PyMOL here
        return

    def render(self, modelConfiguration, pdb_fname, modelName, color_str="green"):
        """Renders a pdb with given configuration in PyMOL.
        modelConfiguration is of type (x,y,z,a,b,g);
        model_fname is of type string and represents the model PDB file."""
        
        x = 10*float(modelConfiguration[0])  ## scale "up" by a factor of 10
        y = 10*float(modelConfiguration[1])  ## because modelConfiguration's x,y,z have units of nm
        z = 10*float(modelConfiguration[2])  ## but PyMOL uses coordinates of Angstroms
        translationVector = "[{},{},{}]".format(x,y,z)
        xDegOfRot = 360*float(modelConfiguration[3])  ## alpha, is in units of a fraction of 360 degrees
        yDegOfRot = 360*float(modelConfiguration[4])
        zDegOfRot = 360*float(modelConfiguration[5])
        
        #cmd.load(pdb_fname, "original_{}".format(modelName))  ## for debug, something to compare against
        cmd.load(pdb_fname, modelName)
        cmd.rotate('x', xDegOfRot, modelName, 0, 0)  ## rotate about x-axis
        cmd.rotate('y', yDegOfRot, modelName, 0, 0)  ## rotate about y-axis
        cmd.rotate('z', zDegOfRot, modelName, 0, 0)  ## rotate about z-axis
        ## note: apply rotations before translation, because it rotates about the origin
        cmd.translate(translationVector, modelName, 0, 0)  ## all states, and _not_ camera coordinates

        cmd.color(color_str, modelName)

    def final_render(self):
        cmd.hide()
        cmd.show('cartoon')

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

    ## my filterested modelIDs in quadrant 3 (-x values, +z values):
    # (Step999,Model1) is: ('-43.3077', '-4.59132', '27.605', '0', '0.122539', '0')
    # (Step999,Model2) is: ('-56.2992', '-4.59132', '13.7883', '0', '0.392702', '0')
    # (Step999,Model4) is: ('-58.7896', '-4.59132', '27.6803', '0', '0.682062', '0')
    # (Step999,Model6) is: ('-41.8501', '-4.59132', '41.8614', '0', '0.964607', '0')
    # (Step999,Model13) is: ('-45.4026', '0.01561', '33.9709', '0', '0.372954', '0')
    # (Step999,Model14) is: ('-51.3592', '0.01561', '24.7466', '0', '0.130871', '0')
    # (Step999,Model16) is: ('-50.752', '0.01561', '16.224', '0', '0.528833', '0')
    # (Step999,Model18) is: ('-57.1006', '0.01561', '22.1472', '0', '0.775501', '0')


def test_ConfigurationRenderer():
    fname = "input/p0.mb1n1c.path.noHeader"
    numModels = 20
    data = PathData(fname, numModels)

    ## note: models 0-9 are IgE receptors
    ##       and models 10-19 are mutant allergens
    receptor_pdb_fname = "./input_pdbs/Rec.pdb"
    allergen_pdb_fname = "./input_pdbs/mutant-MB1N1C-singleFile.yAligned.pdb"
    model_name = "model1"
    modelConfiguration = data.getConfiguration(999,1)
    configRenderer = ConfigurationRenderer()
    configRenderer.render(modelConfiguration, receptor_pdb_fname, model_name)

def test_ConfigurationRenderer_manyModels():
    fname = "input/p0.mb1n1c.path.noHeader"
    numModels = 20
    data = PathData(fname, numModels)

    ## note: models 0-9 are IgE receptors
    ##       and models 10-19 are mutant allergens
    receptor_pdb_fname = "./input_pdbs/Rec.pdb"
    allergen_pdb_fname = "./input_pdbs/mutant-MB1N1C-singleFile.yAligned.pdb"
    modelConfig_1 = data.getConfiguration(999,1)
    modelConfig_2 = data.getConfiguration(999,2)
    modelConfig_4 = data.getConfiguration(999,4)
    modelConfig_6 = data.getConfiguration(999,6)
    modelConfig_13 = data.getConfiguration(999,13)
    modelConfig_14 = data.getConfiguration(999,14)
    modelConfig_16 = data.getConfiguration(999,16)
    modelConfig_18 = data.getConfiguration(999,18)
    receptor_color_str = 'wheat'
    allergen_color_str = 'green'

    configRenderer = ConfigurationRenderer()
    configRenderer.render(modelConfig_1, receptor_pdb_fname, "model1-rec", receptor_color_str)
    configRenderer.render(modelConfig_2, receptor_pdb_fname, "model2-rec", receptor_color_str)
    configRenderer.render(modelConfig_4, receptor_pdb_fname, "model4-rec", receptor_color_str)
    configRenderer.render(modelConfig_6, receptor_pdb_fname, "model6-rec", receptor_color_str)
    configRenderer.render(modelConfig_13, allergen_pdb_fname, "model13-alg", allergen_color_str)
    configRenderer.render(modelConfig_14, allergen_pdb_fname, "model14-alg", allergen_color_str)
    configRenderer.render(modelConfig_16, allergen_pdb_fname, "model16-alg", allergen_color_str)
    configRenderer.render(modelConfig_18, allergen_pdb_fname, "model18-alg", allergen_color_str)
    configRenderer.final_render()  ## hide all, then show as cartoons

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
##test_outputAllConfigsInLastStep()
#test_ConfigurationRenderer()
test_ConfigurationRenderer_manyModels()
