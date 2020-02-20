from megrim.toolbox import MegrimPlugin

class BaseModifications(MegrimPlugin):
    
    def __init__(self):
        super().__init__()
        self.tool = "BaseModifications"