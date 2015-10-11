import Task

class CastName(Task):
    def __init__(self, inname, outname, tpe, *args, **kwargs):
        super(CastName, self).__init__(*args, **kwargs)
        self.local=True
        self.cache=False
        self.tpe = tpe
        self.inname = inname
        self.outname = outname
    def inputs(self):
        return [(self.inname, self.tpe)]
    def outputs(self):
        return [(self.outname, self.tpe)]
    def run(self):
        self.result =  {self.outname : self.input_data[self.inname]}
        return self.result
    
class Exit(Task):
    def __init__(self, *args, **kwargs):
        super(Exit, self).__init__(*args, **kwargs)
        self.EXIT = True
    def inputs(self):
        return []
    def outputs(self):
        return []
