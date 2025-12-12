class GapTransformationOperator:
    
    def __init__(self, name, template, scoring_fn):
        self.name = name
        self.template = template
        self.scoring_fn = scoring_fn
    
    def calculate_gap(self, model, evidence):
        return self.scoring_fn(model, evidence)
    
    def apply(self, model, evidence):
        pass