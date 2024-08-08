from database_operations import get_session, Compound
from optimization_algorithms import optimize_compound
from retrosynthetic_tools import analyze_retrosynthesis

class WorkflowManager:
    def __init__(self):
        self.workflows = {}
        self.session = get_session()

    async def start_workflow(self, compound_id):
        workflow = {
            'status': 'in_progress',
            'steps': ['retrieval', 'optimization', 'retrosynthesis'],
            'current_step': 'retrieval',
            'errors': [],
            'results': {}
        }
        self.workflows[compound_id] = workflow
        await self.execute_workflow(compound_id)

    async def execute_workflow(self, compound_id):
        workflow = self.workflows[compound_id]
        compound = self.session.query(Compound).get(compound_id)

        for step in workflow['steps']:
            try:
                if step == 'retrieval':
                    workflow['results']['original'] = compound.to_dict()
                elif step == 'optimization':
                    optimized_smiles = optimize_compound(compound.smiles)
                    compound.smiles = optimized_smiles
                    self.session.commit()
                    workflow['results']['optimized'] = compound.to_dict()
                elif step == 'retrosynthesis':
                    retro_result = analyze_retrosynthesis(compound.smiles)
                    workflow['results']['retrosynthesis'] = retro_result
                workflow['current_step'] = step
            except Exception as e:
                workflow['errors'].append(f"Error in {step}: {str(e)}")
                workflow['status'] = 'error'
                return
        workflow['status'] = 'completed'

    def get_workflow_status(self, compound_id):
        return self.workflows.get(compound_id, {'status': 'not_found'})

workflow_manager = WorkflowManager()
