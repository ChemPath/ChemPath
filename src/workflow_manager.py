import asyncio
from database_operations import session, Compound
from optimization_algorithms import optimize_compound
from retrosynthetic_tools import analyze_retrosynthesis
from optimization_algorithms import optimize_compound, refine_based_on_retrosynthesis


class WorkflowManager:
    def __init__(self):
        self.workflows = {}

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
        compound = session.query(Compound).get(compound_id)

        for iteration in range(3):  # Allow up to 3 iterations
            try:
                optimized_smiles = optimize_compound(compound.smiles)
                retro_result = analyze_retrosynthesis(optimized_smiles)
                refined_smiles = refine_based_on_retrosynthesis(optimized_smiles, retro_result)
                
                compound.smiles = refined_smiles
                session.commit()
                
                workflow['results'][f'iteration_{iteration}'] = {
                    'optimized': optimized_smiles,
                    'retrosynthesis': retro_result,
                    'refined': refined_smiles
                }
                
                if refined_smiles == optimized_smiles:
                    break  # No further refinement needed
            except Exception as e:
                workflow['errors'].append(f"Error in iteration {iteration}: {str(e)}")
                workflow['status'] = 'error'
                return
        
        workflow['status'] = 'completed'