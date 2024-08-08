from src.database.chempath_database import create_engine_and_session, get_compound_by_smiles
from src.plant_retrosynthesis import PlantRetrosynthesis
from optimization import optimize_compound
from sqlalchemy.orm import Session
from typing import List, Dict, Any, Tuple
from rdkit import Chem
from rdkit.Chem import Descriptors, SDWriter
from src.database.chempath_database import Compound
import csv
import io
from sqlalchemy.exc import SQLAlchemyError
from .models import Compound
from aiohttp import web
import asyncio
from database_operations import session, Compound, update_compound, delete_compound
from sqlalchemy.exc import SQLAlchemyError
from optimization_algorithms import optimize_compound
from retrosynthetic_tools import analyze_retrosynthesis
from workflow_manager import workflow_manager

class ChemPathAPI:
    """
    The ChemPath API class provides a interface to interact with the ChemPath database.
    """

    def __init__(self, db_path: str):
        """
        Initialize the ChemPath API with a database path.

        Args:
            db_path (str): The path to the ChemPath database.
        """
        self.engine, self.session = create_engine_and_session(db_path)
        self.retro = PlantRetrosynthesis()

    # Validation methods

    def validate_smiles(self, smiles: str) -> bool:
        """
        Validate a SMILES string.

        Args:
            smiles (str): The SMILES string to validate.

        Returns:
            bool: True if the SMILES string is valid, False otherwise.
        """
        mol = Chem.MolFromSmiles(smiles)
        return mol is not None

    def validate_compound_data(self, compound_data: Dict[str, Any]) -> bool:
        """
        Validate compound data.

        Args:
            compound_data (Dict[str, Any]): The compound data to validate.

        Returns:
            bool: True if the compound data is valid, False otherwise.
        """
        required_fields = ['name', 'smiles', 'molecular_weight', 'logp']
        return all(field in compound_data for field in required_fields)

    # Compound methods

    def get_compound_by_id(self, compound_id: int) -> Dict[str, Any]:
        """
        Get a compound by its ID.

        Args:
            compound_id (int): The ID of the compound to retrieve.

        Returns:
            Dict[str, Any]: The compound data.
        """
        compound = self.session.query(Compound).get(compound_id)
        if not compound:
            raise ValueError("Compound not found")
        return {
            'id': compound.id,
            'name': compound.name,
            'smiles': compound.smiles,
            'molecular_weight': compound.molecular_weight,
            'logp': compound.logp,
            'plant_source': compound.plant_source
        }

    def batch_add_compounds(self, compounds_data: List[Dict[str, Any]]) -> List[int]:
        """
        Add multiple compounds to the database.

        Args:
            compounds_data (List[Dict[str, Any]]): The list of compound data to add.

        Returns:
            List[int]: The IDs of the added compounds.
        """
        new_compounds = []
        for compound_data in compounds_data:
            new_compound = Compound(**compound_data)
            self.session.add(new_compound)
            new_compounds.append(new_compound)
        self.session.commit()
        return [compound.id for compound in new_compounds]

    def batch_update_compounds(self, compounds_data: List[Dict[str, Any]]) -> List[int]:
        """
        Update multiple compounds in the database.

        Args:
            compounds_data (List[Dict[str, Any]]): The list of compound data to update.

        Returns:
            List[int]: The IDs of the updated compounds.
        """
        updated_ids = []
        for compound_data in compounds_data:
            compound_id = compound_data.pop('id', None)
            if compound_id:
                compound = self.session.query(Compound).get(compound_id)
                if compound:
                    for key, value in compound_data.items():
                        setattr(compound, key, value)
                    updated_ids.append(compound_id)
        self.session.commit()
        return updated_ids

    def batch_delete_compounds(self, compound_ids: List[int]) -> List[int]:
        """
        Delete multiple compounds from the database.

        Args:
            compound_ids (List[int]): The IDs of the compounds to delete.

        Returns:
            List[int]: The IDs of the deleted compounds.
        """
        deleted_ids = []
        for compound_id in compound_ids:
            compound = self.session.query(Compound).get(compound_id)
            if compound:
                self.session.delete(compound)
                deleted_ids.append(compound_id)
        self.session.commit()
        return deleted_ids

    # Search methods

    def search_compounds(self, query: str, page: int = 1, per_page: int = 10) -> Tuple[List[Dict[str, Any]], int]:
        """
        Search for compounds based on a query.

        Args:
            query (str): The search query.
            page (int, optional): The page number. Defaults to 1.
            per_page (int, optional): The number of results per page. Defaults to 10.

        Returns:
            Tuple[List[Dict[str, Any]], int]: The search results and the total count.
        """
        if not query or len(query) < 3:
            raise ValueError("Search query must be at least 3 characters long")

        search = f"%{query}%"
        base_query = self.session.query(Compound).filter(
            (Compound.name.ilike(search)) |
            (Compound.smiles.ilike(search)) |
            (Compound.plant_source.ilike(search))
        )

        total_count = base_query.count()

        offset = (page - 1) * per_page
        compounds = base_query.order_by(Compound.name).offset(offset).limit(per_page).all()

        result = [
            {
                'id': c.id,
                'name': c.name,
                'smiles': c.smiles,
                'molecular_weight': c.molecular_weight,
                'logp': c.logp,
                'plant_source': c.plant_source
            }
            for c in compounds
        ]

        return result, total_count

    # Export methods

    def export_compounds_to_sdf(self, compounds: List[Dict[str, Any]]) -> str:
        """
        Export compounds to an SDF file.

        Args:
            compounds (List[Dict[str, Any]]): The list of compound data to export.

        Returns:
            str: The SDF file contents.
        """
        output = io.StringIO()
        writer = SDWriter(output)
        for compound in compounds:
            mol = Chem.MolFromSmiles(compound['smiles'])
            for key, value in compound.items():
                if key != 'smiles':
                    mol.SetProp(key, str(value))
            writer.write(mol)
        writer.close()
        return output.getvalue()

    def export_compounds_to_csv(self, compounds: List[Dict[str, Any]]) -> str:
        """
        Export compounds to a CSV file.

        Args:
            compounds (List[Dict[str, Any]]): The list of compound data to export.

        Returns:
            str: The CSV file contents.
        """
        output = io.StringIO()
        writer = csv.DictWriter(output, fieldnames=compounds[0].keys())
        writer.writeheader()
        for compound in compounds:
            writer.writerow(compound)
        return output.getvalue()

    # Import methods

    def import_compounds_from_sdf(self, sdf_data: str) -> List[Dict[str, Any]]:
        """
        Import compounds from an SDF file.

        Args:
            sdf_data (str): The SDF file contents.

        Returns:
            List[Dict[str, Any]]: The list of imported compound data.
        """
        compounds = []
        supplier = Chem.SDMolSupplier()
        supplier.SetData(sdf_data)
        for mol in supplier:
            if mol is not None:
                compound = {
                    'smiles': Chem.MolToSmiles(mol),
                    'name': mol.GetProp('_Name') if mol.HasProp('_Name') else '',
                    'molecular_weight': Descriptors.ExactMolWt(mol),
                    'logp': Descriptors.MolLogP(mol)
                }
                for prop_name in mol.GetPropNames():
                    compound[prop_name] = mol.GetProp(prop_name)
                compounds.append(compound)
        return compounds

    def import_compounds_from_csv(self, csv_data: str) -> List[Dict[str, Any]]:
        """
        Import compounds from a CSV file.

        Args:
            csv_data (str): The CSV file contents.

        Returns:
            List[Dict[str, Any]]: The list of imported compound data.
        """
        compounds = []
        reader = csv.DictReader(io.StringIO(csv_data))
        for row in reader:
            if 'smiles' in row:
                mol = Chem.MolFromSmiles(row['smiles'])
                if mol is not None:
                    compound = {
                        'smiles': row['smiles'],
                        'name': row.get('name', ''),
                        'molecular_weight': Descriptors.ExactMolWt(mol),
                        'logp': Descriptors.MolLogP(mol)
                    }
                    for key, value in row.items():
                        if key not in compound:
                            compound[key] = value
                    compounds.append(compound)
        return compounds

    # Optimization methods

    def optimize_structure(self, smiles: str) -> str:
        """
        Optimize the structure of a compound.

        Args:
            smiles (str): The SMILES string of the compound to optimize.

        Returns:
            str: The optimized SMILES string.
        """
        if not self.validate_smiles(smiles):
            raise ValueError("Invalid SMILES string")
        return optimize_compound(self.session, smiles)

    # Retrosynthesis methods

    def perform_retrosynthesis(self, smiles: str) -> Dict[str, Any]:
        """
        Perform retrosynthesis on a compound.

        Args:
            smiles (str): The SMILES string of the compound to perform retrosynthesis on.

        Returns:
            Dict[str, Any]: The retrosynthesis results.
        """
        if not self.validate_smiles(smiles):
            raise ValueError("Invalid SMILES string")
        return self.retro.analyze(smiles)
    # Utility methods

    def get_compound_by_smiles(self, smiles: str) -> Dict[str, Any]:
        """
        Get a compound by its SMILES string.

        Args:
            smiles (str): The SMILES string of the compound to retrieve.

        Returns:
            Dict[str, Any]: The compound data.
        """
        compound = self.session.query(Compound).filter_by(smiles=smiles).first()
        if not compound:
            raise ValueError("Compound not found")
        return {
            'id': compound.id,
            'name': compound.name,
            'smiles': compound.smiles,
            'molecular_weight': compound.molecular_weight,
            'logp': compound.logp,
            'plant_source': compound.plant_source
        }

    def get_compounds_by_plant_source(self, plant_source: str) -> List[Dict[str, Any]]:
        """
        Get compounds by their plant source.

        Args:
            plant_source (str): The plant source to retrieve compounds for.

        Returns:
            List[Dict[str, Any]]: The list of compound data.
        """
        compounds = self.session.query(Compound).filter_by(plant_source=plant_source).all()
        return [
            {
                'id': compound.id,
                'name': compound.name,
                'smiles': compound.smiles,
                'molecular_weight': compound.molecular_weight,
                'logp': compound.logp,
                'plant_source': compound.plant_source
            }
            for compound in compounds
        ]

    def get_plant_sources(self) -> List[str]:
        """
        Get all plant sources.

        Returns:
            List[str]: The list of plant sources.
        """
        plant_sources = self.session.query(Compound.plant_source).distinct().all()
        return [plant_source[0] for plant_source in plant_sources]

    def get_compound_count(self) -> int:
        """
        Get the total number of compounds.

        Returns:
            int: The total number of compounds.
        """
        return self.session.query(Compound).count()

    def get_plant_source_count(self) -> int:
        """
        Get the total number of plant sources.

        Returns:
            int: The total number of plant sources.
        """
        return self.session.query(Compound.plant_source).distinct().count()
    def get_plant_compound_by_smiles(self, smiles: str) -> Dict[str, Any]:
        """
        Get a plant compound by its SMILES string.

        Args:
            smiles (str): The SMILES string of the compound to retrieve.

        Returns:
            Dict[str, Any]: The compound data.
        """
        if not self.validate_smiles(smiles):
            raise ValueError("Invalid SMILES string")
        return get_compound_by_smiles(self.session, smiles)

    def add_plant_compound(self, compound_data: Dict[str, Any]) -> int:
        """
        Add a new plant compound.

        Args:
            compound_data (Dict[str, Any]): The compound data to add.

        Returns:
            int: The ID of the added compound.
        """
        if not self.validate_compound_data(compound_data):
            raise ValueError("Invalid compound data")
        # Implement the add logic here
        pass

    def update_plant_compound(self, compound_id: int, compound_data: Dict[str, Any]) -> bool:
        """
        Update an existing plant compound.

        Args:
            compound_id (int): The ID of the compound to update.
            compound_data (Dict[str, Any]): The updated compound data.

        Returns:
            bool: True if the update was successful, False otherwise.
        """
        if compound_id < 1:
            raise ValueError("Invalid compound ID")
        if not self.validate_compound_data(compound_data):
            raise ValueError("Invalid compound data")
        # Implement the update logic here
        pass

    def get_biosynthetic_pathways(self) -> List[Dict[str, Any]]:
        """
        Get all biosynthetic pathways.

        Returns:
            List[Dict[str, Any]]: The list of pathway data.
        """
        # Implement the retrieval logic here
        pass

    def get_compounds_by_pathway(self, pathway_id: int) -> List[Dict[str, Any]]:
        """
        Get compounds by their pathway ID.

        Args:
            pathway_id (int): The ID of the pathway to retrieve compounds for.

        Returns:
            List[Dict[str, Any]]: The list of compound data.
        """
        if pathway_id < 1:
            raise ValueError("Invalid pathway ID")
        # Implement the retrieval logic here
        pass

    def optimize_plant_compound(self, smiles: str, target_property: str, iterations: int = 100) -> str:
        """
        Optimize a plant compound.

        Args:
            smiles (str): The SMILES string of the compound to optimize.
            target_property (str): The target property to optimize for.
            iterations (int, optional): The number of iterations. Defaults to 100.

        Returns:
            str: The optimized SMILES string.
        """
        if not self.validate_smiles(smiles):
            raise ValueError("Invalid SMILES string")
        if target_property not in ['logP', 'molecular_weight', 'polararea', 'complexity']:
            raise ValueError("Invalid target property")
        if iterations < 1:
            raise ValueError("Number of iterations must be positive")
        return optimize_compound(self.session, smiles, target_property, iterations)

    def perform_plant_retrosynthesis(self, smiles: str) -> Dict[str, Any]:
        """
        Perform retrosynthesis on a plant compound.

        Args:
            smiles (str): The SMILES string of the compound to perform retrosynthesis on.

        Returns:
            Dict[str, Any]: The retrosynthesis results.
        """
        if not self.validate_smiles(smiles):
            raise ValueError("Invalid SMILES string")
        return self.retro.analyze(smiles)

    def suggest_plant_precursors(self, smiles: str) -> List[Dict[str, Any]]:
        """
        Suggest precursors for a plant compound.

        Args:
            smiles (str): The SMILES string of the compound to suggest precursors for.

        Returns:
            List[Dict[str, Any]]: The list of precursor data.
        """
        if not self.validate_smiles(smiles):
            raise ValueError("Invalid SMILES string")
        return self.retro.suggest_precursors(smiles)

    def close(self):
        """
        Close the database session.
        """
        self.session.close()

    def optimize_compounds(self, compound_data: List[Dict[str, str]]) -> List[Dict[str, str]]:
        """
        Optimize multiple compounds.

        Args:
            compound_data (List[Dict[str, str]]): The list of compound data to optimize.

        Returns:
            List[Dict[str, str]]: The list of optimized compound data.
        """
        if not compound_data:
            raise ValueError("Compound data list is empty")
        for compound in compound_data:
            if not self.validate_smiles(compound.get('smiles', '')):
                raise ValueError(f"Invalid SMILES string: {compound.get('smiles', '')}")
        # Implement the optimization logic here
        pass   

class WorkflowManager:
    def __init__(self):
        self.workflows = {}

    async def start_workflow(self, compound_id):
        workflow = {
            'status': 'in_progress',
            'steps': ['retrieval', 'optimization', 'retrosynthesis', 'refinement'],
            'current_step': 'retrieval',
            'errors': []
        }
        self.workflows[compound_id] = workflow
        await self.execute_workflow(compound_id)

    async def execute_workflow(self, compound_id):
        workflow = self.workflows[compound_id]
        for step in workflow['steps']:
            try:
                if step == 'retrieval':
                    await get_compound(web.Request.from_dict({'match_info': {'compound_id': compound_id}}))
                elif step == 'optimization':
                    await optimize_compound_structure(web.Request.from_dict({'match_info': {'compound_id': compound_id}}))
                elif step == 'retrosynthesis':
                    await analyze_compound_retrosynthesis(web.Request.from_dict({'match_info': {'compound_id': compound_id}}))
                elif step == 'refinement':
                    await refine_compound(web.Request.from_dict({'match_info': {'compound_id': compound_id}, 'json': {'refined_smiles': 'C1=CC=CC=C1'}}))
                workflow['current_step'] = step
            except Exception as e:
                workflow['errors'].append(f"Error in {step}: {str(e)}")
                workflow['status'] = 'error'
                return
        workflow['status'] = 'completed'

workflow_manager = WorkflowManager()

routes = web.RouteTableDef()

@routes.get('/compounds/{compound_id}')
async def get_compound(request):
    compound_id = int(request.match_info['compound_id'])
    try:
        compound = session.query(Compound).filter(Compound.id == compound_id).first()
        if compound is None:
            return web.json_response({"error": "Compound not found"}, status=404)
        return web.json_response(compound.to_dict())
    except SQLAlchemyError as e:
        return web.json_response({"error": f"Database error: {str(e)}"}, status=500)

@routes.get('/compounds')
async def get_compounds(request):
    query = request.query.get('query', '')
    page = int(request.query.get('page', 1))
    per_page = int(request.query.get('per_page', 10))
    min_mw = float(request.query.get('min_mw')) if 'min_mw' in request.query else None
    max_mw = float(request.query.get('max_mw')) if 'max_mw' in request.query else None
    min_logp = float(request.query.get('min_logp')) if 'min_logp' in request.query else None
    max_logp = float(request.query.get('max_logp')) if 'max_logp' in request.query else None
    plant_source = request.query.get('plant_source')
    
@routes.post('/compounds/batch')
async def batch_operations(request):
    data = await request.json()
    operation = data.get('operation')
    compounds_data = data.get('compounds', [])
    
    api = ChemPathAPI("chempath_database.db")
    
    if operation == 'add':
        result = api.batch_add_compounds(compounds_data)
        return web.json_response({'added_ids': result})
    elif operation == 'update':
        result = api.batch_update_compounds(compounds_data)
        return web.json_response({'updated_ids': result})
    elif operation == 'delete':
        result = api.batch_delete_compounds(compounds_data)
        return web.json_response({'deleted_ids': result})
    else:
        return web.json_response({'error': 'Invalid operation'}, status=400)
    
    api = ChemPathAPI("chempath_database.db")
    compounds, total_count = api.search_compounds(query, page, per_page, 
                                                  min_mw, max_mw, min_logp, max_logp, plant_source)
    
    return web.json_response({
        'compounds': compounds,
        'total_count': total_count,
        'page': page,
        'per_page': per_page
    })


@routes.post('/compounds')
async def create_compound(request):
    data = await request.json()
    try:
        new_compound = Compound(**data)
        session.add(new_compound)
        session.commit()
        return web.json_response(new_compound.to_dict(), status=201)
    except SQLAlchemyError as e:
        session.rollback()
        return web.json_response({"error": f"Error creating compound: {str(e)}"}, status=400)

@routes.put('/compounds/{compound_id}')
async def update_compound_route(request):
    compound_id = int(request.match_info['compound_id'])
    data = await request.json()
    try:
        success = update_compound(session, compound_id, data)
        if success:
            return web.json_response({"message": "Compound updated successfully"})
        else:
            return web.json_response({"error": "Compound not found"}, status=404)
    except SQLAlchemyError as e:
        return web.json_response({"error": f"Error updating compound: {str(e)}"}, status=400)

@routes.delete('/compounds/{compound_id}')
async def delete_compound_route(request):
    compound_id = int(request.match_info['compound_id'])
    try:
        success = delete_compound(session, compound_id)
        if success:
            return web.json_response({"message": "Compound deleted successfully"})
        else:
            return web.json_response({"error": "Compound not found"}, status=404)
    except SQLAlchemyError as e:
        return web.json_response({"error": f"Error deleting compound: {str(e)}"}, status=400)

@routes.post('/compounds/{compound_id}/optimize')
async def optimize_compound_structure(request):
    compound_id = int(request.match_info['compound_id'])
    compound = session.query(Compound).get(compound_id)
    optimized_structure = optimize_compound(compound.smiles)
    compound.smiles = optimized_structure
    session.commit()
    return web.json_response(compound.to_dict())

@routes.get('/compounds/{compound_id}/retrosynthesis')
async def analyze_compound_retrosynthesis(request):
    compound_id = int(request.match_info['compound_id'])
    compound = session.query(Compound).get(compound_id)
    retrosynthesis_result = analyze_retrosynthesis(compound.smiles)
    return web.json_response(retrosynthesis_result)

@routes.put('/compounds/{compound_id}/refine')
async def refine_compound(request):
    compound_id = int(request.match_info['compound_id'])
    data = await request.json()
    compound = session.query(Compound).get(compound_id)
    compound.smiles = data['refined_smiles']
    session.commit()
    return web.json_response(compound.to_dict())


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

        for step in workflow['steps']:
            try:
                if step == 'retrieval':
                    workflow['results']['original'] = compound.to_dict()
                elif step == 'optimization':
                    optimized_smiles = optimize_compound(compound.smiles)
                    compound.smiles = optimized_smiles
                    session.commit()
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

# Add these routes to your existing routes in the main API file:

@routes.post('/workflows/{compound_id}')
async def start_compound_workflow(request):
    compound_id = int(request.match_info['compound_id'])
    await workflow_manager.start_workflow(compound_id)
    return web.json_response({'message': 'Workflow started'})

@routes.get('/workflows/{compound_id}')
async def get_workflow_status(request):
    compound_id = int(request.match_info['compound_id'])
    status = workflow_manager.get_workflow_status(compound_id)
    return web.json_response(status)

routes.append(start_compound_workflow)
routes.append(get_workflow_status)

async def main():
    app = web.Application()
    app.add_routes(routes)
    return app

if __name__ == "__main__":
    asyncio.run(main())