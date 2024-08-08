from aiohttp import web
from database_operations import get_session, Compound, update_compound, delete_compound
from sqlalchemy.exc import SQLAlchemyError
from optimization_algorithms import optimize_compound
from retrosynthetic_tools import analyze_retrosynthesis
from workflow_manager import workflow_manager

# Create a session
session = get_session()
routes = web.RouteTableDef()

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
