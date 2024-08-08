import sys
import os
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from aiohttp import web
from aiohttp_session import setup, get_session
from aiohttp_session.cookie_storage import EncryptedCookieStorage
import aiohttp_jinja2
import jinja2
from cryptography import fernet
from api_v1 import routes as v1_routes
from auth import check_credentials, get_user_role
from rate_limit import RateLimiter
from flask import Flask, request, jsonify
from src.database.chempath_database import db, PlantCompound
# Generate a random secret key
fernet_key = fernet.Fernet.generate_key()
secret_key = fernet_key[:32]

app = web.Application()

# Set up session middleware
setup(app, EncryptedCookieStorage(secret_key))

# Authentication middleware
@web.middleware
async def auth_middleware(request, handler):
    session = await get_session(request)
    request['user'] = session.get('user')
    return await handler(request)

# Authorization middleware
@web.middleware
async def authorize_middleware(request, handler):
    user = request.get('user')
    if user:
        request['role'] = await get_user_role(user)
    else:
        request['role'] = 'guest'
    return await handler(request)

# Rate limiting middleware
rate_limiter = RateLimiter(rate=5, per=1.0)  # 5 requests per second

@web.middleware
async def rate_limit_middleware(request, handler):
    if not rate_limiter.allow_request():
        raise web.HTTPTooManyRequests(text="Rate limit exceeded")
    return await handler(request)

# Set up middlewares
app.middlewares.extend([auth_middleware, authorize_middleware, rate_limit_middleware])

# Set up Jinja2 templates
aiohttp_jinja2.setup(app, loader=jinja2.FileSystemLoader('templates'))

# Version 1 routes
for route in v1_routes:
    if route.method == 'GET':
        app.router.add_get('/v1' + route.path, route.handler)
    elif route.method == 'POST':
        app.router.add_post('/v1' + route.path, route.handler)
    elif route.method == 'PUT':
        app.router.add_put('/v1' + route.path, route.handler)
    elif route.method == 'DELETE':
        app.router.add_delete('/v1' + route.path, route.handler)

# Login and logout handlers
@aiohttp_jinja2.template('login.html')
async def login_handler(request):
    if request.method == 'POST':
        data = await request.post()
        username = data['username']
        password = data['password']
        if await check_credentials(username, password):
            session = await get_session(request)
            session['user'] = username
            raise web.HTTPFound('/')
        else:
            return {'error': 'Invalid credentials'}
    return {}

async def logout_handler(request):
    session = await get_session(request)
    session.clear()
    raise web.HTTPFound('/')

# Add login and logout routes
app.router.add_get('/login', login_handler)
app.router.add_post('/login', login_handler)
app.router.add_get('/logout', logout_handler)

if __name__ == '__main__':
    web.run_app(app)

app = Flask(__name__)
app.config['SQLALCHEMY_DATABASE_URI'] = 'sqlite:///chempath.db'
db.init_app(app)

@app.route('/api/search', methods=['GET'])
def search():
    query = request.args.get('q', '')
    search_type = request.args.get('type', 'full_text')

    if search_type == 'full_text':
        results = PlantCompound.full_text_search(query)
    elif search_type == 'fuzzy':
        results = PlantCompound.fuzzy_search(query)
    else:
        return jsonify({'error': 'Invalid search type'}), 400

    return jsonify([{
        'id': r.id,
        'name': r.name,
        'plant_source': r.plant_source,
        'biological_activity': r.biological_activity
    } for r in results])

@app.route('/api/filter', methods=['GET'])
def filter():
    filter_params = {
        'min_molecular_weight': request.args.get('min_mw', type=float),
        'max_molecular_weight': request.args.get('max_mw', type=float),
        'min_logp': request.args.get('min_logp', type=float),
        'max_logp': request.args.get('max_logp', type=float),
        'plant_source': request.args.get('plant_source'),
        'biological_activity': request.args.get('biological_activity')
    }
    
    results = PlantCompound.advanced_filter(**{k: v for k, v in filter_params.items() if v is not None})

    return jsonify([{
        'id': r.id,
        'name': r.name,
        'molecular_weight': r.molecular_weight,
        'logp': r.logp,
        'plant_source': r.plant_source,
        'biological_activity': r.biological_activity
    } for r in results])

if __name__ == '__main__':
    app.run(debug=True)