# auth.py

import asyncio

# Simulated user database
users = {
    'admin': {'password': 'admin123', 'role': 'admin'},
    'user': {'password': 'user123', 'role': 'user'},
}

async def check_credentials(username, password):
    await asyncio.sleep(0.1)  # Simulate database lookup
    if username in users and users[username]['password'] == password:
        return True
    return False

async def get_user_role(username):
    await asyncio.sleep(0.1)  # Simulate database lookup
    if username in users:
        return users[username]['role']
    return 'guest'
