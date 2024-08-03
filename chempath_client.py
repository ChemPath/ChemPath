import requests
import concurrent.futures
import random
import asyncio
import aiohttp

async def send_request(session, url, data):
    try:
        async with session.post(url, json=data) as response:
            response.raise_for_status()
            return await response.json()
    except aiohttp.ClientError as e:
        print(f"Error: {e}")
        return None

async def simulate_user(session, user_id):
    url = "http://localhost:5000/api/process_request"
    data = {
        "request_type": "search",
        "query": f"compound_{random.randint(1, 100)}"
    }
    result = await send_request(session, url, data)
    return f"User {user_id}: {result}"

async def simulate_users(num_users):
    async with aiohttp.ClientSession() as session:
        tasks = [simulate_user(session, i) for i in range(num_users)]
        results = await asyncio.gather(*tasks)
        for result in results:
            print(result)

async def main():
    await simulate_users(100)

if __name__ == "__main__":
    asyncio.run(main())
