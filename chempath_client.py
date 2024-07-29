import requests
import concurrent.futures
import random

def send_request(url, data):
    try:
        response = requests.post(url, json=data)
        response.raise_for_status()
        return response.json()
    except requests.exceptions.RequestException as e:
        print(f"Error: {e}")
        print(f"Response content: {response.text}")
        return None

def simulate_user(user_id):
    url = "http://localhost:5000/api/process_request"
    data = {
        "request_type": "search",
        "query": f"compound_{random.randint(1, 100)}"
    }
    result = send_request(url, data)
    return f"User {user_id}: {result}"

def simulate_users(num_users):
    with concurrent.futures.ThreadPoolExecutor(max_workers=10) as executor:
        futures = [executor.submit(simulate_user, i) for i in range(num_users)]
        for future in concurrent.futures.as_completed(futures):
            print(future.result())

if __name__ == "__main__":
    simulate_users(100)
