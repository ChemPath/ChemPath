import time
import csv
import random
import asyncio
import aiohttp
from chempath_api import ChemPathAPI
from chempath_utils import validate_smiles
import psutil

async def generate_large_dataset(num_compounds):
    compounds = []
    for i in range(num_compounds):
        compound = {
            'name': f'Compound_{i}',
            'smiles': 'C' * random.randint(5, 20),  # Simplified SMILES generation
            'molecular_weight': random.uniform(100, 500),
            'logp': random.uniform(-2, 5),
            'plant_source': f'Plant_{i % 100}',
            'biological_activities': f'Activity_{i % 50}',
            'traditional_use': f'Use_{i % 30}',
            'h_bond_donors': random.randint(0, 5),
            'h_bond_acceptors': random.randint(0, 10),
            'polar_surface_area': random.uniform(0, 200),
            'rotatable_bonds': random.randint(0, 10)
        }
        compounds.append(compound)
    return compounds

async def test_high_volume_insert(api, compounds):
    start_time = time.time()
    tasks = [api.insert_compound(compound) for compound in compounds]
    await asyncio.gather(*tasks)
    end_time = time.time()
    print(f"Inserted {len(compounds)} compounds in {end_time - start_time:.2f} seconds")

async def test_concurrent_searches(api, num_searches):
    queries = ['C', 'O', 'N', 'S', 'P']
    filters = [
        {'molecular_weight': (100, 200)},
        {'logp': (0, 2)},
        {'plant_source': 'Plant_1'},
        {'biological_activities': 'Activity_1'},
        {'traditional_use': 'Use_1'}
    ]

    async def search_worker():
        query = random.choice(queries)
        filter_set = random.choice(filters)
        return await api.search_compounds(query=query, filters=filter_set)

    start_time = time.time()
    tasks = [search_worker() for _ in range(num_searches)]
    await asyncio.gather(*tasks)
    end_time = time.time()
    print(f"Completed {num_searches} concurrent searches in {end_time - start_time:.2f} seconds")

async def test_bulk_predictions(api, compounds):
    start_time = time.time()
    tasks = [api.predict_therapeutic_areas(compound['smiles']) for compound in compounds]
    await asyncio.gather(*tasks)
    end_time = time.time()
    print(f"Predicted therapeutic areas for {len(compounds)} compounds in {end_time - start_time:.2f} seconds")

def monitor_resources():
    cpu_percent = psutil.cpu_percent()
    memory_percent = psutil.virtual_memory().percent
    disk_io = psutil.disk_io_counters()
    print(f"CPU Usage: {cpu_percent}%")
    print(f"Memory Usage: {memory_percent}%")
    print(f"Disk I/O - Read: {disk_io.read_bytes / (1024*1024):.2f} MB, Write: {disk_io.write_bytes / (1024*1024):.2f} MB")

async def run_performance_tests():
    api = ChemPathAPI('chempath_test.db')
    
    compounds = await generate_large_dataset(10000)

    print("\nTesting high-volume insert...")
    await test_high_volume_insert(api, compounds[:3000])
    monitor_resources()

    print("\nTesting concurrent searches...")
    await test_concurrent_searches(api, 500)
    monitor_resources()

    print("\nTesting bulk predictions...")
    await test_bulk_predictions(api, compounds[3000:4000])
    monitor_resources()

    print("\nTesting combined high-load scenario...")
    start_time = time.time()
    await asyncio.gather(
        test_high_volume_insert(api, compounds[:3000]),
        test_concurrent_searches(api, 500),
        test_bulk_predictions(api, compounds[3000:4000])
    )
    end_time = time.time()
    print(f"Completed combined high-load scenario in {end_time - start_time:.2f} seconds")
    monitor_resources()

if __name__ == "__main__":
    asyncio.run(run_performance_tests())
