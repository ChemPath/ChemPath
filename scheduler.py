from apscheduler.schedulers.background import BackgroundScheduler
from ml_model import train_model
from your_database_module import fetch_latest_data_from_database

def retrain_model():
    # Fetch the latest data from the database
    df = fetch_latest_data_from_database()
    all_therapeutic_areas = list(set([area for areas in df['therapeutic_areas'] for area in areas]))

    # Train the model
    train_model(df, all_therapeutic_areas)

    print("Model retrained successfully")

def start_scheduler():
    scheduler = BackgroundScheduler()
    scheduler.add_job(retrain_model, 'interval', weeks=1)
    scheduler.start()
