import sys
import os
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from database.chempath_database import main as db_main
from database.visualization import generate_visualizations



if __name__ == "__main__":
    db_main()
    generate_visualizations()

