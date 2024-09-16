from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.chrome.service import Service as ChromeService
from selenium.webdriver.chrome.options import Options
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC

# Configure WebDriver
def configure_webdriver():
    chrome_driver_path = 'C:/WebDriver/chromedriver-win64/chromedriver.exe'
    chrome_options = Options()
    chrome_options.add_argument('--headless')  # Optional: Run in headless mode
    chrome_options.add_argument('--ignore-certificate-errors')
    chrome_options.add_argument('--ignore-ssl-errors')
    
    service = ChromeService(executable_path=chrome_driver_path)
    return webdriver.Chrome(service=service, options=chrome_options)

# Function to scrape herb details
def scrape_herb_details(herb_name):
    driver = configure_webdriver()
    url = f"https://old.tcmsp-e.com/tcmspsearch.php?qr={herb_name}&qsr=herb_en_name&token=b5fe4654ca6f6f1e7cd4ab4f35cd37cb"
    driver.get(url)

    try:
        # Attempt to find and click the 'Ingredients' tab
        try:
            ingredients_tab = WebDriverWait(driver, 10).until(
                EC.element_to_be_clickable((By.XPATH, "//li[@role='tab' and .//a[text()='Ingredients']]"))
            )
            ingredients_tab.click()
        except:
            print(f"No 'Ingredients' tab found for herb: {herb_name}.")
            driver.quit()
            return [], [], []  # Return empty lists if the tab is not found

        # Wait for the table content to load
        table_rows = WebDriverWait(driver, 10).until(
            EC.presence_of_all_elements_located((By.XPATH, "//table//tr"))
        )

        # Extract non-empty table content
        compounds = []
        for row in table_rows:
            columns = row.find_elements(By.TAG_NAME, "td")
            row_data = [col.text.strip() for col in columns]
            
            # Filter out empty rows
            if any(row_data):
                # Add logic to parse row_data into a structured dictionary
                compound_data = {
                    'name': row_data[0],  # Adjust based on table structure
                    'molecular_weight': row_data[2],
                    'alogp': row_data[3],
                    # Add more fields as needed
                }
                compounds.append(compound_data)

        # You would need to implement logic to scrape related targets and diseases
        related_targets = []  # Placeholder
        related_diseases = []  # Placeholder

        return compounds, related_targets, related_diseases

    except Exception as e:
        print(f"An error occurred while scraping: {e}")
        return [], [], []

    finally:
        # Close the WebDriver
        driver.quit()
