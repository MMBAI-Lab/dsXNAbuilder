"""Screenshot of the base-pair grid and warnings for the LNA gapmer example (light theme)."""
import os, subprocess, sys, time
from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.firefox.options import Options
from selenium.webdriver.firefox.service import Service
from selenium.webdriver.support.ui import WebDriverWait
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", ".."))
port = 8792
srv = subprocess.Popen([sys.executable, os.path.join(ROOT, "scripts", "serve.py"), str(port)],
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
time.sleep(1)
o = Options(); o.add_argument("-headless")
d = webdriver.Firefox(options=o, service=Service("/snap/bin/geckodriver"))
d.set_window_size(1500, 1200)
try:
    d.get(f"http://localhost:{port}/")
    d.execute_script("localStorage.setItem('theme','light'); localStorage.setItem('dsxna-lang','en');")
    d.get(f"http://localhost:{port}/")
    WebDriverWait(d, 120).until(lambda d: d.find_element(By.ID, "model").find_elements(By.TAG_NAME, "option"))
    d.find_element(By.CSS_SELECTOR, '[data-ex="lna"]').click()
    WebDriverWait(d, 60).until(lambda d: d.find_element(By.TAG_NAME, "body").get_attribute("data-build-ms"))
    time.sleep(0.5)
    panel = d.find_element(By.ID, "grid").find_element(By.XPATH, "..")
    panel.screenshot(os.path.join(ROOT, "docs", "paper", "figures", "screen_grid_lna.png"))
    d.find_element(By.ID, "top").find_element(By.XPATH, "../..").screenshot(
        os.path.join(ROOT, "docs", "paper", "figures", "screen_input_lna.png"))
finally:
    d.quit(); srv.terminate()
print("ok")
