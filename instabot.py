from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC

# Set the web driver and login URL
webdriver = webdriver.Chrome()
login_url = "https://www.instagram.com/accounts/login/"

# Use Selenium to open the web browser and navigate to the login page
webdriver.get(login_url)

# Use Selenium to enter your username and password into the login form and submit it
username_input = webdriver.find_element_by_name("username")
password_input = webdriver.find_element_by_name("password")
username_input.send_keys("your_username")
password_input.send_keys("your_password")
password_input.submit()

# Wait for the page to load and navigate to the DM inbox
WebDriverWait(webdriver, 10).until(EC.presence_of_element_located((By.CSS_SELECTOR, "a[href='/direct/inbox/']")))
webdriver.find_element_by_css_selector("a[href='/direct/inbox/']").click()

# Scroll through the DM inbox and download the posts that have been sent to you
posts = []
while True:
    # Wait for the page to load
    WebDriverWait(webdriver, 10).until(EC.presence_of_element_located((By.CSS_SELECTOR, "div[role='presentation']")))

    # Extract the data for the posts in the DM inbox
    elements = webdriver.find_elements_by_css_selector("div[role='presentation'] > div > div > div")

    for element in elements:
        # Extract the image or video file and any accompanying text or captions
        image_src = element.find_element_by_css_selector("img").get_attribute("src")
        text = element.find_element_by_css_selector("span").text

        # Save the post data to a list
        posts.append({
            "image_src": image_src,
            "text": text
        })

    # Scroll to the bottom of the page to load more posts
    webdriver.execute_script("window.scrollTo(0, document.body.scrollHeight);")

    # Check if there are more posts to load
    try:
        webdriver.find_element_by_css_selector("a[href='/direct/new/']")
        break
    except:
        pass

# Save the downloaded posts to a file or do something else with them
for post in posts:
    # Save the post to a file or do something else with it
    pass


import instabot

# Authenticate with Instagram using your API key
bot = instabot.Bot()
bot.login(username="your_username", password="your_password", api_key="your_api_key")

# Use the instabot library to access the Instagram API and download the posts that have been sent to you via DM
posts = bot.get_direct_pending_posts()

# Save the downloaded posts to a file or do something else with them
for post in posts:
    # Save the post to a file or do something else with it
    pass



import requests

# Set the API endpoint and your API key
api_endpoint = "https://i.instagram.com/api/v1/"
api_key = "your_api_key"

# Authenticate with Instagram and get an access token
response = requests.post(api_endpoint + "accounts/login/", json={
    "username": "your_username",
    "password": "your_password",
    "api_key": api_key
})

# Extract the access token from the response
access_token = response.json()["access_token"]

# Use the access token to make a request to the direct_v2/inbox/ endpoint to get your DM inbox
response = requests.get(api_endpoint + "direct_v2/inbox/", headers={
    "Authorization": "Bearer " + access_token
})

# Extract the posts from the response
posts = response.json()["inbox"]["threads"]

# Save the downloaded posts to a file or do something else with them
for post in posts:
    # Save the post to a file or do something else with it
    pass



