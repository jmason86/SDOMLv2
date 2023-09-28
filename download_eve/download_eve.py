import subprocess
import calendar

# Specify the base URL
base_url = "https://lasp.colorado.edu/eve/data_access/evewebdataproducts/level2b"

# Specify the years
years = [2023]

# Directory where wget should save the files
download_directory = '/Volumes/Roci Extension/sdo/eve/download/'

# Loop over each year
for year in years:
    # Determine the number of days in the current year
    is_leap_year = calendar.isleap(year)
    days_in_year = 366 if is_leap_year else 365

    # Loop over each day of the year
    for doy in range(261, days_in_year + 1):
        # Create the URL for the current year and day of the year
        url = f"{base_url}/{year}/{doy:03d}/"
        
        # Define the wget command as a list of strings
        wget_command = [
            "wget",
            "--no-parent",
            "-N",
            "-nH",
            "-r",
            "-nd",
            "-R", "*.html",
            "-e", "robots=off",
            "--mirror",
            "-P", download_directory,  # Specify the download directory
            url
        ]

        try:
            # Execute the wget command
            subprocess.run(wget_command, check=True)
        except subprocess.CalledProcessError as e:
            print(f"Error downloading {url}: {e}")
