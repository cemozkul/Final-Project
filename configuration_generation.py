import yaml
import os

# Load the original configuration file
with open('config.yaml') as file:
    config = yaml.load(file, Loader=yaml.FullLoader)

# Define the list of desired numbers of cities
num_cities_list = [10, 20, 30, 40, 50, 75, 100]

# Loop over the list and modify the configuration file for each number of cities
for num_cities in num_cities_list:
    config['cities'][0]['numberOfCities'] = num_cities
    # Save the modified configuration file with the corresponding name
    filename = f'config{num_cities}.yaml'
    with open(os.path.join('configWithCities', filename), 'w') as file:
        yaml.dump(config, file)
