import yaml
import numpy as np
import pandas as pd
from datetime import datetime
import random
import os
import math


# reading yaml file and storing inside dictionary
with open(r'config.yaml') as file:
    dictionary = yaml.load(file, Loader=yaml.FullLoader)
   # for key, value in dictionary.items():
    #    print(key + " : " + str(value))

# taking necessary information from dictionary to generate random instances

# number of cities
numberOfCities = dictionary['cities'][0]['numberOfCities']

# coordinates
latitude_lower = dictionary['coordinates'][0]['latitude_lower']
latitude_upper = dictionary['coordinates'][0]['latitude_upper']
longitude_lower = dictionary['coordinates'][0]['longitude_lower']
longitude_upper = dictionary['coordinates'][0]['longitude_upper']


# information about fuel price related information
fuelPricesDistributions_lower = dictionary['fuelPricesDistributions'][0]['lower']
fuelPricesDistributions_upper = dictionary['fuelPricesDistributions'][0]['upper']
fuelPricesDistributions_mode = dictionary['fuelPricesDistributions'][0]['mode']

# information about toll price distribution
tollPricesDistribution_lower = dictionary['tollPricesDistribution'][0]['lower']
tollPricesDistribution_upper = dictionary['tollPricesDistribution'][0]['upper']

#information about fuel consumption distribution
fuelConsumptionDistribution_lower= dictionary['fuelConsumptionDistribution'][0]['lower']
fuelConsumptionDistribution_mode= dictionary['fuelConsumptionDistribution'][0]['mode']
fuelConsumptionDistribution_upper= dictionary['fuelConsumptionDistribution'][0]['upper']

#information about velocity distribution
averageVelocityDistribution_lower = dictionary['averageVelocityDistribution'][0]['lower']
averageVelocityDistribution_mode = dictionary['averageVelocityDistribution'][0]['mode']
averageVelocityDistribution_upper = dictionary['averageVelocityDistribution'][0]['upper']

#alpha - employee cost
alpha = dictionary['decisions'][0]['alpha']

#F - max fuel capacity
F = dictionary['decisions'][0]['F']

#Sf - starting fuel amount
Sf = dictionary['decisions'][0]['Sf']

#Earth radius in kilometers
R = dictionary['decisions'][0]['R']

#Set the seed value
seed = dictionary['decisions'][0]['seed']
np.random.seed(seed)
random.seed(seed)

#Limitation abput eliminating far cities
distance_cons = dictionary['decisions'][0]['distance_constraint']

def haversine(lat1, lon1, lat2, lon2):

    # Convert latitude and longitude to radians
    phi1 = math.radians(lat1)
    phi2 = math.radians(lat2)
    delta_phi = math.radians(lat2 - lat1)
    delta_lambda = math.radians(lon2 - lon1)

    # Calculate the Haversine formula
    a = math.sin(delta_phi / 2)**2 + \
        math.cos(phi1) * math.cos(phi2) * \
        math.sin(delta_lambda / 2)**2
    c = 2 * math.atan2(math.sqrt(a), math.sqrt(1 - a))

    # Calculate the distance between the two points
    distance = R * c
    return distance


def generate_random_places(n, lat_range, lon_range):
    latitudes = np.random.uniform(low=lat_range[0], high=lat_range[1], size=n)
    longitudes = np.random.uniform(low=lon_range[0], high=lon_range[1], size=n)
    return list(zip(latitudes, longitudes))

def calculate_distances_with_ids(cities):
    city_pairs_distances = []
    for i in range(len(cities)):
        for j in range(len(cities)):
            if i != j:
                city1 = cities[i]
                city2 = cities[j]
                city_id1 = city1['id']
                city_id2 = city2['id']
                distance = haversine(city1['latitude'], city1['longitude'], city2['latitude'], city2['longitude'])
                if distance > distance_cons:
                    distance = -1  # Set distance to -1 for pairs with distance > 1000
                city_pair_distance = [(city_id1, city_id2), distance]
                city_pairs_distances.append(city_pair_distance)

    return city_pairs_distances


cities_with_ids = []  #creates unique id's for every city

ID_list = []
city_id_map = {city['id']: index for index, city in enumerate(cities_with_ids)}
distances_with_city_ids = []
city_pairs_toll = []
toll = []
#paths = [[[] for j in range(numberOfCities)] for i in range(numberOfCities)]
#city_pairs_paths = []
#path_list = []
#pathDistance = [[[] for j in range(numberOfCities)] for i in range(numberOfCities)]
#city_pairs_pathDistance = []
#pathDistance_list = []
fuel_prices_array = []
fuelPrice = []
city_pairs_fuelC = []
fuelConsumption = []
city_pairs_velocity = []
velocity_list = []
city_pairs_time = []
time = []
total_costs = []


def generate_random():
    
    cities = generate_random_places(numberOfCities, (latitude_lower, latitude_upper), (longitude_lower, longitude_upper))

    id_counter = 1
    for city in cities:
        city_with_id = {'id': id_counter, 'latitude': city[0], 'longitude': city[1]}
        cities_with_ids.append(city_with_id)
        id_counter += 1

    for dictionary in cities_with_ids:
        ID_list.append(dictionary['id'])

    distances_with_ids = calculate_distances_with_ids(cities_with_ids)

    # Convert the (i, j) pairs in distances_with_ids to use city ids
    for (city_id1, city_id2), distance in distances_with_ids:
        distances_with_city_ids.append([city_id1, city_id2, round(distance, 2)])
    
    
    # Create the toll arrays
    tolls = np.zeros((numberOfCities, numberOfCities), dtype=object)

    for i in range(numberOfCities):
        for j in range(numberOfCities):
            # Determine the length of the arrays
            array_length2 = 1
            # Create arrays of that length for toll
            tolls[i][j] = np.around(np.random.uniform(tollPricesDistribution_lower, tollPricesDistribution_upper, size=array_length2), decimals=3)

    for i in range(numberOfCities):
        for j in range(numberOfCities):
            # the upper triangle of the array has to be equal to the lower triangle beacuse the toll between two city must be same
            if i < j:
                tolls[j][i] = tolls[i][j]
                
    for i in range(numberOfCities):
        for j in range(numberOfCities):
            if i == j:
                tolls[j][i] = np.zeros(1)



    for i in range(len(cities_with_ids)):
        for j in range(len(cities_with_ids)):
            if i != j:
                city1 = cities_with_ids[i]
                city2 = cities_with_ids[j]
                city_id1 = city1['id']
                city_id2 = city2['id']
                toll_values = tolls[i][j].tolist()
            

                for tollp in toll_values:  # iterate over each toll price separately
                                city_pair_toll = [(city_id1, city_id2), tollp]  # create a tuple with the city pair and a single toll price
                                city_pairs_toll.append(city_pair_toll)

                # Print the toll values for each city pair

    for (city_id1, city_id2), tollp in city_pairs_toll:  # iterate over the updated toll price list
        toll.append([city_id1, city_id2, tollp])  # append each city-toll price tuple to the toll price list
    # create a new 2D array of the same size as distance with 0 values 
                    
    '''                
    for i in range(len(cities_with_ids)):
        for j in range(len(cities_with_ids)):
            if i != j:
                city_id1 = cities_with_ids[i]['id']
                city_id2 = cities_with_ids[j]['id']
                paths_values = paths[i][j]
                for path in paths_values:  # iterate over each path separately
                    city_pair_paths = [(city_id1, city_id2), [x+1 for x in path]]  # create a tuple with the city pair and a single path
                    city_pairs_paths.append(city_pair_paths)


    # Print the toll values for each city pair
    for (city_id1, city_id2), path in city_pairs_paths:  # iterate over the updated city-paths list
        path_list.append([city_id1, city_id2, path])  # append each city-paths tuple to the path list


    for i in range(numberOfCities):
        for j in range(i+1, numberOfCities):
            path_distances = []
            for path in paths[i][j]:
                path_distance = 0
                for k in range(len(path)-1):
                    city1 = cities[path[k]]
                    city2 = cities[path[k+1]]
                    distance = haversine(city1[0], city1[1], city2[0], city2[1])
                    path_distance += distance
                path_distances.append(path_distance)
            pathDistance[i][j] = path_distances
            pathDistance[j][i] = path_distances[::-1]


    for i in range(numberOfCities):
        for j in range(numberOfCities):
            # the upper triangle of the array has to be equal to the lower triangle beacuse the distance between two city must be same
            if i < j:
                pathDistance[j][i] = pathDistance[i][j]
            
                
    for i in range(numberOfCities):
        for j in range(numberOfCities):
            
            if i == j:
                pathDistance[j][i] = np.zeros(1)

    
    for i in range(len(cities_with_ids)):
        for j in range(len(cities_with_ids)):
            if i != j:
                city1 = cities_with_ids[i]
                city2 = cities_with_ids[j]
                city_id1 = city1['id']
                city_id2 = city2['id']
                pathDistance_values = pathDistance[i][j]
                for pathDistances in pathDistance_values:
                    
                    pathDistances = round(pathDistances, 2)
                    city_pair_pathDistance = [(city_id1, city_id2), pathDistances]
                    city_pairs_pathDistance.append(city_pair_pathDistance)

    # Print the toll values for each city pair
    for (city_id1, city_id2), pathDistances in city_pairs_pathDistance:
        pathDistance_list.append([city_id1, city_id2, pathDistances])  
        
    '''              
                    
    
    fuel = np.around(np.random.triangular(fuelPricesDistributions_lower, fuelPricesDistributions_mode, fuelPricesDistributions_upper, size=numberOfCities), decimals=2)

    for i, city in enumerate(cities_with_ids):
        city_id = city['id']
        fuel_price = fuel[i]
        fuel_prices_array.append([city_id, fuel_price])

    
    for city in cities_with_ids:
        city_id = city['id']
        for city_fuel_price in fuel_prices_array:
            if city_fuel_price[0] == city_id:
                fuelPrice.append(city_fuel_price)
                


    #fuel consumption for the truck 
    fuelConsumptions = np.array([[np.around(np.random.triangular(fuelConsumptionDistribution_lower,fuelConsumptionDistribution_mode, fuelConsumptionDistribution_upper), decimals=2)for j in range(numberOfCities)] for k in range(numberOfCities)])
    
    # Create the city_pairs_fuel list
    for i in range(len(cities_with_ids)):
        for j in range(len(cities_with_ids)):
            if i != j:
                city1 = cities_with_ids[i]
                city2 = cities_with_ids[j]
                city_id1 = city1['id']
                city_id2 = city2['id']
                fuelC_values = fuelConsumptions[i][j].tolist()
                city_pair_fuel = [(city_id1, city_id2), fuelC_values]
                city_pairs_fuelC.append(city_pair_fuel)


    for (city_id1, city_id2), fuelC_values in city_pairs_fuelC:
        fuelConsumption.append([city_id1, city_id2, fuelC_values])

    #average velocity from city i to j
    velocity = np.array([[np.around(np.random.triangular(averageVelocityDistribution_lower,averageVelocityDistribution_mode, averageVelocityDistribution_upper), decimals=0)for j in range(numberOfCities)] for k in range(numberOfCities)])

    for i in range(numberOfCities):
        for j in range(numberOfCities):
            # the upper triangle of the array has to be equal to the lower triangle beacuse the distance between two city must be same
            if i < j:
                velocity[j][i] = velocity[i][j]
            
                
    for i in range(len(cities_with_ids)):
        for j in range(len(cities_with_ids)):
            if i != j:
                city1 = cities_with_ids[i]
                city2 = cities_with_ids[j]
                city_id1 = city1['id']
                city_id2 = city2['id']
                velocity_values = velocity[i][j].tolist()
                city_pair_vel = [(city_id1, city_id2), velocity_values]
                city_pairs_velocity.append(city_pair_vel)


    for (city_id1, city_id2), velocity_values in city_pairs_velocity:
        velocity_list.append([city_id1, city_id2, velocity_values])

    
    #time spent on travelling city i to city j (distances/average velocity)
    for i in range(len(distances_with_city_ids)):
        time.append([distances_with_city_ids[i][0], distances_with_city_ids[i][1], round(distances_with_city_ids[i][2]/velocity_list[i][2],2)])
    
  
    # Print the toll values for each city pair
    for i in range(len(toll)):
        total_costs.append([time[i][0], time[i][1] , round((time[i][2]* alpha) + toll[i][2],2)])

      
    
    '''
    for i in range(len(time)):
        city1 = time[i][0]
        city2 = time[i][1]
        distance = distance(city1, city2) 
        fuel_consumption = fuel_consumption(city1, city2) 
        fuel_price = fuel_price(city1) 
        fuel_cost = distance * fuel_consumption * fuel_price
        total_cost = round((time[i][2] * alpha) + toll[i][2] + fuel_cost, 2)
        total_costs.append([city1, city2, total_cost])
    '''

    
# take current time
now = datetime.now()

# save all csv file with one function call
def saveAllAsCSV():
    saveCitiesOutput()
    saveDistancesOutput()
    saveVelocityOutput()
    saveTimeOutput()
    saveTollOutput()
    saveCostOutput()
    saveFuelPriceOutput()
    saveFuelConsumptionOutput()
    #savePathDistanceOutput()
    #savePathsOutput()
    


directory = os.path.join("Outputs" ,  now.strftime("%Y-%m-%d_%H-%M-%S"))
sol_directory = os.path.join("Solutions" ,  now.strftime("%Y-%m-%d_%H-%M-%S"))
sol_nn_directory = os.path.join("Solutions_NN" ,  now.strftime("%Y-%m-%d_%H-%M-%S"))
sol_opt_directory = os.path.join("Solutions_Exact" ,  now.strftime("%Y-%m-%d_%H-%M-%S"))


def saveTimeOutput():
    df_timeMatrix = pd.DataFrame(time)
    # define csv file name
    time_file_name = 'time_' + now.strftime("%Y-%m-%d_%H-%M-%S") + '.csv'
        # set column names
    df_timeMatrix.columns = ['from', 'to', 'time']
    
    if not os.path.exists(directory):
        os.makedirs(directory)
    # writing data frame to a CSV file
    df_timeMatrix.to_csv(os.path.join(directory, time_file_name), sep=';', index=False)

'''
def savePathDistanceOutput():
    df_distanceMatrix = pd.DataFrame(pathDistance_list)
    # define csv file name
    distance_file_name = 'pathDistance_' + now.strftime("%Y-%m-%d_%H-%M-%S") + '.csv'
        # set column names
    df_distanceMatrix.columns = ['from', 'to', 'pathDistance']
    # writing data frame to a CSV file
    if not os.path.exists(directory):
        os.makedirs(directory)
    # writing data frame to a CSV file
    df_distanceMatrix.to_csv(os.path.join(directory, distance_file_name), sep=';', index=False)
'''

def saveFuelPriceOutput():
    df_fuelPriceMatrix = pd.DataFrame(fuelPrice)
    # define csv file name
    fuel_file_name = 'fuel_price_' + now.strftime("%Y-%m-%d_%H-%M-%S") + '.csv'
    # set column names
    df_fuelPriceMatrix.columns = ['cityID', 'fuelPrice']
    if not os.path.exists(directory):
        os.makedirs(directory)
    # writing data frame to a CSV file
    df_fuelPriceMatrix.to_csv(os.path.join(directory, fuel_file_name), sep=';', index=False)

def saveVelocityOutput():
    df_velocityMatrix = pd.DataFrame(velocity_list)
    # define csv file name
    velocity_file_name = 'velocity_' + now.strftime("%Y-%m-%d_%H-%M-%S") + '.csv'
        # set column names
    df_velocityMatrix.columns = ['from', 'to', 'velocity']
    # writing data frame to a CSV file
    if not os.path.exists(directory):
        os.makedirs(directory)
    # writing data frame to a CSV file
    df_velocityMatrix.to_csv(os.path.join(directory, velocity_file_name), sep=';', index=False)

def saveFuelConsumptionOutput():
    df_FuelConsumptionMatrix = pd.DataFrame(fuelConsumption)
    # define csv file name
    fuel_consumption_file_name = 'fuel_consumption_' + now.strftime("%Y-%m-%d_%H-%M-%S") + '.csv'
        # set column names
    df_FuelConsumptionMatrix.columns = ['from', 'to', 'FuelConsumption']
    # writing data frame to a CSV file
    if not os.path.exists(directory):
        os.makedirs(directory)
    # writing data frame to a CSV file
    df_FuelConsumptionMatrix.to_csv(os.path.join(directory, fuel_consumption_file_name), sep=';', index=False)

def saveTollOutput():
    df_TollPriceMatrix = pd.DataFrame(toll)
    # define csv file name
    toll_file_name = 'toll_price_' + now.strftime("%Y-%m-%d_%H-%M-%S") + '.csv'
    # set column names
    df_TollPriceMatrix.columns = ['from', 'to','TollPrice']
    # writing data frame to a CSV file
    if not os.path.exists(directory):
        os.makedirs(directory)
    # writing data frame to a CSV file
    df_TollPriceMatrix.to_csv(os.path.join(directory, toll_file_name), sep=';', index=False)

'''
def savePathsOutput():
    df_paths = pd.DataFrame(path_list)
    # define csv file name
    paths_file_name = 'paths_' + now.strftime("%Y-%m-%d_%H-%M-%S") + '.csv'
        # set column names
    df_paths.columns = ['from', 'to', 'paths']
    # writing data frame to a CSV file
    if not os.path.exists(directory):
        os.makedirs(directory)
    # writing data frame to a CSV file
    df_paths.to_csv(os.path.join(directory, paths_file_name), sep=';', index=False)
'''

def saveDistancesOutput():
    df_DistancesMatrix = pd.DataFrame(distances_with_city_ids)
    # define csv file name
    distances_file_name = 'distances_' + now.strftime("%Y-%m-%d_%H-%M-%S") + '.csv'
    # set column names
    df_DistancesMatrix.columns = ['from', 'to', 'distances']
    # writing data frame to a CSV file
    if not os.path.exists(directory):
        os.makedirs(directory)
    # writing data frame to a CSV file
    df_DistancesMatrix.to_csv(os.path.join(directory, distances_file_name), sep=';', index=False)

def saveCitiesOutput():
    df_CitiesMatrix = pd.DataFrame(cities_with_ids)
    # define csv file name
    cities_file_name = 'cities_' + now.strftime("%Y-%m-%d_%H-%M-%S") + '.csv'
    # set column names
    df_CitiesMatrix.columns = ['City ID','latitude', 'longitude']
    # writing data frame to a CSV file
    if not os.path.exists(directory):
        os.makedirs(directory)
    # writing data frame to a CSV file
    df_CitiesMatrix.to_csv(os.path.join(directory, cities_file_name), sep=';', index=False)

def saveCostOutput():
    df_totalCostMatrix = pd.DataFrame(total_costs)
    # define csv file name
    totalCost_file_name = 'totalCost_' + now.strftime("%Y-%m-%d_%H-%M-%S") + '.csv'
    # set column names
    df_totalCostMatrix.columns = ['from', 'to', 'Total Cost']
    # writing data frame to a CSV file
    if not os.path.exists(directory):
        os.makedirs(directory)
    # writing data frame to a CSV file
    df_totalCostMatrix.to_csv(os.path.join(directory, totalCost_file_name), sep=';', index=False)
