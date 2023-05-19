from random_instance import *
import ast


def wanted_date():
    options = os.listdir(os.path.join("Outputs", ))
    print(options)

    date = input('date: ')
    directory_read = os.path.join("Outputs", date)
    return date, directory_read


all_time = []
all_fuelPrice = []
all_velocity = []
all_FuelConsumption = []
all_TollPrice = []
all_paths = []
all_distances = []
all_cities = []
all_totalCost = []
fuel_needed = [] 


def readtimeOutput(date, directory_read):
    time_file_name = 'time_' + date + '.csv'

    # construct the full file path using the 'directory' variable
    file_path = os.path.join(directory_read, time_file_name)

    # read the CSV file into a DataFrame
    df = pd.read_csv(file_path, sep=';')
    from_read = df['from'].tolist()
    to_read = df['to'].tolist()
    # extract the 'time' column as a list
    time_read = df['time'].tolist()

    for i in range(len(time_read)):
        all_time.append([from_read[i], to_read[i], time_read[i]])



def readfuelPriceOutput(date, directory_read):
    fuelPrice_file_name = 'fuel_price_' + date + '.csv'

    # construct the full file path using the 'directory' variable
    file_path = os.path.join(directory_read, fuelPrice_file_name)

    # read the CSV file into a DataFrame
    df = pd.read_csv(file_path, sep=';')
    City_ID_read = df['cityID'].tolist()
    # extract the 'time' column as a list
    fuelPrice_read = df['fuelPrice'].tolist()

    for i in range(len(fuelPrice_read)):
        all_fuelPrice.append([City_ID_read[i], fuelPrice_read[i]])


def readvelocityOutput(date, directory_read):
    velocity_file_name = 'velocity_' + date + '.csv'

    # construct the full file path using the 'directory' variable
    file_path = os.path.join(directory_read, velocity_file_name)

    # read the CSV file into a DataFrame
    df = pd.read_csv(file_path, sep=';')
    from_read = df['from'].tolist()
    to_read = df['to'].tolist()
    # extract the 'time' column as a list
    velocity_read = df['velocity'].tolist()

    for i in range(len(velocity_read)):
        all_velocity.append([(from_read[i], to_read[i]), velocity_read[i]])


def readfuelConsumptionOutput(date, directory_read):
    FuelConsumption_file_name = 'fuel_consumption_' + date + '.csv'

    # construct the full file path using the 'directory' variable
    file_path = os.path.join(directory_read, FuelConsumption_file_name)

    # read the CSV file into a DataFrame
    df = pd.read_csv(file_path, sep=';')
    from_read = df['from'].tolist()
    to_read = df['to'].tolist()
    # extract the 'time' column as a list
    FuelConsumption_read = df['FuelConsumption'].tolist()

    for i in range(len(FuelConsumption_read)):
        all_FuelConsumption.append([from_read[i], to_read[i], FuelConsumption_read[i]])


def readtollPriceOutput(date, directory_read):
    TollPrice_file_name = 'toll_price_' + date + '.csv'

    # construct the full file path using the 'directory' variable
    file_path = os.path.join(directory_read, TollPrice_file_name)

    # read the CSV file into a DataFrame
    df = pd.read_csv(file_path, sep=';')
    from_read = df['from'].tolist()
    to_read = df['to'].tolist()
    # extract the 'time' column as a list
    TollPrice_read = df['TollPrice'].tolist()

    for i in range(len(TollPrice_read)):
        all_TollPrice.append([from_read[i], to_read[i], TollPrice_read[i]])



def readdistancesOutput(date, directory_read):
    distances_file_name = 'distances_' + date + '.csv'

    # construct the full file path using the 'directory' variable
    file_path = os.path.join(directory_read, distances_file_name)

    # read the CSV file into a DataFrame
    df = pd.read_csv(file_path, sep=';')
    from_read = df['from'].tolist()
    to_read = df['to'].tolist()
    # extract the 'time' column as a list
    distances_read = df['distances'].tolist()

    for i in range(len(distances_read)):
        all_distances.append([from_read[i], to_read[i], distances_read[i]])


def readcitiesOutput(date, directory_read):
    cities_file_name = 'cities_' + date + '.csv'

    # construct the full file path using the 'directory' variable
    file_path = os.path.join(directory_read, cities_file_name)

    # read the CSV file into a DataFrame
    df = pd.read_csv(file_path, sep=';')
    cityId_read = df['City ID'].tolist()
    x_read = df['latitude'].tolist()
    # extract the 'time' column as a list
    y_read = df['longitude'].tolist()

    for i in range(len(cityId_read)):
        all_cities.append([cityId_read[i], x_read[i], y_read[i]])

def readtotalCostOutput(date, directory_read):
    totalCost_file_name = 'totalCost_' + date + '.csv'

    # construct the full file path using the 'directory' variable
    file_path = os.path.join(directory_read, totalCost_file_name)

    # read the CSV file into a DataFrame
    df = pd.read_csv(file_path, sep=';')
    from_read = df['from'].tolist()
    to_read = df['to'].tolist()
    # extract the 'totalCost' column as a list
    totalCost_read = df['Total Cost'].tolist()

    for i in range(len(totalCost_read)):
        all_totalCost.append([from_read[i], to_read[i], totalCost_read[i]])

def finding_fuel_needed():
    for i in range(len(all_distances)):
        fuel_needed.append([all_distances[i][0] ,all_distances[i][1] , all_distances[i][2] * all_FuelConsumption[i][2]])

def readAllAsCSV():
    date, directory_read = wanted_date()
    readcitiesOutput(date, directory_read)
    readdistancesOutput(date, directory_read)
    readvelocityOutput(date, directory_read)
    readtimeOutput(date, directory_read)
    readtollPriceOutput(date, directory_read)
    readtotalCostOutput(date, directory_read)
    readfuelPriceOutput(date, directory_read)
    readfuelConsumptionOutput(date, directory_read)
    finding_fuel_needed()



def dijkstra_path(start_city_id, end_city_id, costs, distances):
    # Create a dictionary to store the minimum costs from start_city to each city
    min_costs = {start_city_id: 0}
    
    # Create a dictionary to store the previous city in the optimal path
    previous_city = {}
    
    # Create a set to keep track of visited cities
    visited = set()
    
    # Start with the start_city as the current city
    current_city = start_city_id
    
    while current_city != end_city_id:
        # Mark the current city as visited
        visited.add(current_city)
        
        # Get the cost to reach the current city
        current_cost = min_costs[current_city]
        
        # Find the neighbors of the current city
        neighbors = [c2 for (c1, c2, _) in costs if c1 == current_city]
        
        for neighbor in neighbors:
            # Check if the distance is -1, skip this city
            if any(c1 == current_city and c2 == neighbor and d == -1 for (c1, c2, d) in distances):
                continue
            
            # Calculate the cost to reach the neighbor from the start_city
            neighbor_cost = current_cost + next(c for (c1, c2, c) in costs if c1 == current_city and c2 == neighbor)
            
            if neighbor not in min_costs or neighbor_cost < min_costs[neighbor]:
                # Update the minimum cost and previous city for the neighbor
                min_costs[neighbor] = neighbor_cost
                previous_city[neighbor] = current_city
        
        # Find the next unvisited city with the minimum cost
        unvisited_cities = {city: min_costs[city] for city in min_costs if city not in visited}
        
        if not unvisited_cities:
            # No unvisited cities remaining, end_city_id is not reachable
            raise ValueError("End city is not reachable from the start city.")
        
        current_city = min(unvisited_cities, key=unvisited_cities.get)
    
    # Reconstruct the optimal path
    path = []
    city = end_city_id
    while city != start_city_id:
        path.append(city)
        city = previous_city[city]
    path.append(start_city_id)
    
    # Reverse the path since we built it backwards
    path.reverse()
    
    return path


def dijkstra_distance(start_city_id, end_city_id, distances):
    # Create a dictionary to store the minimum distances from start_city to each city
    min_distance = {start_city_id: 0}
    
    # Create a dictionary to store the previous city in the optimal path
    previous_city = {}
    
    # Create a set to keep track of visited cities
    visited = set()
    
    # Start with the start_city as the current city
    current_city = start_city_id
    
    while current_city != end_city_id:
        # Mark the current city as visited
        visited.add(current_city)
        
        # Get the cost to reach the current city
        current_distance = min_distance[current_city]
        
        # Find the neighbors of the current city
        neighbors = [c2 for (c1, c2, _) in distances if c1 == current_city]
        
        for neighbor in neighbors:
            # Check if the distance is -1, skip this city
            if any(c1 == current_city and c2 == neighbor and d == -1 for (c1, c2, d) in distances):
                continue
            
            # Calculate the distance to reach the neighbor from the start_city
            neighbor_distance = current_distance + next(c for (c1, c2, c) in distances if c1 == current_city and c2 == neighbor)
            
            if neighbor not in min_distance or neighbor_distance < min_distance[neighbor]:
                # Update the minimum distance and previous city for the neighbor
                min_distance[neighbor] = neighbor_distance
                previous_city[neighbor] = current_city
        
        # Find the next unvisited city with the minimum distance
        unvisited_cities = {city: min_distance[city] for city in min_distance if city not in visited}
        
        if not unvisited_cities:
            # No unvisited cities remaining, end_city_id is not reachable
            raise ValueError("End city is not reachable from the start city.")
        
        current_city = min(unvisited_cities, key=unvisited_cities.get)
    
    # Reconstruct the optimal path
    path = []
    city = end_city_id
    while city != start_city_id:
        path.append(city)
        city = previous_city[city]
    path.append(start_city_id)
    
    # Reverse the path since we built it backwards
    path.reverse()
    
    return path

"""
    Determine where to buy fuel along the given path and calculate the total fuel cost.
    Args:
    path (list): A list of cities in the order they are visited.
    distance (list): A list of distance values for each leg of the journey.
    fuel_consumption (list): A list of fuel consumption values for each leg of the journey.
    fuel_price (list): A list of fuel prices for each city.
    max_fuel_capacity (float): The maximum fuel capacity of the vehicle.
    starting_fuel (float): The starting fuel level in the vehicle.
    Returns:
    A tuple containing a list of tuples (city, fuel), where each tuple represents a fuel purchase,
    and the total fuel cost for the entire journey.
    
"""

def minimize_fuel_cost(path, fuel_needed, fuel_price, max_fuel_capacity, starting_fuel):
    fuel_cost = 0  # Initialize total fuel cost
    current_fuel = starting_fuel  # Initialize current fuel level
    fuel_purchases = []  # List to store fuel purchases
    
    fuel_prices_in_path = {}  # Dictionary to store fuel prices of cities in the path
    
    # Find the fuel prices of cities in the path
    for city, price in fuel_price:
        if city in path:
            fuel_prices_in_path[city] = price
    
    # Create a copy of fuel_prices_in_path
    modified_fuel_prices = fuel_prices_in_path.copy()

    # Set the last element's value to infinity
    last_key = path[-1]
    modified_fuel_prices[last_key] = float('inf')

    # Sort the modified dictionary by values
    sorted_fuel_prices = sorted(modified_fuel_prices.items(), key=lambda x: x[1])
    

       
    # Iterate through each city in the path
    for i in range(len(path)-1):
        source_city = path[i]
        destination_city = path[i+1]
        fuel_needed_current_leg = None
        
        # Find the fuel needed for the current leg of the journey
        for fuel_leg in fuel_needed:
            if fuel_leg[0] == source_city and fuel_leg[1] == destination_city:
                fuel_needed_current_leg = fuel_leg[2]
                break
        
        # Check if the destination is the last city
        if i ==  path[len(path) - 2] and i == path[0]:
            
            # Finish the journey at the destination city (no need to buy fuel)
            if fuel_needed_current_leg <= current_fuel:
                current_fuel = current_fuel - fuel_needed_current_leg
                continue
            else:
                # Fill the tank with the amount of fuel needed for the current leg
                fuel_to_buy = fuel_needed_current_leg - current_fuel
                fuel_cost += fuel_to_buy * fuel_prices_in_path[source_city]
                current_fuel = 0
                fuel_purchases.append((source_city, fuel_to_buy))
                break
                
                
                
                
        # Check if the destination is the minimum fuel price city
        if source_city == sorted_fuel_prices[0][0]:
            # Fill the tank with the amount of fuel on hand to reach the next city
            fuel_to_buy = max_fuel_capacity - current_fuel
            fuel_cost += fuel_to_buy * fuel_prices_in_path[source_city]
            current_fuel = max_fuel_capacity - fuel_needed_current_leg
            fuel_purchases.append((source_city, fuel_to_buy))
            
        else:
            # Check if it's possible to reach the destination with current fuel
            if fuel_needed_current_leg <= current_fuel:
                current_fuel = current_fuel - fuel_needed_current_leg
                # Find the next minimum fuel price city that is reachable with current fuel
                continue       
            else:
                # Fill the tank with the amount of fuel needed for the current leg
                fuel_to_buy = fuel_needed_current_leg - current_fuel
                fuel_cost += fuel_to_buy * fuel_prices_in_path[source_city]
                current_fuel = 0
                fuel_purchases.append((source_city, round(fuel_to_buy, 2)))
    
    # Return the list of fuel purchases and the total fuel cost
    return fuel_purchases, fuel_cost




def calculate_path_distance(path, distances):
    total_distance = 0
    
    for i in range(len(path)-1):
        city_id1 = path[i]
        city_id2 = path[i+1]
        
        # Look up the distance between the two cities
        for (city_1,city_2, distance) in distances:
            if (city_id1, city_id2) == (city_1,city_2) or (city_id2, city_id1) == (city_1,city_2):
                total_distance += distance
                break
                
    return total_distance

def calculate_path_time(path, times):
  
    total_time = 0
    
    for i in range(len(path)-1):
        city_id1 = path[i]
        city_id2 = path[i+1]
        
        # Look up the time between the two cities
        for (city_1,city_2, time) in times:
            if (city_id1, city_id2) == (city_1,city_2) or (city_id2, city_id1) == (city_1,city_2):
                total_time += time
                break
                
    return total_time


def calculate_path_toll(path, tolls):
    
    total_toll = 0
    
    for i in range(len(path)-1):
        city_id1 = path[i]
        city_id2 = path[i+1]
        
        # Look up the toll between the two cities
        for (city_1,city_2, toll) in tolls:
            if (city_id1, city_id2) == (city_1,city_2) or (city_id2, city_id1) == (city_1,city_2):
                total_toll += toll
                break
                
    return total_toll

def calculate_path_total_cost(path, costs):
    
        total_total_cost = 0
        
        for i in range(len(path)-1):
            city_id1 = path[i]
            city_id2 = path[i+1]
            
            # Look up the total_cost between the two cities
            for (city_1,city_2, total_costs) in costs:
                if (city_id1, city_id2) == (city_1,city_2) or (city_id2, city_id1) == (city_1,city_2):
                    total_total_cost += total_costs
                    break
                    
        return total_total_cost
    

def sol_distances_with_ids(i, j):
    start_city_id = int(i)
    end_city_id = int(j)
    path = dijkstra_distance(start_city_id, end_city_id, all_distances)
    total_distance = calculate_path_distance(path, all_distances) 
    total_time = calculate_path_time(path, all_time) 
    total_toll = calculate_path_toll(path, all_TollPrice)
    total_cost = calculate_path_total_cost(path, all_totalCost)
    print(f"Path from {start_city_id} to {end_city_id}: {path}")
    print(f"Total Distance: {total_distance}")
    print(f'Total Time:  {total_time}')
    print(f'Total Toll Price: {total_toll}')
    print(f'Total Cost: {total_cost}')

    return path, total_distance, total_time, total_toll, total_cost

def sol_for_distance():
    IDs = []
    for id in range(len(all_cities)):
        IDs.append(all_cities[id][0])
    for i in IDs:
        for j in IDs:
            if i != j:
                start_city_id = int(i)
                end_city_id = int(j)
                path = dijkstra_distance(start_city_id, end_city_id, all_distances)
                total_distance = calculate_path_distance(path, all_distances) 
                total_time = calculate_path_time(path, all_time) 
                total_toll = calculate_path_toll(path, all_TollPrice)
                total_cost = calculate_path_total_cost(path, all_totalCost)
                print(f"Path from {start_city_id} to {end_city_id}: {path}")
                print(f"Total Distance: {total_distance}")
                print(f'Total Time:  {total_time}')
                print(f'Total Toll Price: {total_toll}')
                print(f'Total Cost: {total_cost}')

all_distance_output = []
def save_for_distance():

    IDs = []
    for id in range(len(all_cities)):
        IDs.append(all_cities[id][0])
    
    for i in IDs:
        for j in IDs:
            if i != j:
                start_city_id = int(i)
                end_city_id = int(j)
                path = dijkstra_distance(start_city_id, end_city_id, all_distances)
                total_distance = calculate_path_distance(path, all_distances) 
                total_time = calculate_path_time(path, all_time) 
                total_toll = calculate_path_toll(path, all_TollPrice)
                total_cost = calculate_path_total_cost(path, all_totalCost)
                fuel_purchase , fuel_cost = minimize_fuel_cost(path, fuel_needed, all_fuelPrice, F, Sf)
                all_distance_output.append([start_city_id, end_city_id, total_cost + fuel_cost, total_distance, path, total_toll, total_time, fuel_purchase, fuel_cost])
                
def saveDistanceSolutionOutput():
    df_SolutionDistMatrix = pd.DataFrame(all_distance_output)
    # define csv file name
    solD_file_name = 'solution_distance_' + now.strftime("%Y-%m-%d_%H-%M-%S") + '.csv'
    # set column names
    df_SolutionDistMatrix.columns = ['from', 'to', 'Total Cost', 'pathDistance', 'path', 'toll', 'time', 'City_Amount', 'fuel cost']

    if not os.path.exists(sol_directory):
        os.makedirs(sol_directory)
    # writing data frame to a CSV file
    df_SolutionDistMatrix.to_csv(os.path.join(sol_directory, solD_file_name), sep=';', index=False)
   

def sol_cost_with_ids(i,j):
    
    start_city_id = int(i)
    end_city_id = int(j)
    path = dijkstra_path(start_city_id, end_city_id, all_totalCost, all_distances)
    total_distance = calculate_path_distance(path, all_distances) 
    total_time = calculate_path_time(path, all_time) 
    total_toll = calculate_path_toll(path, all_TollPrice)
    total_cost = calculate_path_total_cost(path, all_totalCost)
    print(f"Path from {start_city_id} to {end_city_id}: {path}")
    print(f"Total Distance: {total_distance}")
    print(f'Total Time:  {total_time}')
    
    return path, total_distance, total_time, total_toll, total_cost


def sol_for_cost():
    IDs = []
    for id in range(len(all_cities)):
        IDs.append(all_cities[id][0])
    
    for i in IDs:
        for j in IDs:
            if i != j:
                start_city_id = int(i)
                end_city_id = int(j)
                path = dijkstra_path(start_city_id, end_city_id, all_totalCost, all_distances)
                total_distance = calculate_path_distance(path, all_distances) 
                total_time = calculate_path_time(path, all_time) 
                total_cost = calculate_path_total_cost(path, all_totalCost)
                fuel_purchase , fuel_cost = minimize_fuel_cost(path, fuel_needed, all_fuelPrice, F, Sf)
                print(f"Path from {start_city_id} to {end_city_id}: {path}")
                print(f"Total Distance: {total_distance}")
                print(f'Total Time:  {total_time}')
                print(f'Total Cost: {total_cost + fuel_cost}')   
                print(f'Fuel buy in: {fuel_purchase}')
                print(f'Fuel Cost: {fuel_cost}')

all_cost_output = []
def save_for_cost():
    
    IDs = []
    for id in range(len(all_cities)):
        IDs.append(all_cities[id][0])
        
    for i in IDs:
        for j in IDs:
            if i != j:
                start_city_id = int(i)
                end_city_id = int(j)
                path = dijkstra_path(start_city_id, end_city_id, all_totalCost, all_distances)
                total_distance = calculate_path_distance(path, all_distances) 
                total_time = calculate_path_time(path, all_time) 
                total_toll = calculate_path_toll(path, all_TollPrice)
                total_cost = calculate_path_total_cost(path, all_totalCost)
                fuel_purchase , fuel_cost = minimize_fuel_cost(path, fuel_needed, all_fuelPrice, F, Sf)

                all_cost_output.append([start_city_id, end_city_id, total_cost + fuel_cost, total_distance, path, total_toll, total_time, fuel_purchase, fuel_cost])
                                
def saveCostSolutionOutput():
    df_SolutionMatrix = pd.DataFrame(all_cost_output)
    # define csv file name
    sol_file_name = 'solution_cost_' + now.strftime("%Y-%m-%d_%H-%M-%S") + '.csv'
    # set column names
    df_SolutionMatrix.columns = ['from', 'to', 'Total Cost', 'pathDistance', 'path', 'toll', 'time', 'City-Amount', 'fuel cost']

    
    if not os.path.exists(sol_directory):
        os.makedirs(sol_directory)
    # writing data frame to a CSV file
    df_SolutionMatrix.to_csv(os.path.join(sol_directory, sol_file_name), sep=';', index=False)

          
def saveSolution():
    readAllAsCSV()
    save_for_cost()
    saveCostSolutionOutput()

def saveDistanceSolution():
    readAllAsCSV()
    save_for_distance()
    saveDistanceSolutionOutput()
    
#readAllAsCSV()
#sol_cost_with_ids(3,7)
#sol_for_cost()

#saveDistanceSolution()
