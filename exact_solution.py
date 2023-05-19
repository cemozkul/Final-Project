import os
import sys

import gurobipy as gp
import pandas as pd
from gurobipy import GRB
import time as time_calculator
from random_instance import *


def wanted_date():
    options = os.listdir(os.path.join("Outputs", ))
    print(options)

    date = input('date: ')
    directory_read = os.path.join("Outputs", date)
    return date, directory_read


exact_time = []
exact_fuelPrice = []
exact_velocity = []
exact_FuelConsumption = []
exact_TollPrice = []
exact_paths = []
exact_distances = []
exact_cities = []
exact_totalCost = []


def readTimeOutput(date, directory_read):
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
        exact_time.append([from_read[i], to_read[i], time_read[i]])


def readFuelPriceOutput(date, directory_read):
    fuelPrice_file_name = 'fuel_price_' + date + '.csv'

    # construct the full file path using the 'directory' variable
    file_path = os.path.join(directory_read, fuelPrice_file_name)

    # read the CSV file into a DataFrame
    df = pd.read_csv(file_path, sep=';')
    City_ID_read = df['cityID'].tolist()
    # extract the 'time' column as a list
    fuelPrice_read = df['fuelPrice'].tolist()

    for i in range(len(fuelPrice_read)):
        exact_fuelPrice.append([City_ID_read[i], fuelPrice_read[i]])


def readVelocityOutput(date, directory_read):
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
        exact_velocity.append([(from_read[i], to_read[i]), velocity_read[i]])


def readFuelConsumptionOutput(date, directory_read):
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
        exact_FuelConsumption.append([from_read[i], to_read[i], FuelConsumption_read[i]])


def readTollPriceOutput(date, directory_read):
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
        exact_TollPrice.append([from_read[i], to_read[i], TollPrice_read[i]])


def readDistancesOutput(date, directory_read):
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
        exact_distances.append([from_read[i], to_read[i], distances_read[i]])


def readCitiesOutput(date, directory_read):
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
        exact_cities.append([cityId_read[i], x_read[i], y_read[i]])


def readTotalCostOutput(date, directory_read):
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
        exact_totalCost.append([from_read[i], to_read[i], totalCost_read[i]])


def readexactAsCSV():
    date, directory_read = wanted_date()
    readCitiesOutput(date, directory_read)
    readDistancesOutput(date, directory_read)
    readVelocityOutput(date, directory_read)
    readTimeOutput(date, directory_read)
    readTollPriceOutput(date, directory_read)
    readTotalCostOutput(date, directory_read)
    readFuelPriceOutput(date, directory_read)
    readFuelConsumptionOutput(date, directory_read)


#readexactAsCSV()



def exact_solve(source, destination):
    start_time = time_calculator.time()
    nCities = len(exact_cities)
    distances_dict = {(distance[0], distance[1]): distance[2] for distance in exact_distances}
    time_dict = {(time[0], time[1]): time[2] for time in exact_time}
    toll_dict = {(toll[0], toll[1]): toll[2] for toll in exact_TollPrice}
    fuel_dict = {fuel[0]: fuel[1] for fuel in exact_fuelPrice}
    consumption_dict = {(consumption[0], consumption[1]): consumption[2] for consumption in exact_FuelConsumption}

    
    try:

        # Create a new model
        model = gp.Model("InternationalTruckRouting")
        # Create decision variables

        # Initialize vehicle routing decision variables x[i][j] == 1 if vehicle travel city i to city j.
        # define list of tuples for valid index pairs
        pairs = [(i, j) for i in range(1, nCities + 1) for j in range(1, nCities + 1) if i != j]

        # add variables for valid index pairs
        x = model.addVars(pairs, vtype=GRB.BINARY, name='x')

        # amountOfFuelBought[i] is a continues variable which define the amount of fuel bought in city i
        amountOfFuelBought = [model.addVar(0, F, 0, GRB.CONTINUOUS, f"Fuel Bought Amount ")
                              for i in range(1, nCities + 1)]

        # availableFuelAmount[i] is a continuous variable which defines the amount of fuel bought in city i

        availableFuelAmount = [model.addVar(0, F, 0, GRB.CONTINUOUS, f"Available Fuel Amount ")
                               for i in range(1, nCities + 1)]

        # Objective function: minimize cost total travelling cost from source to destination
        model.setObjective(gp.quicksum(
            ((alpha * time_dict[(i, j)] + toll_dict[(i, j)]) * x[i, j]) for i, j in
            [(i, j) for i in range(1, nCities + 1) for j in range(1, nCities + 1) if i != j])
                           + gp.quicksum((fuel_dict[i] * amountOfFuelBought[i - 1]) for i in range(1, nCities + 1)),
                           GRB.MINIMIZE)

        # Creating constraints of Vehicle Routing Problem

        sourceID = source  # this will change

        # Vehicle should leave from source
        # model.addConstr(sum(x[sourceID, i] for i in range(nCities) not in {sourceID}), gp.GRB.EQUAL, 1,
        #                "Vehicle should leave from the source")
        model.addConstr(gp.quicksum(x[sourceID, i] for i in range(1, nCities + 1) if sourceID != i) == 1,
                        "Vehicle should leave from the source")

        destinationID = destination  # this will change

        # Vehicle should enter the destination
        # model.addConstr(sum(x[destinationID, i] for i in range(nCities) not in {destinationID}), gp.GRB.EQUAL, 1,
        #               "Vehicle should enter the destination")
        model.addConstr(gp.quicksum(x[i, destinationID] for i in range(1, nCities + 1) if destinationID != i) == 1,
                        "Vehicle should enter the destination")

        # the cities are not the source and destination if the vehicle enters the city, it must leave from that city
        for i in range(1, nCities + 1):
            if i not in {sourceID, destinationID}:
                model.addConstr(
                    gp.quicksum(x[j, i] for j in range(1, nCities + 1) if j != i) - gp.quicksum(
                        x[i, j] for j in range(1, nCities + 1) if j != i) == 0)

        # sub-tour elimination constraint
        for i in range(1, nCities + 1):
            for j in range(1, nCities + 1):
                if i != j:
                    model.addConstr(x[i, j] + x[j, i] <= 1, name="sub-tour elimination")

        # If going city i to j is labelled -1 you shouldn't go to this city
        for i in range(1, nCities + 1):
            for j in range(1, nCities + 1):
                if i != j:
                    if distances_dict[(i, j)] == -1:
                        model.addConstr(x[i, j] == 0)

        # initial fuel capacity
        model.addConstr(availableFuelAmount[sourceID - 1] == Sf, name="initial fuel")

        # Constrain guarantees the fuel amount should be greater than or equal the city that will be gone.
        for i in range(1, nCities + 1):
            for j in range(1, nCities + 1):
                if i != j and (i, j) in distances_dict and (i, j) in consumption_dict:
                    model.addConstr(availableFuelAmount[j - 1] <= (F * (1 - x[i, j]))
                                    + amountOfFuelBought[i - 1] + availableFuelAmount[i - 1]
                                    - (consumption_dict[(i, j)] * distances_dict[(i, j)]) * x[i, j],
                                    name=f"fuel conservation")

        # The fuel amount that we can buy from each city cannot exceed the fuel capacity of the vehicle
        model.addConstrs((F - availableFuelAmount[i - 1] >= amountOfFuelBought[i - 1]
                          for i in range(1, nCities)), name="amount of fuel can be bought")

        # Set the OutputFlag parameter to 0 to suppress solver output
        model.Params.OutputFlag = 0

        # Solve
        model.optimize()

        # Timer
        end_time = time_calculator.time()
        elapsed_time = end_time - start_time
        #print(f"Elapsed time: {elapsed_time:.2f} seconds")

        # Record
        total_distance = 0.0
        toll_cost = 0.0
        total_time = 0.0
        for i in range(1, nCities + 1):
            for j in range(1, nCities + 1):
                if i != j:
                    if x[i, j].X > 0:
                        total_distance = total_distance + distances_dict[(i, j)] * x[i, j].X
                        toll_cost = toll_cost + time_dict[(i, j)] * x[i, j].X
                        total_time = total_time + toll_dict[(i, j)] * x[i, j].X

        # Record Path
        path = []
        # Start at the source city
        current_city = sourceID
        # Follow the edges of the TSP tour until we reach the destination
        while current_city != destinationID:
            # Find the next city in the tour
            for j in range(1, nCities + 1):
                if j != source and j != current_city:
                    if x[current_city, j].X > 0.0:
                        # Add the current city to the path
                        path.append(current_city)
                        # Move to the next city
                        current_city = j
                        break

            # Add the destination city to the path
        path.append(destination)

        # save Fuel
        total_fuel_bought = 0
        for i in range(1, nCities + 1):
            if amountOfFuelBought[i - 1].X > 0:
                total_fuel_bought = total_fuel_bought + amountOfFuelBought[i - 1].X

        # save Fuel
        fuel_cities = {}
        total_fuel_bought = 0
        for i in range(1, nCities + 1):
            if amountOfFuelBought[i - 1].X > 0:
                total_fuel_bought = total_fuel_bought + amountOfFuelBought[i - 1].X
                fuel_cities[i] = amountOfFuelBought[i - 1].X

        status = model.Status
        if status == GRB.UNBOUNDED:
            # print('The model cannot be solved because it is unbounded')
            return 'The model cannot be solved because it is unbounded'
            # sys.exit(0)
        if status == GRB.OPTIMAL:
            # print('The optimal objective is %g' % model.ObjVal)
            return model.ObjVal, total_distance, toll_cost, total_time, path, total_fuel_bought, fuel_cities
        if status != GRB.INF_OR_UNBD and status != GRB.INFEASIBLE:
            # print('Optimization was stopped with status %d' % status)
            return 'Optimization was stopped with status %d' % status
            # sys.exit(0)

    except gp.GurobiError as e:
        print('Error code ' + str(e.errno) + ': ' + str(e))

    except AttributeError:
        print('Encountered an attribute error')


def solve_for_all():
    nCities = len(exact_cities)
    pairs = [(source, destination) for source in range(1, nCities + 1) for destination in range(1, nCities + 1) if
             source != destination]
    for source, destination in pairs:
        print(exact_solve(source, destination))


#solve_for_all()

def solve_for_ids(source, destination):
    obj, distance, time, toll, path, fuel, fuel_cities = exact_solve(source, destination)
    return  obj, distance, time, toll, path, fuel, fuel_cities

#solve_for_ids(4, 10)


def saveOptimalSolutionOutput():
    nCities = len(exact_cities)
    pairs = [(source, destination) for source in range(1, nCities + 1) for destination in range(1, nCities + 1) if
             source != destination]

    results = []
    for source, destination in pairs:
        optimal_cost, optimal_path, toll, time, path, fuel, city_amount = exact_solve(source, destination)
        results.append([source, destination, optimal_cost, optimal_path, toll, time, path, fuel, city_amount])

    # Convert the list to a DataFrame
    df_ExactSolutionMatrix = pd.DataFrame(results)

    # define csv file name
    sol_file_name = 'optimal_cost_' + now.strftime("%Y-%m-%d_%H-%M-%S") + '.csv'
    df_ExactSolutionMatrix.columns = ['Source', 'Destination', 'Optimal Cost', 'Distance', 'Toll', 'Time', 'Path', 'Fuel', 'City Amount']

    if not os.path.exists(sol_opt_directory):
        os.makedirs(sol_opt_directory)
    # writing data frame to a CSV file
    df_ExactSolutionMatrix.to_csv(os.path.join(sol_opt_directory, sol_file_name), sep=';', index=False)

def saveOptimalSolution():
    readexactAsCSV()
    saveOptimalSolutionOutput()

#saveOptimalSolution()