import csv
import folium
import pandas as pd
import os
import webbrowser
from dijkstra_solution import *


def wanted_date():
    options = os.listdir(os.path.join("Outputs", ))
    print(options)

    date = input('date: ') 
    directory_read = os.path.join("Outputs" ,  date) 
    return date,directory_read

def readcitiesOutput(date, directory_read):
    cities_file_name = 'cities_' + date + '.csv'

    # construct the full file path using the 'directory' variable
    file_path = os.path.join(directory_read, cities_file_name)

    # read the CSV file into a DataFrame
    df = pd.read_csv(file_path, sep=';')
    # Create a map centered on the mean latitude and longitude
    m = folium.Map(location=[df['latitude'].mean(), df['longitude'].mean()], zoom_start=2)

    start_city = int(input('Starting from: '))
    end_city = int(input('Going to: '))
    #readPathDistanceOutput(date, directory_read)
    readfuelPriceOutput(date, directory_read)
    readvelocityOutput(date, directory_read)
    readfuelConsumptionOutput(date, directory_read)
    readtollPriceOutput(date, directory_read)
    #readpathsOutput(date, directory_read)
    readdistancesOutput(date, directory_read)
    readtimeOutput(date, directory_read)
    readtotalCostOutput(date, directory_read)
    '''
    def calculate_city2city_cost(path, costs):
    
        cost_list=[]
        
        for i in range(len(path)-1):
            city_id1 = path[i]
            city_id2 = path[i+1]
            
            # Look up the total_cost between the two cities
            for (city_1,city_2, total_costs) in costs:
                if (city_id1, city_id2) == (city_1,city_2) or (city_id2, city_id1) == (city_1,city_2):
                    cost_list.append(total_costs)
        return cost_list
    
    def calculate_city2city_distance(path, distances):
        distance_list = []
        
        for i in range(len(path)-1):
            city_id1 = path[i]
            city_id2 = path[i+1]
            
            # Look up the distance between the two cities
            for (city_1,city_2, distance) in distances:
                if (city_id1, city_id2) == (city_1,city_2) or (city_id2, city_id1) == (city_1,city_2):
                    distance_list.append(distance)
        return distance_list       

    def calculate_city2city_time(path, times):
    
        time_list = []
        
        for i in range(len(path)-1):
            city_id1 = path[i]
            city_id2 = path[i+1]
            
            # Look up the time between the two cities
            for (city_1,city_2, time) in times:
                if (city_id1, city_id2) == (city_1,city_2) or (city_id2, city_id1) == (city_1,city_2):
                    time_list.append(time)
        return time_list

    def calculate_city2city_toll(path, tolls):
        
        toll_list = []
        
        for i in range(len(path)-1):
            city_id1 = path[i]
            city_id2 = path[i+1]
            
            # Look up the toll between the two cities
            for (city_1,city_2, toll) in tolls:
                if (city_id1, city_id2) == (city_1,city_2) or (city_id2, city_id1) == (city_1,city_2):
                    toll_list.append(toll)
        return toll_list                        
       
    def city2city(i,j):
        start_city_id = int(i)
        end_city_id = int(j)
        path = dijkstra_path(start_city_id, end_city_id, all_distances)
        cost = calculate_city2city_cost(path, all_totalCost)
        distance = calculate_city2city_distance(path, all_distances)
        time = calculate_city2city_time(path, all_time)
        toll = calculate_city2city_toll(path, all_TollPrice)
        print(cost)
        print(distance)
        
    city2city(start_city,end_city)
    '''
    path, total_distance, total_time, total_toll, total_cost = sol_cost_with_ids(start_city,end_city)

    # Add markers for each city
    for _, row in df.iterrows():
        color = 'blue'  # set the default marker color to blue
        if row['City ID'] == start_city:
            color = 'green'  # set the marker color to green for the start city
        elif row['City ID'] == end_city:
            color = 'red'  # set the marker color to red for the end city
        elif row['City ID'] in path:
            color = 'purple'  # set the marker color to purple for cities in the path


        folium.Marker(
            location=[row['latitude'], row['longitude']],
            popup=int(row['City ID']),
            icon=folium.Icon(color=color)
        ).add_to(m)
        
    '''    
    distances = []

    # read the distances CSV file into a dictionary
    with open(os.path.join(directory_read, 'distances_' + date + '.csv'), mode='r') as file:
        reader = csv.DictReader(file)
        for row in reader:
            city1 = row['from']
            city2 = row['to']
            distance = row['distances']
            distances[(city1, city2)] = distance
    '''

    # Add lines between cities in the path
    for i in range(len(path)-1):
        city1 = path[i]
        city2 = path[i+1]
        location1 = (df.loc[df['City ID'] == city1]['latitude'].iloc[0], df.loc[df['City ID'] == city1]['longitude'].iloc[0])
        location2 = (df.loc[df['City ID'] == city2]['latitude'].iloc[0], df.loc[df['City ID'] == city2]['longitude'].iloc[0])
        
        folium.PolyLine(locations=[location1, location2], color='purple').add_to(m)
        
        #distance = distances[(city1, city2)]
        #popup_text = "Distance: {:.2f} km".format(distance)
        #folium.PolyLine(locations=[location1, location2], color='purple', popup=popup_text).add_to(m)

    # Display the map
    filename = "map.html"
    m.save(filename)
    webbrowser.open_new_tab(filename)



def readCitiesOutput_distance(date, directory_read):
    cities_file_name = 'cities_' + date + '.csv'

    # construct the full file path using the 'directory' variable
    file_path = os.path.join(directory_read, cities_file_name)

    # read the CSV file into a DataFrame
    df = pd.read_csv(file_path, sep=';')
    # Create a map centered on the mean latitude and longitude
    m = folium.Map(location=[df['latitude'].mean(), df['longitude'].mean()], zoom_start=2)

    start_city = int(input('Starting from: '))
    end_city = int(input('Going to: '))
    readfuelPriceOutput(date, directory_read)
    readvelocityOutput(date, directory_read)
    readfuelConsumptionOutput(date, directory_read)
    readtollPriceOutput(date, directory_read)
    readdistancesOutput(date, directory_read)
    readtimeOutput(date, directory_read)
    readtotalCostOutput(date, directory_read)

    path, total_distance, total_time, total_toll, total_cost = sol_distances_with_ids(start_city,end_city)

    # Add markers for each city
    for _, row in df.iterrows():
        color = 'blue'  # set the default marker color to blue
        if row['City ID'] == start_city:
            color = 'green'  # set the marker color to green for the start city
        elif row['City ID'] == end_city:
            color = 'red'  # set the marker color to red for the end city
        elif row['City ID'] in path:
            color = 'purple'  # set the marker color to purple for cities in the path


        folium.Marker(
            location=[row['latitude'], row['longitude']],
            popup=int(row['City ID']),
            icon=folium.Icon(color=color)
        ).add_to(m)

    # Add lines between cities in the path
    for i in range(len(path)-1):
        city1 = path[i]
        city2 = path[i+1]
        location1 = (df.loc[df['City ID'] == city1]['latitude'].iloc[0], df.loc[df['City ID'] == city1]['longitude'].iloc[0])
        location2 = (df.loc[df['City ID'] == city2]['latitude'].iloc[0], df.loc[df['City ID'] == city2]['longitude'].iloc[0])
        
        folium.PolyLine(locations=[location1, location2], color='purple').add_to(m)
    
    # Display the map
    filename = "map.html"
    m.save(filename)
    webbrowser.open_new_tab(filename)



def readCitiesFromCSV():
    date,directory_read = wanted_date()
    readcitiesOutput(date,directory_read) 

#readCitiesFromCSV()      

def readCitiesFromDistance():
    date,directory_read = wanted_date()
    readCitiesOutput_distance(date,directory_read) 
 
#readCitiesFromDistance()



def guiMapDijkstraCost(date, directory_read, start_city, end_city):
    cities_file_name = 'cities_' + date + '.csv'

    # construct the full file path using the 'directory' variable
    file_path = os.path.join(directory_read, cities_file_name)

    # read the CSV file into a DataFrame
    df = pd.read_csv(file_path, sep=';')
    # Create a map centered on the mean latitude and longitude
    m = folium.Map(location=[df['latitude'].mean(), df['longitude'].mean()], zoom_start=2)

    
    #readPathDistanceOutput(date, directory_read)
    readfuelPriceOutput(date, directory_read)
    readvelocityOutput(date, directory_read)
    readfuelConsumptionOutput(date, directory_read)
    readtollPriceOutput(date, directory_read)
    #readpathsOutput(date, directory_read)
    readdistancesOutput(date, directory_read)
    readtimeOutput(date, directory_read)
    readtotalCostOutput(date, directory_read)
    '''
    def calculate_city2city_cost(path, costs):
    
        cost_list=[]
        
        for i in range(len(path)-1):
            city_id1 = path[i]
            city_id2 = path[i+1]
            
            # Look up the total_cost between the two cities
            for (city_1,city_2, total_costs) in costs:
                if (city_id1, city_id2) == (city_1,city_2) or (city_id2, city_id1) == (city_1,city_2):
                    cost_list.append(total_costs)
        return cost_list
    
    def calculate_city2city_distance(path, distances):
        distance_list = []
        
        for i in range(len(path)-1):
            city_id1 = path[i]
            city_id2 = path[i+1]
            
            # Look up the distance between the two cities
            for (city_1,city_2, distance) in distances:
                if (city_id1, city_id2) == (city_1,city_2) or (city_id2, city_id1) == (city_1,city_2):
                    distance_list.append(distance)
        return distance_list       

    def calculate_city2city_time(path, times):
    
        time_list = []
        
        for i in range(len(path)-1):
            city_id1 = path[i]
            city_id2 = path[i+1]
            
            # Look up the time between the two cities
            for (city_1,city_2, time) in times:
                if (city_id1, city_id2) == (city_1,city_2) or (city_id2, city_id1) == (city_1,city_2):
                    time_list.append(time)
        return time_list

    def calculate_city2city_toll(path, tolls):
        
        toll_list = []
        
        for i in range(len(path)-1):
            city_id1 = path[i]
            city_id2 = path[i+1]
            
            # Look up the toll between the two cities
            for (city_1,city_2, toll) in tolls:
                if (city_id1, city_id2) == (city_1,city_2) or (city_id2, city_id1) == (city_1,city_2):
                    toll_list.append(toll)
        return toll_list                        
       
    def city2city(i,j):
        start_city_id = int(i)
        end_city_id = int(j)
        path = dijkstra_path(start_city_id, end_city_id, all_distances)
        cost = calculate_city2city_cost(path, all_totalCost)
        distance = calculate_city2city_distance(path, all_distances)
        time = calculate_city2city_time(path, all_time)
        toll = calculate_city2city_toll(path, all_TollPrice)
        print(cost)
        print(distance)
        
    city2city(start_city,end_city)
    '''
    path, total_distance, total_time, total_toll, total_cost = sol_cost_with_ids(start_city,end_city)

    # Add markers for each city
    for _, row in df.iterrows():
        color = 'blue'  # set the default marker color to blue
        if row['City ID'] == start_city:
            color = 'green'  # set the marker color to green for the start city
        elif row['City ID'] == end_city:
            color = 'red'  # set the marker color to red for the end city
        elif row['City ID'] in path:
            color = 'purple'  # set the marker color to purple for cities in the path


        folium.Marker(
            location=[row['latitude'], row['longitude']],
            popup=int(row['City ID']),
            icon=folium.Icon(color=color)
        ).add_to(m)
        
    '''    
    distances = []

    # read the distances CSV file into a dictionary
    with open(os.path.join(directory_read, 'distances_' + date + '.csv'), mode='r') as file:
        reader = csv.DictReader(file)
        for row in reader:
            city1 = row['from']
            city2 = row['to']
            distance = row['distances']
            distances[(city1, city2)] = distance
    '''

    # Add lines between cities in the path
    for i in range(len(path)-1):
        city1 = path[i]
        city2 = path[i+1]
        location1 = (df.loc[df['City ID'] == city1]['latitude'].iloc[0], df.loc[df['City ID'] == city1]['longitude'].iloc[0])
        location2 = (df.loc[df['City ID'] == city2]['latitude'].iloc[0], df.loc[df['City ID'] == city2]['longitude'].iloc[0])
        
        folium.PolyLine(locations=[location1, location2], color='purple').add_to(m)
        
        #distance = distances[(city1, city2)]
        #popup_text = "Distance: {:.2f} km".format(distance)
        #folium.PolyLine(locations=[location1, location2], color='purple', popup=popup_text).add_to(m)

    
    # Display the map
    filename = "map.html"
    m.save(filename)
    webbrowser.open_new_tab(filename)



def guiMapDijkstraDistance(date, directory_read, start_city, end_city):
    cities_file_name = 'cities_' + date + '.csv'

    # construct the full file path using the 'directory' variable
    file_path = os.path.join(directory_read, cities_file_name)

    # read the CSV file into a DataFrame
    df = pd.read_csv(file_path, sep=';')
    # Create a map centered on the mean latitude and longitude
    m = folium.Map(location=[df['latitude'].mean(), df['longitude'].mean()], zoom_start=2)


    readfuelPriceOutput(date, directory_read)
    readvelocityOutput(date, directory_read)
    readfuelConsumptionOutput(date, directory_read)
    readtollPriceOutput(date, directory_read)
    readdistancesOutput(date, directory_read)
    readtimeOutput(date, directory_read)
    readtotalCostOutput(date, directory_read)

    path, total_distance, total_time, total_toll, total_cost = sol_distances_with_ids(start_city,end_city)

    # Add markers for each city
    for _, row in df.iterrows():
        color = 'blue'  # set the default marker color to blue
        if row['City ID'] == start_city:
            color = 'green'  # set the marker color to green for the start city
        elif row['City ID'] == end_city:
            color = 'red'  # set the marker color to red for the end city
        elif row['City ID'] in path:
            color = 'purple'  # set the marker color to purple for cities in the path


        folium.Marker(
            location=[row['latitude'], row['longitude']],
            popup=int(row['City ID']),
            icon=folium.Icon(color=color)
        ).add_to(m)

    # Add lines between cities in the path
    for i in range(len(path)-1):
        city1 = path[i]
        city2 = path[i+1]
        location1 = (df.loc[df['City ID'] == city1]['latitude'].iloc[0], df.loc[df['City ID'] == city1]['longitude'].iloc[0])
        location2 = (df.loc[df['City ID'] == city2]['latitude'].iloc[0], df.loc[df['City ID'] == city2]['longitude'].iloc[0])
        
        folium.PolyLine(locations=[location1, location2], color='purple').add_to(m)
    
    # Display the map
    filename = "map.html"
    m.save(filename)
    webbrowser.open_new_tab(filename)