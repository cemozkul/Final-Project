import csv
import folium
import pandas as pd
import os
import webbrowser
from solution_nn import *

def wanted_date():
    options = os.listdir(os.path.join("Outputs", ))
    print(options)

    date = input('date: ') 
    directory_read = os.path.join("Outputs" ,  date) 
    return date,directory_read


def readcitiesOutput_distance(date, directory_read):
    cities_file_name = 'cities_' + date + '.csv'

    # construct the full file path using the 'directory' variable
    file_path = os.path.join(directory_read, cities_file_name)

    # read the CSV file into a DataFrame
    df = pd.read_csv(file_path, sep=';')
    # Create a map centered on the mean latitude and longitude
    m = folium.Map(location=[df['latitude'].mean(), df['longitude'].mean()], zoom_start=2)

    start_city = int(input('Starting from: '))
    end_city = int(input('Going to: '))
    
    #readvelocityOutput(date, directory_read)
    readNNdistancesOutput(date, directory_read)
    #readTimeOutput(date, directory_read)
    readNNtotalCostOutput(date, directory_read)
    
    path= sol_with_ids(start_city,end_city)
        
    
    

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
        
        #distance = distances[(city1, city2)]
        #popup_text = "Distance: {:.2f} km".format(distance)
        #folium.PolyLine(locations=[location1, location2], color='purple', popup=popup_text).add_to(m)

    # Display the map
    filename = "map_nn.html"
    m.save(filename)
    webbrowser.open_new_tab(filename)



def readcitiesOutput_cost(date, directory_read):
    cities_file_name = 'cities_' + date + '.csv'

    # construct the full file path using the 'directory' variable
    file_path = os.path.join(directory_read, cities_file_name)

    # read the CSV file into a DataFrame
    df = pd.read_csv(file_path, sep=';')
    # Create a map centered on the mean latitude and longitude
    m = folium.Map(location=[df['latitude'].mean(), df['longitude'].mean()], zoom_start=2)

    start_city = int(input('Starting from: '))
    end_city = int(input('Going to: '))
    
    #readvelocityOutput(date, directory_read)
    readNNdistancesOutput(date, directory_read)
    #readTimeOutput(date, directory_read)
    readNNtotalCostOutput(date, directory_read)
    
    path= sol_with_ids_cost(start_city,end_city)
        
    
    

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
        
        #distance = distances[(city1, city2)]
        #popup_text = "Distance: {:.2f} km".format(distance)
        #folium.PolyLine(locations=[location1, location2], color='purple', popup=popup_text).add_to(m)

    # Display the map
    filename = "map_nn.html"
    m.save(filename)
    webbrowser.open_new_tab(filename)

def readCitiesFromCSV_nn_distance():
    date,directory_read = wanted_date()
    readcitiesOutput_distance(date,directory_read) 
    
#readCitiesFromCSV_nn_distance()

def readCitiesFromCSV_nn_cost():
    date,directory_read = wanted_date()
    readcitiesOutput_cost(date,directory_read) 
    
#readCitiesFromCSV_nn_cost()
def guiMapNNCost(date, directory_read, start_city, end_city):
    
    cities_file_name = 'cities_' + date + '.csv'

    # construct the full file path using the 'directory' variable
    file_path = os.path.join(directory_read, cities_file_name)

    # read the CSV file into a DataFrame
    df = pd.read_csv(file_path, sep=';')
    # Create a map centered on the mean latitude and longitude
    m = folium.Map(location=[df['latitude'].mean(), df['longitude'].mean()], zoom_start=2)

    
    #readvelocityOutput(date, directory_read)
    readNNdistancesOutput(date, directory_read)
    #readTimeOutput(date, directory_read)
    readNNtotalCostOutput(date, directory_read)
    
    path= sol_with_ids_cost(start_city,end_city)
        
    
    

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
        
        #distance = distances[(city1, city2)]
        #popup_text = "Distance: {:.2f} km".format(distance)
        #folium.PolyLine(locations=[location1, location2], color='purple', popup=popup_text).add_to(m)

    # Display the map
    filename = "map_nn.html"
    m.save(filename)
    webbrowser.open_new_tab(filename)
    
def guiMapNNDistance(date, directory_read, start_city, end_city):
    cities_file_name = 'cities_' + date + '.csv'

    # construct the full file path using the 'directory' variable
    file_path = os.path.join(directory_read, cities_file_name)

    # read the CSV file into a DataFrame
    df = pd.read_csv(file_path, sep=';')
    # Create a map centered on the mean latitude and longitude
    m = folium.Map(location=[df['latitude'].mean(), df['longitude'].mean()], zoom_start=2)

    
    #readvelocityOutput(date, directory_read)
    readNNdistancesOutput(date, directory_read)
    #readTimeOutput(date, directory_read)
    readNNtotalCostOutput(date, directory_read)
    
    path= sol_with_ids(start_city,end_city)
        
    
    

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
        
        #distance = distances[(city1, city2)]
        #popup_text = "Distance: {:.2f} km".format(distance)
        #folium.PolyLine(locations=[location1, location2], color='purple', popup=popup_text).add_to(m)

    # Display the map
    filename = "map_nn.html"
    m.save(filename)
    webbrowser.open_new_tab(filename)