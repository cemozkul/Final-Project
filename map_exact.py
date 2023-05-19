import csv
import folium
import pandas as pd
import os
import webbrowser
from exact_solution import *


def wanted_date():
    options = os.listdir(os.path.join("Outputs", ))
    print(options)

    date = input('date: ')
    directory_read = os.path.join("Outputs", date)
    return date, directory_read


def readcitiesOutput(date, directory_read):
    cities_file_name = 'cities_' + date + '.csv'

    # construct the full file path using the 'directory' variable
    file_path = os.path.join(directory_read, cities_file_name)

    # read the CSV file into a DataFrame
    df = pd.read_csv(file_path, sep=';')
    # Create a map centered on the mean latitude and longitude
    m = folium.Map(location=[df['latitude'].mean(), df['longitude'].mean()], zoom_start=2)

    source = int(input('Starting from: '))
    destination = int(input('Going to: '))
    readFuelPriceOutput(date, directory_read)
    readVelocityOutput(date, directory_read)
    readFuelConsumptionOutput(date, directory_read)
    readTollPriceOutput(date, directory_read)
    readDistancesOutput(date, directory_read)
    readTimeOutput(date, directory_read)
    readTotalCostOutput(date, directory_read)

    path = solve_for_ids(source, destination)

    # Add markers for each city
    for _, row in df.iterrows():
        color = 'blue'  # set the default marker color to blue
        if row['City ID'] == source:
            color = 'green'  # set the marker color to green for the start city
        elif row['City ID'] == destination:
            color = 'red'  # set the marker color to red for the end city
        elif row['City ID'] in path:
            color = 'purple'  # set the marker color to purple for cities in the path

        folium.Marker(
            location=[row['latitude'], row['longitude']],
            popup=int(row['City ID']),
            icon=folium.Icon(color=color)
        ).add_to(m)

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
    filename = "map_exact.html"
    m.save(filename)
    webbrowser.open_new_tab(filename)


def ExactMap():
    date, directory_read = wanted_date()
    readcitiesOutput(date, directory_read)
    


def guiMapExact(date, directory_read, source, destination):
    cities_file_name = 'cities_' + date + '.csv'

    # construct the full file path using the 'directory' variable
    file_path = os.path.join(directory_read, cities_file_name)

    # read the CSV file into a DataFrame
    df = pd.read_csv(file_path, sep=';')
    # Create a map centered on the mean latitude and longitude
    m = folium.Map(location=[df['latitude'].mean(), df['longitude'].mean()], zoom_start=2)

    
    readFuelPriceOutput(date, directory_read)
    readVelocityOutput(date, directory_read)
    readFuelConsumptionOutput(date, directory_read)
    readTollPriceOutput(date, directory_read)
    readDistancesOutput(date, directory_read)
    readTimeOutput(date, directory_read)
    readTotalCostOutput(date, directory_read)

    path = solve_for_ids(source, destination)

    # Add markers for each city
    for _, row in df.iterrows():
        color = 'blue'  # set the default marker color to blue
        if row['City ID'] == source:
            color = 'green'  # set the marker color to green for the start city
        elif row['City ID'] == destination:
            color = 'red'  # set the marker color to red for the end city
        elif row['City ID'] in path:
            color = 'purple'  # set the marker color to purple for cities in the path

        folium.Marker(
            location=[row['latitude'], row['longitude']],
            popup=int(row['City ID']),
            icon=folium.Icon(color=color)
        ).add_to(m)

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
    filename = "map_exact.html"
    m.save(filename)
    webbrowser.open_new_tab(filename)