import tkinter
import customtkinter
from random_instance import *
from business_definition import *
from solution_nn import *
from solution import *
from map_dijkstra import*
from map_exact import*
from map_solnn import*
from model_solution import *
from dijkstra_solution import *
from exact_solution import *
from tkinter import messagebox 
import os
from tkinter import filedialog
import tkinter as tk
import time

#Generate
def generate_callback():
    generate_random()    
    saveAllAsCSV()  
    messagebox.showinfo("Success", "All files saved as CSV in Outputs folder.")

#General Information
def help_callback():
    help()

def inputs_callback():
    inputs()
    
def outputs_callback():
    outputs()




#Dijkstra Solution    
#COST
def dij_cost_solution_callback():
    # Create a new Tkinter window
    input_window = tk.Tk()
    
    # Create and pack the labels and entry fields for start_city and end_city
    start_city_label = tk.Label(input_window, text="Start City:")
    start_city_label.pack()
    start_city_entry = tk.Entry(input_window)
    start_city_entry.pack()

    end_city_label = tk.Label(input_window, text="End City:")
    end_city_label.pack()
    end_city_entry = tk.Entry(input_window)
    end_city_entry.pack()

    # Create and pack the submit button
    submit_button = tk.Button(input_window, text="Submit", command=lambda: dij_cost_input(int(start_city_entry.get()), int(end_city_entry.get()), input_window))
    submit_button.pack()

    # Start the Tkinter event loop
    input_window.mainloop()
    
def dij_cost_input(source, destination, input_window):
    
    
    
    # Close the input window
    input_window.destroy()

    # Open filedialog window for user to select a date
    root = tk.Tk()
    root.withdraw()  # makes the popup window not visible
    directory_read = filedialog.askdirectory(initialdir="Outputs", title="Select a date")

    
    # If no date is selected, terminate the process
    if not directory_read:
        return
    date = directory_read[-19:]
    
    readcitiesOutput(date, directory_read)
    readdistancesOutput(date, directory_read)
    readvelocityOutput(date, directory_read)
    readtimeOutput(date, directory_read)
    readtollPriceOutput(date, directory_read)
    readtotalCostOutput(date, directory_read)
    readfuelPriceOutput(date, directory_read)
    readfuelConsumptionOutput(date, directory_read)
    finding_fuel_needed()
    start_city_id = source
    end_city_id = destination    
    path = dijkstra_path(start_city_id, end_city_id, all_totalCost, all_distances)
    total_distance = calculate_path_distance(path, all_distances) 
    total_time = calculate_path_time(path, all_time) 
    total_toll = calculate_path_toll(path, all_TollPrice)
    total_cost = calculate_path_total_cost(path, all_totalCost)
    fuel_purchase , fuel_cost = minimize_fuel_cost(path, fuel_needed, all_fuelPrice, F, Sf)
    

    
    
    messagebox.showinfo("Info", f"Path: {path} \n Total Distance: {total_distance:.2f} km \n Total Time: {total_time} hours \n Total Toll: €{total_toll:.2f} \n Total Cost: €{total_cost + fuel_cost:,.2f} \n Fuel Buy In/Amount: {fuel_purchase} \n Fuel Cost: €{fuel_cost:.2f}")

def map_dij_cost_callback():
    # Create a new Tkinter window
    input_window = tk.Tk()
    
    # Create and pack the labels and entry fields for start_city and end_city
    start_city_label = tk.Label(input_window, text="Start City:")
    start_city_label.pack()
    start_city_entry = tk.Entry(input_window)
    start_city_entry.pack()

    end_city_label = tk.Label(input_window, text="End City:")
    end_city_label.pack()
    end_city_entry = tk.Entry(input_window)
    end_city_entry.pack()

    # Create and pack the submit button
    submit_button = tk.Button(input_window, text="Submit", command=lambda: map_dij_cost_input(int(start_city_entry.get()), int(end_city_entry.get()), input_window))
    submit_button.pack()

    # Start the Tkinter event loop
    input_window.mainloop()

def map_dij_cost_input(start_city, end_city, input_window):
    # Close the input window
    input_window.destroy()

    # Open filedialog window for user to select a date
    root = tk.Tk()
    root.withdraw()  # makes the popup window not visible
    directory_read = filedialog.askdirectory(initialdir="Outputs", title="Select a date")

    # If no date is selected, terminate the process
    if not directory_read:
        return
    
    
    date = directory_read[-19:]
    
    # Call the function with the user input and selected date
    guiMapDijkstraCost(date, directory_read, start_city, end_city)
    
    messagebox.showinfo("Success", "You can find the Dijkstra map for distance in map.html")

#DISTANCE    
def dij_distance_solution_callback():
    # Create a new Tkinter window
    input_window = tk.Tk()
    
    # Create and pack the labels and entry fields for start_city and end_city
    start_city_label = tk.Label(input_window, text="Start City:")
    start_city_label.pack()
    start_city_entry = tk.Entry(input_window)
    start_city_entry.pack()

    end_city_label = tk.Label(input_window, text="End City:")
    end_city_label.pack()
    end_city_entry = tk.Entry(input_window)
    end_city_entry.pack()

    # Create and pack the submit button
    submit_button = tk.Button(input_window, text="Submit", command=lambda: dij_distance_input(int(start_city_entry.get()), int(end_city_entry.get()), input_window))
    submit_button.pack()

    # Start the Tkinter event loop
    input_window.mainloop()
    
def dij_distance_input(source, destination, input_window):
    
    
    
    # Close the input window
    input_window.destroy()

    # Open filedialog window for user to select a date
    root = tk.Tk()
    root.withdraw()  # makes the popup window not visible
    directory_read = filedialog.askdirectory(initialdir="Outputs", title="Select a date")
    
    
    # If no date is selected, terminate the process
    if not directory_read:
        return
    date = directory_read[-19:]
    
    readcitiesOutput(date, directory_read)
    readdistancesOutput(date, directory_read)
    readvelocityOutput(date, directory_read)
    readtimeOutput(date, directory_read)
    readtollPriceOutput(date, directory_read)
    readtotalCostOutput(date, directory_read)
    readfuelPriceOutput(date, directory_read)
    readfuelConsumptionOutput(date, directory_read)
    finding_fuel_needed()
    start_city_id = source
    end_city_id = destination
    path = dijkstra_distance(start_city_id, end_city_id, all_distances)
    total_distance = calculate_path_distance(path, all_distances) 
    total_time = calculate_path_time(path, all_time) 
    total_toll = calculate_path_toll(path, all_TollPrice)
    total_cost = calculate_path_total_cost(path, all_totalCost)
    fuel_purchase , fuel_cost = minimize_fuel_cost(path, fuel_needed, all_fuelPrice, F, Sf)
    
    
    messagebox.showinfo("Info", f"Path: {path} \n Total Distance: {total_distance} km \n Total Time: {total_time} hours \n Total Toll: €{total_toll} \n Total Cost: €{total_cost + fuel_cost:,.2f} \n Fuel Buy In/Amount: {fuel_purchase} \n Fuel Cost: €{fuel_cost}")

def map_dij_distance_callback():
    # Create a new Tkinter window
    input_window = tk.Tk()
    
    # Create and pack the labels and entry fields for start_city and end_city
    start_city_label = tk.Label(input_window, text="Start City:")
    start_city_label.pack()
    start_city_entry = tk.Entry(input_window)
    start_city_entry.pack()

    end_city_label = tk.Label(input_window, text="End City:")
    end_city_label.pack()
    end_city_entry = tk.Entry(input_window)
    end_city_entry.pack()

    # Create and pack the submit button
    submit_button = tk.Button(input_window, text="Submit", command=lambda: map_distance_cost_input(int(start_city_entry.get()), int(end_city_entry.get()), input_window))
    submit_button.pack()

    # Start the Tkinter event loop
    input_window.mainloop()

def map_distance_cost_input(start_city, end_city, input_window):
    # Close the input window
    input_window.destroy()

    # Open filedialog window for user to select a date
    root = tk.Tk()
    root.withdraw()  # makes the popup window not visible
    directory_read = filedialog.askdirectory(initialdir="Outputs", title="Select a date")

    # If no date is selected, terminate the process
    if not directory_read:
        return
    
    
    date = directory_read[-19:]
    
    # Call the function with the user input and selected date
    guiMapDijkstraDistance(date, directory_read, start_city, end_city)   
    messagebox.showinfo("Success", "You can find the Dijkstra map for distance in map.html") 

#SAVE
def save_dij_solution_callback():
    # Open filedialog window for user to select a date
    root = tk.Tk()
    root.withdraw()  # makes the popup window not visible
    directory_read = filedialog.askdirectory(initialdir="Outputs", title="Select a date")
    
    # If no date is selected, terminate the process
    if not directory_read:
        return
    date = directory_read[-19:]
    
    
    
    readcitiesOutput(date, directory_read)
    readdistancesOutput(date, directory_read)
    readvelocityOutput(date, directory_read)
    readtimeOutput(date, directory_read)
    readtollPriceOutput(date, directory_read)
    readtotalCostOutput(date, directory_read)
    readfuelPriceOutput(date, directory_read)
    readfuelConsumptionOutput(date, directory_read)
    finding_fuel_needed()
    # Save CSV files of solutions found on selected date
    all_cost_output = []
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
                fuel_purchase, fuel_cost = minimize_fuel_cost(path, fuel_needed, all_fuelPrice, F, Sf)

                all_cost_output.append([start_city_id, end_city_id, total_cost + fuel_cost, total_distance, path,
                                        total_toll, total_time, fuel_purchase, fuel_cost])
        

    
    
    def Costsaver():
        df_SolutionMatrix = pd.DataFrame(all_cost_output)
        # define csv file name
        sol_file_name = 'solution_cost_' + now.strftime("%Y-%m-%d_%H-%M-%S") + '.csv'
        # set column names
        df_SolutionMatrix.columns = ['from', 'to', 'Total Cost', 'pathDistance', 'path', 'toll', 'time', 'City-Amount', 'fuel cost']

        
        if not os.path.exists(sol_directory):
            os.makedirs(sol_directory)
        # writing data frame to a CSV file
        df_SolutionMatrix.to_csv(os.path.join(sol_directory, sol_file_name), sep=';', index=False)
        
    
    Costsaver()
    
    all_distance_output = []
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
                fuel_purchase, fuel_cost = minimize_fuel_cost(path, fuel_needed, all_fuelPrice, F, Sf)

                all_distance_output.append([start_city_id, end_city_id, total_cost + fuel_cost, total_distance, path,
                                        total_toll, total_time, fuel_purchase, fuel_cost])
        

    
    
    def Distancesaver():
        df_SolutionMatrix = pd.DataFrame(all_distance_output)
        # define csv file name
        sol_file_name = 'solution_distance_' + now.strftime("%Y-%m-%d_%H-%M-%S") + '.csv'
        # set column names
        df_SolutionMatrix.columns = ['from', 'to', 'Total Cost', 'pathDistance', 'path', 'toll', 'time', 'City-Amount', 'fuel cost']

        
        if not os.path.exists(sol_directory):
            os.makedirs(sol_directory)
        # writing data frame to a CSV file
        df_SolutionMatrix.to_csv(os.path.join(sol_directory, sol_file_name), sep=';', index=False)
        messagebox.showinfo("Info", "Solutions are saved as csv files in Solution folder.")
    
    Distancesaver()





#NN Solution
#COST
def NN_cost_solution_callback():
    # Create a new Tkinter window
    input_window = tk.Tk()
    
    # Create and pack the labels and entry fields for start_city and end_city
    start_city_label = tk.Label(input_window, text="Start City:")
    start_city_label.pack()
    start_city_entry = tk.Entry(input_window)
    start_city_entry.pack()

    end_city_label = tk.Label(input_window, text="End City:")
    end_city_label.pack()
    end_city_entry = tk.Entry(input_window)
    end_city_entry.pack()

    # Create and pack the submit button
    submit_button = tk.Button(input_window, text="Submit", command=lambda: NN_cost_input(int(start_city_entry.get()), int(end_city_entry.get()), input_window))
    submit_button.pack()

    # Start the Tkinter event loop
    input_window.mainloop()
    
def NN_cost_input(source, destination, input_window):
    
    
    
    # Close the input window
    input_window.destroy()

    # Open filedialog window for user to select a date
    root = tk.Tk()
    root.withdraw()  # makes the popup window not visible
    directory_read = filedialog.askdirectory(initialdir="Outputs", title="Select a date")


    
    # If no date is selected, terminate the process
    if not directory_read:
        return
    date = directory_read[-19:]
    
    readNNcitiesOutput(date, directory_read)
    readNNdistancesOutput(date, directory_read)
    readNNtotalCostOutput(date, directory_read)
    readNNFuelPriceOutput(date, directory_read)
    readNNFuelConsumptionOutput(date, directory_read)
    readNNtimeOutput(date, directory_read)
    readNNtollPriceOutput(date, directory_read)
    NN_finding_fuel_needed()
    start_city_id  = source
    end_city_id = destination
    path = nearest_neighbor_path_cost(start_city_id, end_city_id, NN_totalCost, NN_distances)
    total_distance = calculate_path_distance(path, NN_distances)
    total_time = calculate_path_time(path, NN_time) 
    total_toll = calculate_path_toll(path, NN_toll) 
    total_cost = calculate_path_total_cost(path, NN_totalCost)  
    fuel_purchase , fuel_cost = minimize_fuel_cost(path, NN_fuel_needed, NN_fuelPrice, F, Sf) 
    
    
    
    messagebox.showinfo("Info", f"Path: {path} \n Total Distance: {total_distance} km \n Total Time: {total_time} hours \n Total Toll: €{total_toll} \n Total Cost: €{total_cost + fuel_cost:,.2f} \n Fuel Buy In/Amount: {fuel_purchase} \n Fuel Cost: €{fuel_cost}")

def map_NN_cost_callback():
    # Create a new Tkinter window
    input_window = tk.Tk()
    
    # Create and pack the labels and entry fields for start_city and end_city
    start_city_label = tk.Label(input_window, text="Start City:")
    start_city_label.pack()
    start_city_entry = tk.Entry(input_window)
    start_city_entry.pack()

    end_city_label = tk.Label(input_window, text="End City:")
    end_city_label.pack()
    end_city_entry = tk.Entry(input_window)
    end_city_entry.pack()

    # Create and pack the submit button
    submit_button = tk.Button(input_window, text="Submit", command=lambda: map_NN_cost_input(int(start_city_entry.get()), int(end_city_entry.get()), input_window))
    submit_button.pack()

    # Start the Tkinter event loop
    input_window.mainloop()    

def map_NN_cost_input(start_city, end_city, input_window):
    # Close the input window
    input_window.destroy()

    # Open filedialog window for user to select a date
    root = tk.Tk()
    root.withdraw()  # makes the popup window not visible
    directory_read = filedialog.askdirectory(initialdir="Outputs", title="Select a date")

    # If no date is selected, terminate the process
    if not directory_read:
        return
    
    
    date = directory_read[-19:]
    
    # Call the function with the user input and selected date
    guiMapNNCost(date, directory_read, start_city, end_city)
    messagebox.showinfo("Success", "You can find the Nearest Neigbor map for distance in map_nn.html")

#DISTANCE
def NN_distance_solution_callback():
    # Create a new Tkinter window
    input_window = tk.Tk()
    
    # Create and pack the labels and entry fields for start_city and end_city
    start_city_label = tk.Label(input_window, text="Start City:")
    start_city_label.pack()
    start_city_entry = tk.Entry(input_window)
    start_city_entry.pack()

    end_city_label = tk.Label(input_window, text="End City:")
    end_city_label.pack()
    end_city_entry = tk.Entry(input_window)
    end_city_entry.pack()

    # Create and pack the submit button
    submit_button = tk.Button(input_window, text="Submit", command=lambda: NN_distance_input(int(start_city_entry.get()), int(end_city_entry.get()), input_window))
    submit_button.pack()

    # Start the Tkinter event loop
    input_window.mainloop()
    
def NN_distance_input(source, destination, input_window):
    
    
    
    # Close the input window
    input_window.destroy()

    # Open filedialog window for user to select a date
    root = tk.Tk()
    root.withdraw()  # makes the popup window not visible
    directory_read = filedialog.askdirectory(initialdir="Outputs", title="Select a date")

    
    # If no date is selected, terminate the process
    if not directory_read:
        return
    date = directory_read[-19:]
    
    readNNcitiesOutput(date, directory_read)
    readNNdistancesOutput(date, directory_read)
    readNNtotalCostOutput(date, directory_read)
    readNNFuelPriceOutput(date, directory_read)
    readNNFuelConsumptionOutput(date, directory_read)
    readNNtimeOutput(date, directory_read)
    readNNtollPriceOutput(date, directory_read)
    NN_finding_fuel_needed()
    start_city_id  = source
    end_city_id = destination 
    path = nearest_neighbor_path(start_city_id, end_city_id, NN_distances)
    total_distance = calculate_path_distance(path, NN_distances) 
    total_time = calculate_path_time(path, NN_time) 
    total_toll = calculate_path_toll(path, NN_toll)
    total_cost = calculate_path_total_cost(path, NN_totalCost)  
    fuel_purchase , fuel_cost = minimize_fuel_cost(path, NN_fuel_needed, NN_fuelPrice, F, Sf)
    
    
    
    
    
    messagebox.showinfo("Info", f"Path: {path} \n Total Distance: {total_distance} km \n Total Time: {total_time} hours \n Total Toll: €{total_toll} \n Total Cost: €{total_cost + fuel_cost:,.2f} \n Fuel Buy In/Amount: {fuel_purchase} \n Fuel Cost: €{fuel_cost}")

def map_NN_distance_callback():
    # Create a new Tkinter window
    input_window = tk.Tk()
    
    # Create and pack the labels and entry fields for start_city and end_city
    start_city_label = tk.Label(input_window, text="Start City:")
    start_city_label.pack()
    start_city_entry = tk.Entry(input_window)
    start_city_entry.pack()

    end_city_label = tk.Label(input_window, text="End City:")
    end_city_label.pack()
    end_city_entry = tk.Entry(input_window)
    end_city_entry.pack()

    # Create and pack the submit button
    submit_button = tk.Button(input_window, text="Submit", command=lambda: map_NN_distance_input(int(start_city_entry.get()), int(end_city_entry.get()), input_window))
    submit_button.pack()

    # Start the Tkinter event loop
    input_window.mainloop()    

def map_NN_distance_input(start_city, end_city, input_window):
    # Close the input window
    input_window.destroy()

    # Open filedialog window for user to select a date
    root = tk.Tk()
    root.withdraw()  # makes the popup window not visible
    directory_read = filedialog.askdirectory(initialdir="Outputs", title="Select a date")
 
    # If no date is selected, terminate the process

    
    if not directory_read:
        return
    
    
    date = directory_read[-19:]
    
    # Call the function with the user input and selected date
    guiMapNNDistance(date, directory_read, start_city, end_city)
    messagebox.showinfo("Success", "You can find the Nearest Neigbor map for distance in map_nn.html")

#SAVE
def save_NN_solution_callback():
   # Open filedialog window for user to select a date
    root = tk.Tk()
    root.withdraw()  # makes the popup window not visible
    directory_read = filedialog.askdirectory(initialdir="Outputs", title="Select a date")

    
    # If no date is selected, terminate the process
    if not directory_read:
        return
    date = directory_read[-19:]
    
    
    
    readNNcitiesOutput(date, directory_read)
    readNNdistancesOutput(date, directory_read)
    readNNtotalCostOutput(date, directory_read)
    readNNFuelPriceOutput(date, directory_read)
    readNNFuelConsumptionOutput(date, directory_read)
    readNNtimeOutput(date, directory_read)
    readNNtollPriceOutput(date, directory_read)
    NN_finding_fuel_needed()
    # Save CSV files of solutions found on selected date
    all_cost_NN_output = []
    IDs = []
    for id in range(len(NN_cities)):
        IDs.append(NN_cities[id][0])

    for i in IDs:
        for j in IDs:
            if i != j:
                start_city_id = int(i)
                end_city_id = int(j)
                path = nearest_neighbor_path_cost(start_city_id, end_city_id, NN_totalCost, NN_distances)
                total_time = calculate_path_time(path, NN_time) 
                total_toll = calculate_path_toll(path, NN_toll)
                total_costs = calculate_path_distance(path, NN_totalCost)  
                total_distance = calculate_path_distance(path, NN_distances)
                fuel_purchase , fuel_cost = minimize_fuel_cost(path, NN_fuel_needed, NN_fuelPrice, F, Sf)
                all_cost_NN_output.append([start_city_id, end_city_id, total_costs + fuel_cost, total_distance, path,
                                        total_toll, total_time, fuel_purchase, fuel_cost])
                
    

    def CostsaverNN():
        df_SolutionNNMatrix = pd.DataFrame(all_cost_NN_output)
        # define csv file name
        solutionNN_file_name = 'solution_NN_cost' + now.strftime("%Y-%m-%d_%H-%M-%S") + '.csv'
        # set column names
        df_SolutionNNMatrix.columns = ['from', 'to', 'Total Cost', 'pathDistance', 'path', 'toll', 'time', 'City-Amount', 'fuel cost']
        # writing data frame to a CSV file
        if not os.path.exists(sol_nn_directory):
            os.makedirs(sol_nn_directory)
        # writing data frame to a CSV file
        df_SolutionNNMatrix.to_csv(os.path.join(sol_nn_directory, solutionNN_file_name), sep=';', index=False)
        
    CostsaverNN()

    all_distance_NN_output = []
    IDs = []
    for id in range(len(NN_cities)):
        IDs.append(NN_cities[id][0])

    for i in IDs:
        for j in IDs:
            if i != j:
                start_city_id = int(i)
                end_city_id = int(j)
                path = nearest_neighbor_path(start_city_id, end_city_id, NN_distances)
                total_time = calculate_path_time(path, NN_time) 
                total_toll = calculate_path_toll(path, NN_toll)
                total_distance = calculate_path_distance(path, NN_distances)  
                total_costs = calculate_path_distance(path, NN_totalCost) 
                fuel_purchase , fuel_cost = minimize_fuel_cost(path, NN_fuel_needed, NN_fuelPrice, F, Sf)
                all_distance_NN_output.append([start_city_id, end_city_id, total_costs + fuel_cost, total_distance, path,
                                        total_toll, total_time, fuel_purchase, fuel_cost])

        

    def Distancesaver():
        df_SolutionNNMatrix = pd.DataFrame(all_distance_NN_output)
        # define csv file name
        solutionNN_file_name = 'solution_NN_distance' + now.strftime("%Y-%m-%d_%H-%M-%S") + '.csv'
        # set column names
        df_SolutionNNMatrix.columns = ['from', 'to', 'Total Cost', 'pathDistance', 'path', 'toll', 'time', 'City-Amount', 'fuel cost']
        # writing data frame to a CSV file
        if not os.path.exists(sol_nn_directory):
            os.makedirs(sol_nn_directory)
        # writing data frame to a CSV file
        df_SolutionNNMatrix.to_csv(os.path.join(sol_nn_directory, solutionNN_file_name), sep=';', index=False)
        messagebox.showinfo("Info", "Solutions are saved as csv files in Solution_NN folder.")
        
    Distancesaver()






#Exact Solution
#SOLUTION
def exact_callback():
     # Create a new Tkinter window
    input_window = tk.Tk()
    
    # Create and pack the labels and entry fields for start_city and end_city
    start_city_label = tk.Label(input_window, text="Start City:")
    start_city_label.pack()
    start_city_entry = tk.Entry(input_window)
    start_city_entry.pack()

    end_city_label = tk.Label(input_window, text="End City:")
    end_city_label.pack()
    end_city_entry = tk.Entry(input_window)
    end_city_entry.pack()

    # Create and pack the submit button
    submit_button = tk.Button(input_window, text="Submit", command=lambda: Exact_input(int(start_city_entry.get()), int(end_city_entry.get()), input_window))
    submit_button.pack()

    # Start the Tkinter event loop
    input_window.mainloop()

def Exact_input(source, destination, input_window):
    
    
    
    # Close the input window
    input_window.destroy()

    # Open filedialog window for user to select a date
    root = tk.Tk()
    root.withdraw()  # makes the popup window not visible
    directory_read = filedialog.askdirectory(initialdir="Outputs", title="Select a date")

    
    # If no date is selected, terminate the process
    if not directory_read:
        return
    date = directory_read[-19:]
    
    readCitiesOutput(date, directory_read)
    readDistancesOutput(date, directory_read)
    readVelocityOutput(date, directory_read)
    readTimeOutput(date, directory_read)
    readTollPriceOutput(date, directory_read)
    readTotalCostOutput(date, directory_read)
    readFuelPriceOutput(date, directory_read)
    readFuelConsumptionOutput(date, directory_read)
    start_city_id  = source
    end_city_id = destination
    obj, distance, time, toll, path, fuel, fuel_cities = exact_solve(start_city_id, end_city_id)
    messagebox.showinfo("Info", f"Path: {path} \n Total Distance: {distance} km \n Total Time: {time} hours \n Total Toll: €{toll} \n Total Cost: €{obj:,.2f} \n Fuel Buy In/Amount: {fuel_cities} \n Fuel Cost: €{fuel}")
    return path
    
    
           

#MAP
def map_Exact_callback():
    # Create a new Tkinter window
    input_window = tk.Tk()
    
    # Create and pack the labels and entry fields for start_city and end_city
    start_city_label = tk.Label(input_window, text="Start City:")
    start_city_label.pack()
    start_city_entry = tk.Entry(input_window)
    start_city_entry.pack()

    end_city_label = tk.Label(input_window, text="End City:")
    end_city_label.pack()
    end_city_entry = tk.Entry(input_window)
    end_city_entry.pack()

    # Create and pack the submit button
    submit_button = tk.Button(input_window, text="Submit", command=lambda: map_Exact_input(int(start_city_entry.get()), int(end_city_entry.get()), input_window))
    submit_button.pack()

    # Start the Tkinter event loop
    input_window.mainloop()    

def map_Exact_input(start_city, end_city, input_window):
    # Close the input window
    input_window.destroy()

    # Open filedialog window for user to select a date
    root = tk.Tk()
    root.withdraw()  # makes the popup window not visible
    directory_read = filedialog.askdirectory(initialdir="Outputs", title="Select a date")

    
    # If no date is selected, terminate the process
    if not directory_read:
        return
    
    
    date = directory_read[-19:]
    
    # Call the function with the user input and selected date
    guiMapExact(date, directory_read, start_city, end_city)
    messagebox.showinfo("Success", "You can find the Exact Solution map for distance in map_exact.html")

#SAVE
def save_Exact_solution_callback():
    # Open filedialog window for user to select a date
    root = tk.Tk()
    root.withdraw()  # makes the popup window not visible
    directory_read = filedialog.askdirectory(initialdir="Outputs", title="Select a date")

    
    # If no date is selected, terminate the process
    if not directory_read:
        return
    date = directory_read[-19:]
    
    
    
    readCitiesOutput(date, directory_read)
    readDistancesOutput(date, directory_read)
    readVelocityOutput(date, directory_read)
    readTimeOutput(date, directory_read)
    readTollPriceOutput(date, directory_read)
    readTotalCostOutput(date, directory_read)
    readFuelPriceOutput(date, directory_read)
    readFuelConsumptionOutput(date, directory_read)
    nCities = len(exact_cities)

    pairs = [(source, destination) for source in range(1, nCities + 1) for destination in range(1, nCities + 1) if
                source != destination]

    results = []
    for source, destination in pairs:
        optimal_cost, optimal_path, time, toll, path, fuel, city_amount = exact_solve(source, destination)
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
    
    messagebox.showinfo("Info", "Solutions are saved as csv files in Solution_Exact folder.")



# System defaults
customtkinter.set_appearance_mode(mode_string="light")
customtkinter.set_default_color_theme(color_string="blue")

# App defaults
app = customtkinter.CTk()
app.geometry(geometry_string="720x720")
app.title("International Truck Routing")

# Adding UI elements
title = customtkinter.CTkLabel(app, text="Welcome! Please select an option that suits your research.")
title.pack(padx=10, pady=10)

# Adding frame for buttons
button_frame = customtkinter.CTkFrame(app)
button_frame.pack(padx=10, pady=10)

# Adding options
option_frame = customtkinter.CTkFrame(app)
option_frame.pack(padx=10, pady=10)

option_var = tkinter.StringVar(value="General Information")
option1 = customtkinter.CTkRadioButton(option_frame, text="General Information", variable=option_var, value="General Information")
option1.pack(side="left", padx=5, pady=5)

option2 = customtkinter.CTkRadioButton(option_frame, text="Dijkstra Solution", variable=option_var, value="Dijkstra Solution")
option2.pack(side="left", padx=5, pady=5)

option3 = customtkinter.CTkRadioButton(option_frame, text="NN Solution", variable=option_var, value="NN Solution")
option3.pack(side="left", padx=5, pady=5)

option4 = customtkinter.CTkRadioButton(option_frame, text="Exact Solution", variable=option_var, value="Exact Solution")
option4.pack(side="left", padx=5, pady=5)

option5 = customtkinter.CTkRadioButton(option_frame, text="Random Generation", variable=option_var, value="Random Generation")
option5.pack(side="left", padx=5, pady=5)

# Adding buttons for General Information
general_info_buttons = [
    customtkinter.CTkButton(button_frame, text="Help", command=help_callback),
    customtkinter.CTkButton(button_frame, text="Inputs", command=inputs_callback),
    customtkinter.CTkButton(button_frame, text="Outputs", command=outputs_callback), 
]

# Adding buttons for Solution
Dijkstra_solution_buttons = [
    customtkinter.CTkButton(button_frame, text="Cost Solution", command=dij_cost_solution_callback),
    customtkinter.CTkButton(button_frame, text="Distance Solution", command=dij_distance_solution_callback),
    customtkinter.CTkButton(button_frame, text="Save Solution", command=save_dij_solution_callback),
    customtkinter.CTkButton(button_frame, text="Cost Solution Map", command=map_dij_cost_callback),
    customtkinter.CTkButton(button_frame, text="Distance Solution Map", command=map_dij_distance_callback)
]

NN_solution_buttons = [
    customtkinter.CTkButton(button_frame, text="Cost Solution", command=NN_cost_solution_callback),
    customtkinter.CTkButton(button_frame, text="Distance Solution", command=NN_distance_solution_callback),
    customtkinter.CTkButton(button_frame, text="Save Solution", command=save_NN_solution_callback),
    customtkinter.CTkButton(button_frame, text="Cost Solution Map", command=map_NN_cost_callback),
    customtkinter.CTkButton(button_frame, text="Distance Solution Map", command=map_NN_distance_callback)
]

Exact_solution_buttons = [
    customtkinter.CTkButton(button_frame, text="Solution", command=exact_callback),
    customtkinter.CTkButton(button_frame, text="Save Solution", command=save_Exact_solution_callback),
    customtkinter.CTkButton(button_frame, text="Solution Map", command=map_Exact_callback)
]

Generate_buttons = [
    customtkinter.CTkButton(button_frame, text="Generate", command=generate_callback)
]


# Function to handle option changes
def handle_option_change(new_option):
    if new_option == "General Information":
        for button in Dijkstra_solution_buttons:
            button.grid_forget()
        for button in NN_solution_buttons:
            button.grid_forget()
        for button in Exact_solution_buttons:
            button.grid_forget()
        for button in Generate_buttons:
            button.grid_forget()
        for i, button in enumerate(general_info_buttons):
            button.grid(row=i//3, column=i%3, padx=10, pady=10)
    elif new_option == "Dijkstra Solution":
        for button in general_info_buttons:
            button.grid_forget()
        for button in NN_solution_buttons:
            button.grid_forget()
        for button in Exact_solution_buttons:
            button.grid_forget()
        for button in Generate_buttons:
            button.grid_forget()
        for i, button in enumerate(Dijkstra_solution_buttons):
            button.grid(row=i//3, column=i%3, padx=10, pady=10)
    elif new_option == "NN Solution":
        for button in general_info_buttons:
            button.grid_forget()
        for button in Dijkstra_solution_buttons:
            button.grid_forget()
        for button in Exact_solution_buttons:
            button.grid_forget()
        for button in Generate_buttons:
            button.grid_forget()
        for i, button in enumerate(NN_solution_buttons):
            button.grid(row=i//3, column=i%3, padx=10, pady=10)
    elif new_option == "Exact Solution":
        for button in general_info_buttons:
            button.grid_forget()
        for button in Dijkstra_solution_buttons:
            button.grid_forget()
        for button in NN_solution_buttons:
            button.grid_forget()
        for button in Generate_buttons:
            button.grid_forget()
        for i, button in enumerate(Exact_solution_buttons):
            button.grid(row=i//3, column=i%3, padx=10, pady=10)
    elif new_option == "Random Generation":
        for button in general_info_buttons:
            button.grid_forget()
        for button in Dijkstra_solution_buttons:
            button.grid_forget()
        for button in NN_solution_buttons:
            button.grid_forget()
        for button in Exact_solution_buttons:
            button.grid_forget()
        for i, button in enumerate(Generate_buttons):
            button.grid(row=i//3, column=i%3, padx=10, pady=10)

# Set initial option to General Information
handle_option_change("General Information")

# Bind option changes to function
option_var.trace("w", lambda *args: handle_option_change(option_var.get()))

# Test to run
app.mainloop()