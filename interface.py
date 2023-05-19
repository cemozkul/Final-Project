from random_instance import*
from dijkstra_solution import*
import tkinter as tk
import pandas as pd
from tkinter import ttk

def wanted_date():
    options = os.listdir(os.path.join("Outputs", ))
    print(options)

    date = input('date: ')
    directory_read = os.path.join("Outputs", date)
    return date, directory_read

def readOutput(date, directory_read):
    output_file_name = 'solution_cost_' + date + '.csv'

    # construct the full file path using the 'directory' variable
    file_path = os.path.join(directory_read, output_file_name)

    # read the CSV file into a DataFrame
    df = pd.read_csv(file_path, sep=';')
    
    return df

root = tk.Tk()
root.title("International Truck Routing")

source_label = ttk.Label(root, text="Source:")
source_label.pack()

source_city = ttk.Entry(root)
source_city.pack()

dest_label = ttk.Label(root, text="Destination:")
dest_label.pack()

dest_city = ttk.Entry(root)
dest_city.pack()

result_label = ttk.Label(root, text="")
result_label.pack()

def calculate_cost():
    try:
        global df
        source = int(source_city.get())
        dest = int(dest_city.get())
        row = df[(df["from"] == source) & (df["to"] == dest)]
        if len(row) == 0:
            result_label.config(text="No way.")
        else:
            cost = row.iloc[0]["Total Cost"]
            path_distance = row.iloc[0]["pathDistance"]
            path = row.iloc[0]["path"]
            result_label.config(text=f"Total Cost: {cost}\nPath Distance: {path_distance}\nPath: {path}")
    except ValueError:
        result_label.config(text="Please enter valid numbers.")
    
calculate_button = ttk.Button(root, text="Calculate", command=calculate_cost)
calculate_button.pack()

# call the function to read the output file and store the DataFrame
date, directory_read = wanted_date()
df = readOutput(date, directory_read)

root.mainloop()

