from tkinter import *
from random_instance import *
from business_definition import *
from solution_nn import *
from solution import *
from map_dijkstra import*
from model_solution import *
from dijkstra_solution import *

def generate_callback():
    generate_random()    
    saveAllAsCSV()
    result_label.config(text='All of them are saved as csv files in Outputs folder.')

def help_callback():
    help()
    result_label.config(text='')

def inputs_callback():
    inputs()
    result_label.config(text='')

def outputs_callback():
    outputs()
    result_label.config(text='')

def cost_solution_callback():
    readAllAsCSV()
    sol_for_cost()
    result_label.config(text='')

def save_solution_callback():
    readAllAsCSV()
    save_for_cost()
    result_label.config(text='Solutions are saved as csv files in Solution folder.')

def map_callback():
    readCitiesFromCSV()
    result_label.config(text='You can find the map in map.html')

# create the main window
root = Tk()
root.title("Logistics Interface")

# create buttons for each function
generate_button = Button(root, text="Generate", command=generate_callback)
generate_button.pack()

help_button = Button(root, text="Help", command=help_callback)
help_button.pack()

inputs_button = Button(root, text="Inputs", command=inputs_callback)
inputs_button.pack()

outputs_button = Button(root, text="Outputs", command=outputs_callback)
outputs_button.pack()

cost_solution_button = Button(root, text="Cost Solution", command=cost_solution_callback)
cost_solution_button.pack()

save_solution_button = Button(root, text="Save Solution", command=save_solution_callback)
save_solution_button.pack()

map_button = Button(root, text="Map", command=map_callback)
map_button.pack()

# create a label to display the results
result_label = Label(root, text='')
result_label.pack()

# start the event loop
root.mainloop()