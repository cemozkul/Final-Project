<p align="center">
  <h1 align="center"><br></br>🚚<br>International Truck Routing<br></br></h1>

  <p align="center">
The algorithm we proposed was aimed to offer an optimal route for all logistics companies by optimizing in terms of cost and time. The solution to be created is optimized not only for the routes given but also for every route in the world that can be reached by road.
    <br>
  </p>

# Table of Contents
* [Installation](#Installation)
* [Usage](#Usage)
* [Contributors](#Contributors)
* [License](#License)

## Installation
### Prerequisites
- Python version 2.7+

### Install Package
To use map, you should install given package.
```bash
 pip install folium
```
To use exact solution, you should install Gurobi. To install a specific version (eg 9.1.1), use this command.
```bash
 python -m pip install gurobipy==9.1.1
```
## Usage
`HOW TO REACH INTERNATIONAL TRUCK ROUTING REPOSITORY`

Option 1 to reach inputs list:
- Open GitHub Desktop.
- Select International Truck Routing repository.
- Right-click on the Current Repository and select the 'Open in Terminal' option for Mac and select the 'Open in Cmd Prompt' for Windows. This will take you to Cmd Screen/Terminal. You are now inside InternationalTruckRouting.

Option 2 to reach inputs list:
- Open Command Prompt for Windows or Terminal for Mac.
- Please navigate to the folder location where International Truck Routing folder is located on your computer using the "cd" (change directory) command. For example, if you keep International Truck Routing folder in your Desktop, first navigate to Desktop using ''cd Desktop'' command and then when you get inside Desktop navigate to International Truck Routing folder using ''cd InternationalTruckRouting'' command. You are now inside InternationalTruckRouting.


`logistics.py`

To run logistics.py please follow these steps:
1. Please follow the steps provided under HOW TO REACH INTERNATIONAL TRUCK ROUTING REPOSITORY first.
2. Please use below commands according to your need. (If python3 does not work please try python logistics.py) 

 ```text
To generate and save all the csv file                : python3 logistics.py generate 
To see business definition                           : python3 logistics.py help
To see outputs of the problem                        : python3 logistics.py outputs
To see inputs of the problem                         : python3 logistics.py inputs
To see the result of Dijkstra algorithm for cost     : python3 logistics.py costSolution
To see the result of Dijkstra algorithm for distance : python3 logistics.py distanceSolution
To save the result of Dijkstra algorithm             : python3 logistics.py saveSolution
To visualize the result of algorithm                 : python3 logistics.py map
``` 
3. If you don't type anything, or you enter invalid command after python3 logistics.py `Usage: python3 logistics.py <generate/help/outputs/inputs/costSolution/saveSolution/map` -in form of text-will appear in the screen, and you will turn back to step2.

`config.yaml`

The code reads the configuration details from a YAML file named config.yaml. It uses the information from the configuration file to generate random instances of International Truck Routing problem. Please follow the steps provided below to get inside config.yaml and make changes:

Option 1:
1. Please follow the steps provided under HOW TO REACH INTERNATIONAL TRUCK ROUTING REPOSITORY first.
2. Please type nano config.yaml and press Enter. This will guide you inside the config.yaml file.
3. By using arrow keys on your keyboard, navigate to the location where you want to make changes. For example, if you want to change the number of cities from 30 to 20. Please go to number of cities section using arrow keys.
4. Now, you can delete integer values or make changes and updates. 
5. To save your updates and exit please press Control(^) and X at the same time. 
6. The following question will pop up ''Save modified buffer (ANSWERING "No" WILL DESTROY CHANGES)?'' please press Y on your keyboard.
7. Now your changes are saved.

Option 2: (nano command may not work all the time so option 2 is created for an alternative solution)
1. Please follow the steps provided under HOW TO REACH INTERNATIONAL TRUCK ROUTING REPOSITORY first.
2. Please type notepad config.yaml and press Enter. This will guide you inside the notepad that keeps the values of config.yaml file.
3. You can directly make changes on notepad such as deleting and changing numbers. 
4. To save your updates you can just click on the cross button on the top right. 
5. The following option will pop up ''Save the following changes'' please click on ''Save''.
6. Now your changes are saved.

! Changes you made on config.yaml will directly affect the outputs of logistics.py.

`RANDOM INSTANCE GENERATION`

The purpose of this Random Instance Generator is to create a set of random instances of the International Truck Routing Problem, which those random instances will be used for testing and evaluation purposes.


## Contributors

#### Contributors of the International Truck Routing Project

- Doğa Ağdaş
- Bilge Zorgör
- Damla Barın
- Cem Özkul
- Furkan Yılmaz
- Betül Seyhan

## Lisances 
- [The MIT License (MIT)](https://github.com/OEKundakcioglu/InternationalTruckRouting/blob/main/LICENSE)
