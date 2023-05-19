from map_exact import ExactMap
from random_instance import *
import sys
from business_definition import *
from solution_nn import *
from solution import *
from map_dijkstra import*
from model_solution import *
from dijkstra_solution import *
from map_solnn import*
from exact_solution import*

def userInput():
    if len(sys.argv) < 2:
        print("Usage: python3 logistics.py <generate/help/outputs/inputs/solutionNN/saveNN/solution/saveExactSolution/mapNNCost/mapNNDistance/mapDistance/mapCost/mapExact >")
        print("Type 'help' for business definition \nType 'inputs' for inputs of the problem \nType 'outputs' for outsputs of the problem \nType 'generate' to save all output files as csv files \nType 'costSolution' to see the solution that calculates the minimum cost using Dijkstra algorithm \nType 'saveSolution' to save the solution that calculates the minimum cost using Dijkstra algorithm \nType 'map' to see cities on the map.")
        print("Type 'distanceSolution' to see the solution that calculates the minimum distance using Dijkstra algorithm \nType 'saveSolutionDistance' to save the solution that calculates the minimum distance using Dijkstra algorithm \nType 'saveNN' to save the solutions that calculates the minimum cost and minimum distance using Nearest Neighbor algorithm \nType 'mapNNDistance' to see solution using Nearest Neighbor algorithm on the map.\nType 'mapNNCost' to see solution using Nearest Neighbor algorithm on the map.")
        sys.exit(1)
    
    userInput1 = sys.argv[1]
    
    if userInput1 == 'generate' :
        if len(sys.argv) != 2:
            print('Usage: python3 logistics.py <save>')
            sys.exit(1)
            
        generate_random()    
        saveAllAsCSV()
        print('All of them are saved as csv files in Outputs folder.')
      
    elif userInput1 == 'help':
        if len(sys.argv) != 2:
            print('Usage: python3 logistics.py <help>')
            sys.exit(1)
        help()

    elif userInput1 == 'inputs':
        if len(sys.argv) != 2:
            print('Usage: python3 logistics.py <inputs>')
            sys.exit(1)
        inputs()

    elif userInput1 == 'outputs':
        if len(sys.argv) != 2:
            print('Usage: python3 logistics.py <outputs>')
            sys.exit(1)
        outputs()

    elif userInput1 == 'costSolution':
        if len(sys.argv) != 2:
            print('Usage: python3 logistics.py <costSolution>')
            sys.exit(1)
        readAllAsCSV()
        sol_for_cost()
    
    elif userInput1 == 'saveSolution':
        if len(sys.argv) != 2:
            print('Usage: python3 logistics.py <saveSolution>')
            sys.exit(1)
        saveSolution()
        print('Solutions are saved as csv files in Solution folder.')   

    elif userInput1 == 'distanceSolution':
        if len(sys.argv) != 2:
            print('Usage: python3 logistics.py <distanceSolution>')
            sys.exit(1)
        readAllAsCSV()
        sol_for_distance()

    elif userInput1 == 'saveSolutionDistance':
        if len(sys.argv) != 2:
            print('Usage: python3 logistics.py <saveSolutionDistance>')
            sys.exit(1)
        saveDistanceSolution()
        print('Solutions for distance are saved as csv files in Solution folder.') 

    elif userInput1 == 'mapCost':
        if len(sys.argv) != 2:
            print('Usage: python3 logistics.py <mapCost>')
            sys.exit(1)
        readCitiesFromCSV()
        print('You can find the map for cost in map.html')

    elif userInput1 == 'mapDistance':
        if len(sys.argv) != 2:
            print('Usage: python3 logistics.py <mapDistance>')
            sys.exit(1)
        readCitiesFromDistance()
        print('You can find the map for distance in map.html')
        
    elif userInput1 == 'nearestNeighbor':
        if len(sys.argv) == 5:
            
            userInput2 = sys.argv[2]   
            i          = sys.argv[3]
            j          = sys.argv[4]
            
            if not sys.argv[3].isdigit() or not sys.argv[4].isdigit():
                    print("Invalid input. Please enter an integer.")
                    sys.exit(1)
            
            if userInput2 == 'distance':
                readNNAsCSV()
                if int(i) > len(ID_List2) or int(j) > len(ID_List2):
                    print(f"Invalid input. Please enter a valid ID number.")
                else:
                   sol_with_ids(i,j)

            if userInput2 == 'cost' :
                readNNAsCSV()
                if int(i) > len(ID_List2) or int(j) > len(ID_List2):
                    print(f"Invalid input. Please enter a valid ID number.")
                else:
                   sol_with_ids_cost(i,j)

            else:
                print("Usage: python3 logistics.py <solution> <NN> <i> <j>") 
                sys.exit(1)
        else:
            print("Usage: python3 logistics.py <solution> <NN> <i> <j>") 
            sys.exit(1)
    
    elif userInput1 == 'saveNN':
        if len(sys.argv) != 2:
            print('Usage: python3 logistics.py <saveNN>')
            sys.exit(1)
        saveNN()
        print('Solutions for NN are saved as csv files in Solutions_NN folder.') 
    
    elif userInput1 == 'mapNNDistance':
        if len(sys.argv) != 2:
            print('Usage: python3 logistics.py <mapNNDistance>')
            sys.exit(1)
        readCitiesFromCSV_nn_distance()
        print('You can find the nearest neighbor map for distance in map_nn.html')
    
    elif userInput1 == 'mapNNCost':
        if len(sys.argv) != 2:
            print('Usage: python3 logistics.py <mapNNCost>')
            sys.exit(1)
        readCitiesFromCSV_nn_cost()
        print('You can find the nearest neighbor map for distance in map_nn.html')

    elif userInput1 == 'saveExactSolution':
        if len(sys.argv) != 2:
            print('Usage: python3 logistics.py <saveExactSolution>')
            sys.exit(1)
        saveOptimalSolution()
        print('Solutions are saved as csv files in Solutions_Exact folder.')

    elif userInput1 == 'mapExact':
        if len(sys.argv) != 2:
            print('Usage: python3 logistics.py <mapExact>')
            sys.exit(1)
        ExactMap()
        print('You can find the Gurobi solution map in map_exact.html')

    else:
        print("Usage: python3 logistics.py <generate/help/outputs/inputs/solutionNN/saveNN/solution/saveExactSolution/mapNNCost/mapNNDistance/mapDistance/mapCost>")      
        sys.exit(1)             
userInput()
    
    
    
    

    
