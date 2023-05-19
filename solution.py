from random_instance import *
import ast


def wanted_date():
    options = os.listdir(os.path.join("Outputs", ))
    print(options)

    date = input('date: ') 
    directory_read = os.path.join("Outputs" ,  date) 
    return date,directory_read

all_time = []
all_pathDistance = []
all_fuelPrice = []
all_velocity = []
all_FuelConsumption = []
all_TollPrice = []
all_paths = []
all_NumberOfPaths = []    
all_distances = []
all_cities = []


def readTimeOutput(date,directory_read):
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
        all_time.append([(from_read[i],to_read[i]),time_read[i]])

def readPathDistanceOutput(date,directory_read):
    pathDistance_file_name = 'pathDistance_' + date + '.csv'

    # construct the full file path using the 'directory' variable
    file_path = os.path.join(directory_read, pathDistance_file_name)

    # read the CSV file into a DataFrame
    df = pd.read_csv(file_path, sep=';')
    from_read = df['from'].tolist()
    to_read = df['to'].tolist()
    # extract the 'time' column as a list
    pathDistance_read = df['pathDistance'].tolist()
    
    for i in range(len(pathDistance_read)):
        all_pathDistance.append([(from_read[i],to_read[i]),pathDistance_read[i]])

def readFuelPriceOutput(date,directory_read):
    fuelPrice_file_name = 'fuel_price_' + date + '.csv'

    # construct the full file path using the 'directory' variable
    file_path = os.path.join(directory_read, fuelPrice_file_name)

    # read the CSV file into a DataFrame
    df = pd.read_csv(file_path, sep=';')
    City_ID_read = df['cityID'].tolist()
    # extract the 'time' column as a list
    fuelPrice_read = df['fuelPrice'].tolist()
    
    for i in range(len(fuelPrice_read)):
        all_fuelPrice.append([City_ID_read[i],fuelPrice_read[i]])

def readvelocityOutput(date,directory_read):
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
        all_velocity.append([(from_read[i],to_read[i]),velocity_read[i]])

def readFuelConsumptionOutput(date,directory_read):
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
        all_FuelConsumption.append([from_read[i],to_read[i],FuelConsumption_read[i]])

def readTollPriceOutput(date,directory_read):
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
        all_TollPrice.append([(from_read[i],to_read[i]),TollPrice_read[i]])

def readpathsOutput(date,directory_read):
    paths_file_name = 'paths_' + date + '.csv'

    # construct the full file path using the 'directory' variable
    file_path = os.path.join(directory_read, paths_file_name)

    # read the CSV file into a DataFrame
    df = pd.read_csv(file_path, sep=';')
    from_read = df['from'].tolist()
    to_read = df['to'].tolist()
    # extract the 'time' column as a list
    paths_read = df['paths'].tolist()
    
    for i in range(len(paths_read)):
        all_paths.append([(from_read[i], to_read[i]),  ast.literal_eval(paths_read[i])])

def readNumberOfPathsOutput(date,directory_read):
    NumberOfPaths_file_name = 'number_of_paths_' + date + '.csv'

    # construct the full file path using the 'directory' variable
    file_path = os.path.join(directory_read, NumberOfPaths_file_name)

    # read the CSV file into a DataFrame
    df = pd.read_csv(file_path, sep=';')
    from_read = df['from'].tolist()
    to_read = df['to'].tolist()
    # extract the 'time' column as a list
    NumberOfPaths_read = df['NumberOfPaths'].tolist()
    
    for i in range(len(NumberOfPaths_read)):
        all_NumberOfPaths.append([(from_read[i],to_read[i]),NumberOfPaths_read[i]])

def readdistancesOutput(date,directory_read):
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
        all_distances.append([(from_read[i],to_read[i]),distances_read[i]])

def readcitiesOutput(date,directory_read):
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
        all_cities.append([cityId_read[i],x_read[i],y_read[i]])
 
def readAllAsCSV():
    date,directory_read = wanted_date()
    readcitiesOutput(date,directory_read)
    readPathDistanceOutput(date,directory_read)
    readFuelPriceOutput(date,directory_read)
    readvelocityOutput(date,directory_read)
    readFuelConsumptionOutput(date,directory_read)
    readTollPriceOutput(date,directory_read)
    readpathsOutput(date,directory_read)
    readNumberOfPathsOutput(date,directory_read)
    readdistancesOutput(date,directory_read)
    readTimeOutput(date,directory_read)
    

def solution():
    #readAllAsCSV()
    city_1 = []
    city_2 = []
    pathDistance2 = []
    path2 = []
    toll2 = []
    time2 = []
    all_output = []    
    IDs = []
    for id in range(len(all_cities)):
        IDs.append(all_cities[id][0])
        

        
    for (city_id1, city_id2), pathDistances in  all_pathDistance:   
        city_1.append(city_id1)   
        city_2.append(city_id2)  
        pathDistance2.append([pathDistances])  
    for (city_id1, city_id2), path in  all_paths:            
        path2.append([path])  
    for (city_id1, city_id2), tollp in all_TollPrice:
        toll2.append([tollp])
    for (city_id1, city_id2), times in all_time:        
        time2.append([times]) 

    for i in range(len(path2)):
        all_output.append([city_1[i], city_2[i], pathDistance2[i][0], path2[i][0], toll2[i][0], time2[i][0]])
    #print(all_output)

    def get_min_pathDistance(all_output , a, b):
        sublists = [sublist for sublist in all_output if sublist[0] == a and sublist[1] == b]
        if len(sublists) == 0:
            return None
        else:
            return min(sublists, key=lambda x: x[2])
        
    def get_min_toll(all_output, a, b):
        sublists = [sublist for sublist in all_output if sublist[0] == a and sublist[1] == b]
        if len(sublists) == 0:
            return None
        else:
            return min(sublists, key=lambda x: x[4])

    def get_min_time(all_output, a, b):
        sublists = [sublist for sublist in all_output if sublist[0] == a and sublist[1] == b]
        if len(sublists) == 0:
            return None
        else:
            return min(sublists, key=lambda x: x[5])


    min_pathDistance_list = []        
    for i in IDs:
            for j in IDs:
                if i != j:
                    min_pathDistance = get_min_pathDistance(all_output, i, j)
                    min_pathDistance_list.append(min_pathDistance)
    
    
    min_toll_list = []        
    for i in IDs:
            for j in IDs:
                if i != j:
                    min_toll = get_min_toll(all_output, i, j)
                    min_toll_list.append(min_toll)
                    
  
                    
    min_time_list = []        
    for i in IDs:
            for j in IDs:
                if i != j:
                    min_time = get_min_time(all_output, i, j)
                    min_time_list.append(min_time)
                       

    def savepathSolutionOutput():
        df_SolutionMatrix = pd.DataFrame(min_pathDistance_list)
        # define csv file name
        sol_file_name = 'solution_pathDistance_' + now.strftime("%Y-%m-%d_%H-%M-%S") + '.csv'
            # set column names
        df_SolutionMatrix.columns = ['from', 'to', 'pathDistance', 'path', 'toll', 'time']

        
        if not os.path.exists(sol_directory):
            os.makedirs(sol_directory)
        # writing data frame to a CSV file
        df_SolutionMatrix.to_csv(os.path.join(sol_directory, sol_file_name), sep=';', index=False)

    def savetollSolutionOutput():
        df_SolutionMatrix = pd.DataFrame(min_toll_list)
        # define csv file name
        sol_file_name = 'solution_toll_' + now.strftime("%Y-%m-%d_%H-%M-%S") + '.csv'
            # set column names
        df_SolutionMatrix.columns = ['from', 'to', 'pathDistance', 'path', 'toll', 'time']


        if not os.path.exists(sol_directory):
            os.makedirs(sol_directory)
        # writing data frame to a CSV file
        df_SolutionMatrix.to_csv(os.path.join(sol_directory, sol_file_name), sep=';', index=False)

    def savetimeSolutionOutput():
        df_SolutionMatrix = pd.DataFrame(min_time_list)
        # define csv file name
        sol_file_name = 'solution_time_' + now.strftime("%Y-%m-%d_%H-%M-%S") + '.csv'
            # set column names
        df_SolutionMatrix.columns = ['from', 'to', 'pathDistance', 'path', 'toll', 'time']

        
        if not os.path.exists(sol_directory):
            os.makedirs(sol_directory)
        # writing data frame to a CSV file
        df_SolutionMatrix.to_csv(os.path.join(sol_directory, sol_file_name), sep=';', index=False)
        
    def AllSol():
        savepathSolutionOutput()
        savetollSolutionOutput()
        savetimeSolutionOutput()   
    
    AllSol() 

def Data():
    #readAllAsCSV()
    city_1 = []
    city_2 = []
    pathDistance2 = []
    path2 = []
    toll2 = []
    time2 = []
    all_output = []    
    IDs = []
    for id in range(len(all_cities)):
        IDs.append(all_cities[id][0])
        

        
    for (city_id1, city_id2), pathDistances in  all_pathDistance:   
        city_1.append(city_id1)   
        city_2.append(city_id2)  
        pathDistance2.append([pathDistances])  
    for (city_id1, city_id2), path in  all_paths:            
        path2.append([path])  
    for (city_id1, city_id2), tollp in all_TollPrice:
        toll2.append([tollp])
    for (city_id1, city_id2), times in all_time:        
        time2.append([times]) 

    for i in range(len(path2)):
        all_output.append([city_1[i], city_2[i], pathDistance2[i][0], path2[i][0], toll2[i][0], time2[i][0]])
    #print(all_output)

    def get_min_pathDistance(all_output , a, b):
        sublists = [sublist for sublist in all_output if sublist[0] == a and sublist[1] == b]
        if len(sublists) == 0:
            return None
        else:
            return min(sublists, key=lambda x: x[2])
        
    def get_min_toll(all_output, a, b):
        sublists = [sublist for sublist in all_output if sublist[0] == a and sublist[1] == b]
        if len(sublists) == 0:
            return None
        else:
            return min(sublists, key=lambda x: x[4])

    def get_min_time(all_output, a, b):
        sublists = [sublist for sublist in all_output if sublist[0] == a and sublist[1] == b]
        if len(sublists) == 0:
            return None
        else:
            return min(sublists, key=lambda x: x[5])


    min_pathDistance_list = []        
    for i in IDs:
            for j in IDs:
                if i != j:
                    min_pathDistance = get_min_pathDistance(all_output, i, j)
                    min_pathDistance_list.append(min_pathDistance)
    
    
    min_toll_list = []        
    for i in IDs:
            for j in IDs:
                if i != j:
                    min_toll = get_min_toll(all_output, i, j)
                    min_toll_list.append(min_toll)
                    
  
                    
    min_time_list = []        
    for i in IDs:
            for j in IDs:
                if i != j:
                    min_time = get_min_time(all_output, i, j)
                    min_time_list.append(min_time)
                
    def saveDataOutput():
        df_dataMatrix = pd.DataFrame(all_output)
        # define csv file name
        data_file_name = 'data_' + now.strftime("%Y-%m-%d_%H-%M-%S") + '.csv'
            # set column names
        df_dataMatrix.columns = ['from', 'to', 'pathDistance', 'path', 'toll', 'time']
        
        if not os.path.exists(directory_outputs):
            os.makedirs(directory_outputs)
        # writing data frame to a CSV file
        df_dataMatrix.to_csv(os.path.join(directory_outputs, data_file_name), sep=';', index=False)
    
    saveDataOutput()