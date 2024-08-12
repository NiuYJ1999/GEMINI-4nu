import os
import glob
import argparse
import re
import ROOT

parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('Z', type=int, help='an integer for the atomic number')
parser.add_argument('A', type=int, help='an integer for the mass number')
parser.add_argument('specifyPath', type=str, help='a string for the specify path')
parser.add_argument('suppress', type=str, help='a string for the specify suppress')
args = parser.parse_args()

my_array = []
#B11
my_array_B11 = [0, 0, 0, 0, 0, 0, 0, 665, 7610, 6145, 3330, 5700, 6210, 9140, 14170, 18630, 20920, 19330, 17865, 18630, 18630, 18055, 17675, 16145, 14360, 13340, 12450, 12195, 11305, 10350, 9710, 8120, 7225, 6780, 6210, 5700, 5505, 4935, 4235, 3850, 3340, 3085, 2705, 2260, 1750, 1110,  660, 570, 490, 240]
#n15
my_array_N15 = [0, 0, 0, 0, 0, 0, 0, 1636, 5981, 19526, 22460, 1975, 4006, 2821, 2821, 2031, 5474, 5022, 4966, 4345, 4006, 4796, 5530, 5812, 6997, 7562, 7957, 8295, 7844, 8408, 7054, 7562, 7449, 6997, 6884, 5812, 5191, 4796, 4176, 3837, 3386, 3216, 2708, 2878, 2483, 2200, 1975, 1918, 1975, 1354]
#print(sum(my_array[15:35]))

Z = args.Z
A = args.A
J = 0.5
Ex = 0

if Z == 5 and A == 11:
    my_array = my_array_B11
elif Z == 7 and A == 15:
    my_array = my_array_N15
    
specifyPath = args.specifyPath
suppress = args.suppress
GROOT = os.getenv('GROOT')
GINPUT = os.getenv('GINPUT')
pathF = GROOT + specifyPath
os.makedirs(pathF, exist_ok=True)
#os.system("source $GTOP/bashrc.sh")
for i in my_array:
    Ex = Ex + 1
    os.system(f"deexG --batch {Z} {A} {Ex} {J} {i} {pathF} --suppress {suppress}")

while True:
    error_occurred = False
    files = glob.glob(os.path.join(pathF, f"*.root"))
    for file_path in files:
        try:
            file = ROOT.TFile(file_path, "READ")
            if file.IsZombie():
                print(f"This file {file_path} is a Zombie! Try again!")
                match = re.search(r'Ex_(\d+\.\d+)', file_path)
                if match:
                    number = int(float(match.group(1)))
                    times = my_array[number]
                    os.system(f"deexG --batch {Z} {A} {number} {J} {times} {pathF} --suppress {suppress}")
                    error_occurred = True
            else:
                keys = file.GetListOfKeys()
                if len(keys) == 0:
                    print("File has no keys.")
                    match = re.search(r'Ex_(\d+\.\d+)', file_path)
                    if match:
                        number = int(float(match.group(1)))
                        times = my_array[number]
                        os.system(f"deexG --batch {Z} {A} {number} {J} {times} {pathF} --suppress {suppress}")
                        error_occurred = True
                else:
                    # for key in keys:
                    #     print(f"Key name: {key.GetName()}, Key class: {key.GetClassName()}")
                    continue  
                file.Close()
        except Exception as e:
            print(f"Failed to open file {file_path}: {e}")
            match = re.search(r'Ex_(\d+\.\d+)', file_path)
            if match:
                number = int(float(match.group(1)))
                times = my_array[number]
                os.system(f"deexG --batch {Z} {A} {number} {J} {times} {pathF} --suppress {suppress}")
                error_occurred = True
    if not error_occurred:
        break
