import os
import argparse
import parser

parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('Z', type=int, help='an integer for the atomic number')
parser.add_argument('A', type=int, help='an integer for the mass number')
parser.add_argument('specifyPath', type=str, help='a string for the specify path')
parser.add_argument('Ex', type=float, help='Ex')
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
Ex = args.Ex
i = int(Ex)-1

if Z == 5 and A == 11:
    my_array = my_array_B11
elif Z == 7 and A == 15:
    my_array = my_array_N15
   
number = my_array[i] 
print(number)
specifyPath = args.specifyPath
path = "/junofs/users/niuyujie/GEMINI++4nu/ROOT/" + specifyPath
os.makedirs(path, exist_ok=True)
os.system("source /junofs/users/niuyujie/GEMINI++4nu/bashrc")
os.system(f"/junofs/users/niuyujie/GEMINI++4nu/gemini/./deexGen --batch {Z} {A} {Ex} {J} {number} {path}")

