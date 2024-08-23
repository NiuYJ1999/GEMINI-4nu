# GEMINI++4$\nu$
* A modification version GEMINI++ code designed for de-excitation. 
* Most of codes are the original version from https://github.com/jdfrankland/gemini/tree/upstream
* Reference:
  * R. J. Charity, Physics Review C 82 (2010)014610  (evaporation parameters)
  * D. Mancusi, R. J. Charity, J. Cugnon, Physics Review C 82 (2010) 044610 (fission parameter)
  * R. J. Charity, inJoint ICTP-AIEA Advanced Workshop on Model Codes for Spallation Reactions (IAEA Vienna, 2008), Report INDC(NDC)-0530

# How to Use?
After cloneing this code to your workspace, just in the path of 'bashrc.sh' 
```
source bashrc.sh 
```
Then, you can use it easily.


This code has two modes, event mode and batch mode
* For the event mode (--event), after specify the charge and mass number for Z, A and the excited energy and spin of the parent nucleus, generate only one event
```
# mode 1: event mode
deexG --event {Z} {A} {excited energy} {spin} --suppress {$}
### example: for 11B* with 30 MeV excited energy and spin 0.5
deexG --event 5 11 30 0.5 --suppress 1
```
* For the batch mode (--batch), similar to the event mode, but the times of simulation and path to save the root file need to be specify
```
# mode 2: batch mode
deexG --batch {Z} {A} {excited energy} {spin} {Simulation times} {path} --suppress {$}
# example: simulate 1000 times, and put .root file in current path
deexG --event 5 11 30 0.5 1000 . --suppress 1
```
* Another need to be mentioned is the suppress factor, this need to be used for 3 cases
  * "1" means do not modify the emission probability of different particles
  * "0.5" means that apply a suppression factor 0.5 to charged particles for nucleus whose mass number less than 16 (Larger need to be studied further)
  * "free" menas that you can modify the factor on your own taste
  * Those results below 16 MeV need to be study further, gamma decay study is on-going
