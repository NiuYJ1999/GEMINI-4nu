# GEMINI-4nu
A modification version GEMINI++ code designed for de-excitation

# How to Use?
After cloneing this code to your workspace, just
```
# in the 'bashrc.sh' path
source bashrc.sh 
```
Then, you can use it easily.


This code has two modes, event mode and batch mode
* For the event mode (--event), after specify the charge and mass number for Z, A and the excited energy and spin of the parent nucleus, generate only one event
```
deexG --event/batch {Z} {A} {excited energy} {spin} --suppress {$}
# example:
deexG --event 5 11 30 0.5 --suppress 1 # for 11B* with 30MeV excited energy and spin==0.5
```
* For the batch mode (--batch), similar to the event mode, but the times of simulation and path to save the root file need to be specify
```
deexG --event/batch {Z} {A} {excited energy} {spin} {Simulation times} {path}--suppress {$}
# example:
deexG --event 5 11 30 0.5 1000 . --suppress 1 # simulate 1000 times, and put .root file in current path
```
* Another need to be mentioned is the suppress factor, this need to be used for 3 cases
  * "1" means do not modify the emission probability of different particles
  * "0.5" means that apply a suppression factor 0.5 to charged particles for nucleus whose mass number less than 16 (Larger need to be studied further)
  * "free" menas that you can modify the factor on your own taste
