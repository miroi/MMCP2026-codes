===============
Hg on Gold(111)
===============

ChatGPT
-------
https://chatgpt.com/share/6a988d80-7f7c-83ed-9c43-68aa4b3f59ee

run
---
python run_hg_au111.py

working versions
----------------
e32594b0a409a2fa4f300efaa8c231b48707b87d
1d1676e4b95400781602dfc421132bb39a2bcefa  <--- error in Hads calc !!!
f300ea8bed672dda214041ffdf252a0806f7c9c8  <---- current calculation is numerically converged, but the thermodynamic correction is overestimated because the ZPE treatment is not atom-projected.

e8771ea458e3fb84863d9c6f5b061af2c55a0291  
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
positive Gads !
add full Hg gas entropy
make pressure dependence explicit
calculate ΔG(P,T) correctly

last 0error , provide fix of phonons.py

TODO: The only remaining check will be whether thermodynamics.py is using:

$$ \Delta ZPE = ZPE(HgAu)-ZPE(Au) $$

and not the Hg-projected modes for the adsorption enthalpy. That is the key correction for your previous unphysical entropy values.

d5f21e46c39c304a0c18f9eddc0f0ee839de99f2
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
T = 298.15 K
ΔE_ads  -0.501323 eV
ΔZPE    0.008909 eV
ΔH_ads  -0.423675 eV
ΔS_ads  -2.154333e-04 eV/K
ΔG_ads  -0.359443 eV

9271f711c53afcedbad2b9bb013a86d2fdb9ea01
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
T = 298.15 K
p = 1 bar
ΔE_ads          = -0.50132288 eV
ΔZPE            = 0.00890894 eV
ΔU_vib          = 0.07764820 eV
H_gas thermal   = 0.06423145 eV
ΔH_ads          = -0.48790613 eV
S_gas           = 1.81344702e-03 eV/K
ΔS_vib          = 7.23468926e-04 eV/K
ΔS_ads          = -1.08997809e-03 eV/K
ΔF_vib          = -0.13805406 eV
μ_gas thermal   = -0.47644778 eV
ΔG_ads          = -0.16292917 eV

TO CONTINUE:
A better treatment would be to use a gas reference including Hg translational/rotational entropy, not just Hg vibrational modes.
I would only replace:
thermodynamics.py (final consistent convention)
small workflow.py result handling

I will generate the two complete files as downloadable .py replacements. I need one more thing first: please upload your current thermodynamics.py and workflow.py.

Chat paused until usage resets at 11:58 PM
You’ve reached the limit for chats that include files or images. Start a new text-only chat or upgrade to continue now.


Continuing link
---------------
https://chatgpt.com/s/t_6a99e662377c8191ac1237422af7398b

d1661f75aa213142e94d646bd465fe92a067d287  <--- working for single

https://chatgpt.com/s/t_6a99f7933f288191bd95a141e39821d1  ... last link, to fix workflow.py


TEMPERATURE WHERE ΔG_ads = 0
No ΔG_ads = 0 crossing found.
Search interval     : [200.00, 400.00] K
ΔG_ads at limits    : -0.27065327, -0.05256659 eV

a68ae6eb1f8a26ab77fc27b141f997d542770711 <--- working for scan,  ΔG_ads = 0


https://chatgpt.com/share/6a9a4ad1-5470-83eb-97ce-a1a1783439dd  <--- best analysis of working code, for iliasmiro01@gmail.com, waiting for Hg_on_Au111_modular_FIXED.zip to be downloaded

for p = 1e-09 bar
TEMPERATURE WHERE ΔG_ads = 0
No ΔG_ads = 0 crossing found.
Search interval     : [200.00, 600.00] K
ΔG_ads at limits    : 0.08553233, 1.22972515 eV

 p = 1bar  Tdep = 450K
 p = 1e-5bar Tdep = 235 K
 p=1e-9 bar Tdep=170.4K

commit a9980313bf6ecb3d3aa729ad71c9eda1fb308e9c    <--- working after some update


0f7721d7b1df64964cd625db931be54805c2ff4f <--- working with ASEvibrations

https://chatgpt.com/s/t_6a9f35c90e4081919e9b45381d872c6b

commit a4630fe5fe4f58d5119f185d21f8b65d4a016cb3  <---  working ASEvibr + ASEthermo

control number of frozen layers
-------------------------------
new change - control number of frozen layers for phonons, see https://chatgpt.com/s/t_6aa00dd2de8c8191a85f4adfba730c2a  
works 04bb9e87b5cd42d946537d583efdd440b2ab78ff

restarting
----------
https://chatgpt.com/s/t_6aa036ec51408191a1e5b09070870abd

source files in https://disk.yandex.ru/d/cKgMbcBqufsccQ


Theory
------
theory - see https://chatgpt.com/s/t_6aa06f27c2bc8191b1a392cd00f359b8
Download ABOUT_THE_CODE_Supplementary_Information_EXPANDED.rst


last change iliasmiro07 , https://chatgpt.com/s/t_6aa15743085881919ff84296b7e63b7e
 TODO: get theory   

