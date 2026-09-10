====================
Visualise vibrations
====================

(molmatmodel) milias@DESKTOP-7OTLCGO:~/work/projects/papers/mmcp2026-codes/surfaces/gold/Hg_on_gold111/thermal_corrections/Hg_on_Au111_modular/.python visualize_hg_vibrations_cache_auto.py --structure HgAu_opt.traj  --phonon-dir HgAu_phonons

Hg vibration visualization
===========================
Structure       : HgAu_opt.traj
Phonon cache    : HgAu_phonons
Output directory: Hg_vibrations
Full atoms      : 97
Frozen atoms    : 64
Frozen indices  : [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63]
Mobile indices  : [64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96]
Mobile atoms    : 33
Hg mobile index : [32]
Positive modes  : 99
Positive range  : 8.749–138.777 cm^-1

Selected modes:
  mode   3:     19.641 cm^-1, Hg projection = 0.8693
  mode   2:     19.581 cm^-1, Hg projection = 0.8681
  mode   6:     29.868 cm^-1, Hg projection = 0.2719

Wrote: Hg_vibrations/Hg_mode_019.641_cm-1.traj
       Hg_vibrations/Hg_mode_019.641_cm-1.xyz
Wrote: Hg_vibrations/Hg_mode_019.581_cm-1.traj
       Hg_vibrations/Hg_mode_019.581_cm-1.xyz
Wrote: Hg_vibrations/Hg_mode_029.868_cm-1.traj
       Hg_vibrations/Hg_mode_029.868_cm-1.xyz

Equilibrium: Hg_vibrations/HgAu_equilibrium.traj
Diagnostics: Hg_vibrations/visualization_diagnostics.json

Open an animation with ASE GUI, for example:
  ase gui Hg_vibrations/Hg_mode_019.641_cm-1.traj

No new force calculations were performed.

whole slab frozen
-----------------

(molmatmodel) milias@DESKTOP-7OTLCGO:~/work/projects/papers/mmcp2026-codes/surfaces/gold/Hg_on_gold111/thermal_corrections/Hg_on_Au111_modular/.python visualize_hg_vibrations.py  --structure HgAu_opt.traj  --phonon-dir HgAu_phonons

Hg vibration visualization
===========================
Structure       : HgAu_opt.traj
Phonon cache    : HgAu_phonons
Output directory: Hg_vibrations
Full atoms      : 97
Frozen atoms    : 96
Frozen indices  : [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95]
Mobile indices  : [96]
Mobile atoms    : 1
Hg mobile index : [0]
Positive modes  : 3
Positive range  : 22.213–58.330 cm^-1

Selected modes:
  mode   2:     58.330 cm^-1, Hg projection = 1.0000
  mode   0:     22.213 cm^-1, Hg projection = 1.0000
  mode   1:     22.335 cm^-1, Hg projection = 1.0000

Wrote: Hg_vibrations/Hg_mode_058.330_cm-1.traj
       Hg_vibrations/Hg_mode_058.330_cm-1.xyz
Wrote: Hg_vibrations/Hg_mode_022.213_cm-1.traj
       Hg_vibrations/Hg_mode_022.213_cm-1.xyz
Wrote: Hg_vibrations/Hg_mode_022.335_cm-1.traj
       Hg_vibrations/Hg_mode_022.335_cm-1.xyz

Equilibrium: Hg_vibrations/HgAu_equilibrium.traj
Diagnostics: Hg_vibrations/visualization_diagnostics.json

Open an animation with ASE GUI, for example:
  ase gui Hg_vibrations/Hg_mode_058.330_cm-1.traj

No new force calculations were performed.
