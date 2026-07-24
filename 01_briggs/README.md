# Briggs contour-deformation solver

## Folder structure

- `src/` — Julia solver scripts (`briggs*.jl`, `contour_test_v3.1.jl`)
- `postprocessing/` — plotting scripts: `plots_json.py`, `plot_contour_deformation*.m` (MATLAB)
- `results/json/` — saved contour-iteration output (`contour_iteration*.json`)
- `results/videos/` — rendered contour animations (git-ignored, local only)
- `hpc/` — `slrum.sh` Slurm job script
- `vib_ribbon/` — vibrating-ribbon (signalling) variant, see its own README

**Note on paths:** each solver writes its output JSON to the directory it is run from (e.g. `briggsv5.jl` → `contour_iteration_v5.json`), and the postprocessing scripts read the JSON filename from their working directory as well. After a run, move the output JSON into `results/json/` and generated videos into `results/videos/`. The Slurm script `hpc/slrum.sh` still points to the old cluster path `.../Eigenvalue/Briggs` — update it to `.../Eigenvalue/01_briggs/src` after pulling this structure on the cluster.

## Recommended Testing Order

1. `briggsv2.jl` first.  
   This is the main adaptive-zeta version. It uses the adaptive `zeta_alpha` setup that gave the best visual results, but adds a branch-tracking guard before accepting an omega-contour step. The purpose is to keep the good contour movement while preventing fake pinches caused by the upper and lower branches accidentally selecting the same eigenvalue. If this version works well, it is the strongest candidate for the main run.

2. `briggsv3.jl` as the fixed-zeta baseline.  
   This version removes the adaptive `zeta_alpha` mechanism and instead uses one fixed common zeta value for both the omega and alpha contours. It is meant as a clean comparison case. If this version behaves better than `briggsv2.jl`, then the adaptive-zeta logic may be part of the problem. If it behaves worse, then adaptive zeta is likely helping the contour deformation.

3. `briggsv3.1.jl` if `briggsv3.jl` is too weak, too slow, or does not move the contours enough.  
   This is the more aggressive fixed-zeta version. It uses a larger common zeta and stronger contour movement settings to test whether the fixed-zeta approach can be improved by pushing the contours harder. This version is useful for checking whether a larger zeta helps the method reach the pinch region faster, or whether it causes instability, branch collapse, or bad contour motion.

4.  `briggs_agr.jl` only for deeper debugging.  
   This is the aggressive adaptive-zeta diagnostics version. It is not mainly meant as the clean production run, but as a debugging tool. It prints and checks more information about branch distances, contour distances, directional movement, and spacing of `alpha_L_u` and `alpha_L_l`. Use it when the other versions show suspicious behavior and you need to understand whether the issue comes from tracking, spacing, zeta choice, or the alpha contour moving too close to a branch.
   
5. briggsv3.2.jl This version uses a fixed common repulsion strength, zeta_common = 4e-4, for both the omega and alpha contour potentials, but replaces the earlier reseeded branch-continuation logic with memory-based branch tracking. Instead of re-identifying the upper and lower branches from dominant_eigvals(...) at each new trial L-contour, the code now tracks each new upper/lower eigenvalue from the previously accepted upper/lower branch point using contour_alpha_L_memory(L_trial, alpha_L_u, alpha_L_l). Each trial omega step is accepted only if the memory-tracked branches remain finite, separated, continuous, and do not exceed a maximum allowed tracking jump. The output also prints trackU, trackL, and failed so that branch-tracking stability can be monitored directly during the run. This version is meant to test whether the earlier collapse near the pinch was caused by upper/lower branch misidentification rather than by the fixed ? value itself.

6. version 3.3 uses a fixed common repulsion strength, zeta_common = 4e-4, for both the omega and alpha contour potentials. The branches are still constructed with the original continuation routine contour_alpha_L_conti(L_trial), but this version adds a new overlap guard for branch tracking. After each trial omega step, the code checks whether any discrete point on the upper branch is also selected by the lower branch, even if the matching points occur at different indices. If such an overlap is detected, the trial step is rejected, omega_dt is reduced by half, and the code retries with a smaller step. This is meant to prevent false branch collapse caused by accidentally assigning the same eigenvalue to both upper and lower branches. The output is saved to contour_iteration_commonzeta_4e-4_overlapguard.json.

7. `briggsv4.jl` continues from version 3.3 and keeps the fixed common potential strength, `zeta_common = 4e-4`, for both contours. The main improvement is a stricter acceptance check for the omega contour. Before a new omega-contour position is accepted, the code rebuilds the upper and lower spatial branches and checks whether they remain separated. The branch classification is based on the signed projection relative to the alpha contour `F`: eigenvalues on one side of `F` are treated as the upper branch, and eigenvalues on the other side are treated as the lower branch. If the two branches overlap or collapse onto the same eigenvalue, the step is rejected and the omega time step is reduced.

8. `briggsv4.1.jl` is a small improvement of `briggsv4.jl`. It prevents the omega contour from accepting steps that are too small to be useful. This avoids cases where the algorithm technically accepts a step, but the contour almost does not move. It should be used if `briggsv4.jl` stagnates or makes only very small progress.

9. `briggsv5.jl` stil working on it
