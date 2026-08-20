#!/usr/bin/env python3
"""Patch briggsv5_ribbon.jl in place: frame-to-frame continuity in the tracker.

Diagnosis (v5 run, 15 frames):
  - frame 14 (omega_i=-0.04979) and frame 15 (-0.05379) are BOTH smooth
    (max step 0.0037 / 0.0039 vs DA_MAX 0.012) and both made of genuine
    eigenvalues -- neither jumps.
  - the two branches involved never collide: minimum separation 0.036 over
    omega_r in [0.12,0.20] x omega_i in [-0.046,-0.060]. No branch point.
  - vertical continuation from frame 14's value down to frame 15's level:
        omega_r=0.12375 (seed) : 0.26367-0.04141i -> 0.26815-0.05078i,
                                 frame 15 has 0.26725-0.05057i   (agrees, 9e-4)
        omega_r=0.15050        : 0.30095+0.00000i -> 0.30329+0.00108i,
                                 frame 15 has 0.34548-0.01060i   (differs 0.044)
        omega_r=0.24000        : 0.40660+0.00568i -> 0.40725-0.00015i,
                                 frame 15 has 0.65700+0.07157i   (differs 0.260)
    max step on that vertical at omega_r=0.24 is 0.0008 -- utterly unambiguous.
  => the SEED is right, the ENDPOINT is wrong. Frame 15 drifts onto the
     neighbouring branch between the seed and omega_r ~ 0.15, gradually, with
     every individual step well inside the guard.

The information the tracker is not using: between frames every node moves by at
most |d omega_i| ~ 4e-3, and the branch moves by <= 1.1e-2, against a branch
separation of >= 3.6e-2 -- a factor of 3 margin. The previous frame's alpha at
the SAME node index is therefore a far better predictor than the along-contour
secant, everywhere on the contour (the horizontals shift by d omega_i, the riser
nodes by less, the arc not at all).
"""
import sys, pathlib

P = pathlib.Path("/sessions/magical-upbeat-hamilton/mnt/Eigenvalue/01_briggs/"
                 "vib_ribbon/src/briggsv5_ribbon.jl")
s = P.read_text()
n = 0


def sub(old, new, tag):
    global s, n
    c = s.count(old)
    if c != 1:
        sys.exit(f"PATCH '{tag}' matched {c} times, expected 1")
    s = s.replace(old, new)
    n += 1
    print("  ok ", tag)


# ---------------------------------------------------------------- header note
sub("""#   Plus, cheap and optional: a GRID-CONVERGENCE check (grid_check_every) that""",
    """#   5. FRAME-TO-FRAME CONTINUITY (added after the first v5 run, which still
#      hopped at iteration 14 even at N=300 -- refining the contour was not
#      enough).  Diagnosis, from the 15-frame v5 JSON:
#        - frames 14 and 15 are BOTH smooth (max step 0.0037 / 0.0039 against
#          DA_MAX = 0.012) and both consist of genuine eigenvalues.  Neither
#          jumps: the two curves drift apart gradually from the seed onwards.
#        - the two branches involved NEVER collide.  Minimum separation 0.036
#          over omega_r in [0.12,0.20] x omega_i in [-0.046,-0.060].  So there
#          is no branch point and the continuation is not ambiguous in
#          principle -- the tracker simply loses it.
#        - continuing alpha VERTICALLY from frame 14's value down to frame 15's
#          level settles which is right:
#             omega_r = 0.12375  0.26367-0.04141i -> 0.26815-0.05078i
#                                frame 15: 0.26725-0.05057i   agrees to 9e-4
#             omega_r = 0.15050  0.30095+0.00000i -> 0.30329+0.00108i
#                                frame 15: 0.34548-0.01060i   differs by 0.044
#             omega_r = 0.24000  0.40660+0.00568i -> 0.40725-0.00015i
#                                frame 15: 0.65700+0.07157i   differs by 0.260
#          with a max step of 0.0008 on that vertical.  The seed is right, the
#          endpoint is wrong.
#      The fix uses the one piece of information the tracker had and ignored:
#      each frame re-seeded from scratch and tracked along the contour with no
#      memory of the previous frame.  But between frames every node moves by at
#      most |d omega_i| ~ 4e-3 (horizontals shift by exactly that, riser nodes by
#      less, the arc not at all), and the branch moves by <= 1.1e-2 -- against a
#      branch separation of >= 3.6e-2.  So the previous frame's alpha at the SAME
#      node index is a much better predictor than the along-contour secant, with
#      a factor-3 margin.  continue_branch() now uses it when it exists and falls
#      back to the secant only for frame 1.  This also removes the v4.4-style
#      seed-flip risk, since the seed node is predicted the same way.
#      Reported per frame: framedev_u/l = max |alpha_new - alpha_prev| over all
#      nodes.  NOTE this is informational, NOT a check -- with prev as the
#      predictor it is small by construction.  The independent check is still
#      fixedw_u/l, which compares against frame 1 at nodes whose omega is pinned.
#
#   Plus, cheap and optional: a GRID-CONVERGENCE check (grid_check_every) that""",
    "header")

# ------------------------------------------------------ DFRAME_MAX + predictor
sub("""const DA_MAX = 0.012

function continue_branch(nodes, spectra, sides, seed, branch_side;
                         da_max = DA_MAX, sep_frac = 0.6, side_tol = 0.0)""",
    """const DA_MAX = 0.012
# v5 CHANGE 5: how far a node's alpha may move between two consecutive frames
# before it is worth flagging. Measured on the v5 run: the true motion for
# d omega_i = 4e-3 is <= 1.1e-2, while the nearest rival branch is >= 3.6e-2
# away, so 0.02 sits comfortably between the two.
const DFRAME_MAX = 0.02

function continue_branch(nodes, spectra, sides, seed, branch_side;
                         da_max = DA_MAX, sep_frac = 0.6, side_tol = 0.0,
                         prev = nothing)""",
    "dframe const")

sub("""        evs = spectra[j]; sd = sides[j]
        if isempty(evs) || !isfin(pred)""",
    """        # v5 CHANGE 5: the previous frame's value at THIS node beats the
        # along-contour secant. Every node moves by at most |d omega_i| between
        # frames, so prev[j] is within ~1.1e-2 of the truth, while the nearest
        # rival branch is >= 3.6e-2 away. Frame 1 has no prev and falls back to
        # the secant, which is what defines the branch for the whole run.
        if prev !== nothing && isfin(prev[j])
            pred = prev[j]
        end
        evs = spectra[j]; sd = sides[j]
        if isempty(evs) || !isfin(pred)""",
    "predictor")

# ------------------------------------------------------------ track_path: prev
sub("""function track_path(nodes0, spectra0, sides0, seed, branch_side;
                    max_refine = 4, node_cap = 4000, da_max = DA_MAX, sep_frac = 0.6)
    nodes   = copy(nodes0)
    spectra = Vector{Vector{ComplexF64}}(spectra0)
    sides   = Vector{Vector{Float64}}(sides0)
    is_orig = trues(length(nodes))
    n_added = 0
    vals, flag, nfb = continue_branch(nodes, spectra, sides, seed, branch_side;
                                      da_max = da_max, sep_frac = sep_frac)""",
    """function track_path(nodes0, spectra0, sides0, seed, branch_side;
                    max_refine = 4, node_cap = 4000, da_max = DA_MAX, sep_frac = 0.6,
                    prev = nothing)
    nodes   = copy(nodes0)
    spectra = Vector{Vector{ComplexF64}}(spectra0)
    sides   = Vector{Vector{Float64}}(sides0)
    pv      = prev === nothing ? nothing : ComplexF64.(collect(prev))
    is_orig = trues(length(nodes))
    n_added = 0
    vals, flag, nfb = continue_branch(nodes, spectra, sides, seed, branch_side;
                                      da_max = da_max, sep_frac = sep_frac, prev = pv)""",
    "track_path head")

sub("""        nn = ComplexF64[]; ss = Vector{Vector{ComplexF64}}()
        vv = Vector{Vector{Float64}}(); oo = Bool[]
        p = 1
        for j in 1:length(nodes)
            while p <= length(ins_at) && ins_at[p] == j
                push!(nn, ins_w[p]); push!(ss, new_spec[p]); push!(vv, new_sides[p]); push!(oo, false); p += 1
            end
            push!(nn, nodes[j]); push!(ss, spectra[j]); push!(vv, sides[j]); push!(oo, is_orig[j])
        end
        nodes = nn; spectra = ss; sides = vv; is_orig = oo
        vals, flag, nfb = continue_branch(nodes, spectra, sides, seed, branch_side;
                                          da_max = da_max, sep_frac = sep_frac)""",
    """        nn = ComplexF64[]; ss = Vector{Vector{ComplexF64}}()
        vv = Vector{Vector{Float64}}(); oo = Bool[]
        pp = pv === nothing ? nothing : ComplexF64[]
        p = 1
        for j in 1:length(nodes)
            while p <= length(ins_at) && ins_at[p] == j
                push!(nn, ins_w[p]); push!(ss, new_spec[p]); push!(vv, new_sides[p]); push!(oo, false)
                # inserted nodes sit midway, so predict midway too
                pv === nothing || push!(pp, 0.5 * (pv[j-1] + pv[j]))
                p += 1
            end
            push!(nn, nodes[j]); push!(ss, spectra[j]); push!(vv, sides[j]); push!(oo, is_orig[j])
            pv === nothing || push!(pp, pv[j])
        end
        nodes = nn; spectra = ss; sides = vv; is_orig = oo; pv = pp
        vals, flag, nfb = continue_branch(nodes, spectra, sides, seed, branch_side;
                                          da_max = da_max, sep_frac = sep_frac, prev = pv)""",
    "track_path refine")

# --------------------------------------------------------- track_outward: prev
sub("""function track_outward(nodes, spectra, sides, i_seed, seed, branch_side)
    vf, af, nf = track_path(nodes[i_seed:end],    spectra[i_seed:end],    sides[i_seed:end],    seed, branch_side)
    vb, ab, nb = track_path(nodes[i_seed:-1:1],   spectra[i_seed:-1:1],   sides[i_seed:-1:1],   seed, branch_side)""",
    """function track_outward(nodes, spectra, sides, i_seed, seed, branch_side; prev = nothing)
    pf = prev === nothing ? nothing : prev[i_seed:end]
    pb = prev === nothing ? nothing : prev[i_seed:-1:1]
    vf, af, nf = track_path(nodes[i_seed:end],    spectra[i_seed:end],    sides[i_seed:end],    seed, branch_side; prev = pf)
    vb, ab, nb = track_path(nodes[i_seed:-1:1],   spectra[i_seed:-1:1],   sides[i_seed:-1:1],   seed, branch_side; prev = pb)""",
    "track_outward")

# ------------------------------------------------------ track_on_contour: prev
sub("""function track_on_contour(Lc, i_seed)""",
    """function track_on_contour(Lc, i_seed; prev_u = nothing, prev_l = nothing)""",
    "toc head")

sub("""    au, add_u, nfb_u = track_outward(Lc, spectra, sides, i_seed, seed_u, :upper)
    al, add_l, nfb_l = track_outward(Lc, spectra, sides, i_seed, seed_l, :lower)
    return au, al, spectra, sides, (add_u, nfb_u, add_l, nfb_l)""",
    """    au, add_u, nfb_u = track_outward(Lc, spectra, sides, i_seed, seed_u, :upper; prev = prev_u)
    al, add_l, nfb_l = track_outward(Lc, spectra, sides, i_seed, seed_l, :lower; prev = prev_l)
    return au, al, spectra, sides, (add_u, nfb_u, add_l, nfb_l)""",
    "toc body")

sub("""function ribbon_branches(Lc)
    au, al, spectra, sides, diag = track_on_contour(Lc, i_seed_coarse)""",
    """# v5 CHANGE 5: the previous accepted frame's branches, used as the predictor.
# Empty on frame 1 (falls back to the along-contour secant); ignored on a
# refined contour, where the lengths do not match.
prev_alpha_u = ComplexF64[]
prev_alpha_l = ComplexF64[]
framedev_u = 0.0
framedev_l = 0.0

function ribbon_branches(Lc)
    pu = length(prev_alpha_u) == length(Lc) ? prev_alpha_u : nothing
    pl = length(prev_alpha_l) == length(Lc) ? prev_alpha_l : nothing
    au, al, spectra, sides, diag = track_on_contour(Lc, i_seed_coarse; prev_u = pu, prev_l = pl)""",
    "ribbon_branches")

# ----------------------------------------------- frame deviation + commit hook
sub("""fixedw_ref_u = ComplexF64[]""",
    """# How far the branches moved since the last accepted frame. INFORMATIONAL
# ONLY -- with prev as the predictor this is small by construction, so it is a
# report, not a test. The independent test is fixedw_u/l below, which compares
# against frame 1 at nodes whose omega never moves.
function frame_deviation!()
    if !isempty(prev_alpha_u) && length(prev_alpha_u) == length(alpha_L_u)
        global framedev_u = maximum(abs.(ComplexF64.(alpha_L_u) .- prev_alpha_u))
        global framedev_l = maximum(abs.(ComplexF64.(alpha_L_l) .- prev_alpha_l))
    else
        global framedev_u = 0.0
        global framedev_l = 0.0
    end
    return framedev_u, framedev_l
end
function commit_frame_branches!()
    global prev_alpha_u = ComplexF64.(alpha_L_u)
    global prev_alpha_l = ComplexF64.(alpha_L_l)
    return nothing
end

fixedw_ref_u = ComplexF64[]""",
    "frame_deviation")

# ------------------------------------------------------------ frame reporting
sub("""    @printf("        v5    : fixed-omega=%.3e/%.3e | gridconv=%.3e/%.3e | alpha(omega_f)=%+.6f%+.6fi / %+.6f%+.6fi | Im-crossings u=%d l=%d\\n",
            fixedw_u, fixedw_l, gridconv_u, gridconv_l,
            real(a_u_wf), imag(a_u_wf), real(a_l_wf), imag(a_l_wf), ncr_u, ncr_l)""",
    """    @printf("        v5    : fixed-omega=%.3e/%.3e | gridconv=%.3e/%.3e | framedev=%.3e/%.3e | alpha(omega_f)=%+.6f%+.6fi / %+.6f%+.6fi | Im-crossings u=%d l=%d\\n",
            fixedw_u, fixedw_l, gridconv_u, gridconv_l, framedev_u, framedev_l,
            real(a_u_wf), imag(a_u_wf), real(a_l_wf), imag(a_l_wf), ncr_u, ncr_l)""",
    "frame print")

sub("""        "gridconv_u"  => jsonnum(gridconv_u), "gridconv_l" => jsonnum(gridconv_l),""",
    """        "gridconv_u"  => jsonnum(gridconv_u), "gridconv_l" => jsonnum(gridconv_l),
        "framedev_u"  => jsonnum(framedev_u), "framedev_l" => jsonnum(framedev_l),""",
    "frame dict")

# ------------------------------------------------------- startup: commit frame
sub("""    fixed_omega_check!()          # v5: take the fixed-omega reference from frame 1
    grid_convergence_check!()     # v5: and one convergence check straight away""",
    """    frame_deviation!()            # v5: zero on frame 1 (there is no previous)
    fixed_omega_check!()          # v5: take the fixed-omega reference from frame 1
    grid_convergence_check!()     # v5: and one convergence check straight away""",
    "startup dev")

sub("""    last_frame = make_ribbon_frame(iteration_step)
    write_frame!(frame_writer, last_frame)
    iteration_step += 1
end""",
    """    last_frame = make_ribbon_frame(iteration_step)
    write_frame!(frame_writer, last_frame)
    commit_frame_branches!()      # v5: frame 1 becomes the predictor for frame 2
    iteration_step += 1
end""",
    "startup commit")

# ------------------------------------------------------ main loop: dev/commit
sub("""        load_on_workers()
        global last_frame = make_ribbon_frame(iteration_step)
        write_frame!(frame_writer, last_frame)
        global iteration_step += 1""",
    """        load_on_workers()
        frame_deviation!()        # v5: measure BEFORE prev is overwritten
        global last_frame = make_ribbon_frame(iteration_step)
        write_frame!(frame_writer, last_frame)
        commit_frame_branches!()  # v5: this frame becomes the next predictor
        global iteration_step += 1""",
    "loop commit")

# --------------------------------------------------------------- final report
sub("""    @printf("verification: fixed-omega %.3e / %.3e (tol 1e-8), grid convergence %.3e / %.3e\\n",
            fixedw_u, fixedw_l, gridconv_u, gridconv_l)""",
    """    @printf("verification: fixed-omega %.3e / %.3e (tol 1e-8, the independent check),\\n"
            * "              grid convergence %.3e / %.3e, frame deviation %.3e / %.3e (informational)\\n",
            fixedw_u, fixedw_l, gridconv_u, gridconv_l, framedev_u, framedev_l)""",
    "final report")

P.write_text(s)
print(f"\n{n} patches applied -> {P}")
