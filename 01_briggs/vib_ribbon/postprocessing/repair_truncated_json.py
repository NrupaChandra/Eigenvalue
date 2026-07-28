#!/usr/bin/env python3
"""
repair_truncated_json.py -- salvage a contour_iteration_*.json that was cut off
mid-write.

Up to v4.3 the Julia scripts re-serialised and rewrote the WHOLE frame array on
every iteration. If a write ever came up short (disk quota, killed job, network
filesystem hiccup), the file was left ending in the middle of an object and no
reader could open it -- MATLAB's jsondecode and Julia's JSON.parse both fail on
the whole file, so every frame appeared lost.

They are not lost. Because each write restarted from frame 1, a truncated file
still contains a valid prefix: frames 1..k of the last write. This script finds
the last complete frame, closes the array, and writes a loadable file.

Usage:
    python repair_truncated_json.py broken.json                 # -> broken_repaired.json
    python repair_truncated_json.py broken.json fixed.json
    python repair_truncated_json.py broken.json --inplace

(v4.3.1 onwards streams frames and keeps the closing "]" on disk at all times,
so this should not be needed again.)
"""

import json
import os
import sys


def find_last_complete_frame(text):
    """Scan the top-level array and return the end offset of the last complete
    element, tracking brace depth while ignoring braces inside strings."""
    i = text.find("[")
    if i < 0:
        raise ValueError("no opening '[' found -- this is not a frame array")

    depth = 0
    in_string = False
    escaped = False
    last_good = None
    n_frames = 0

    for pos in range(i + 1, len(text)):
        c = text[pos]

        if in_string:
            if escaped:
                escaped = False
            elif c == "\\":
                escaped = True
            elif c == '"':
                in_string = False
            continue

        if c == '"':
            in_string = True
        elif c == "{":
            depth += 1
        elif c == "}":
            depth -= 1
            if depth == 0:
                last_good = pos
                n_frames += 1

    if last_good is None:
        raise ValueError("no complete frame found -- the file is too damaged")

    return last_good, n_frames


def repair(src, dst):
    size = os.path.getsize(src)
    with open(src, "r", encoding="utf-8") as f:
        text = f.read()

    stripped = text.rstrip()
    if stripped.endswith("]"):
        try:
            data = json.loads(text)
            print(f"{src} is already valid JSON ({len(data)} frames, "
                  f"{size/1048576:.1f} MB) -- nothing to repair.")
            return 0
        except json.JSONDecodeError:
            print("File ends with ']' but does not parse; repairing anyway.")

    end, n_frames = find_last_complete_frame(text)
    repaired = text[:end + 1] + "]"

    data = json.loads(repaired)          # verify before writing anything
    iters = [d.get("iteration") for d in data]

    with open(dst, "w", encoding="utf-8") as f:
        f.write(repaired)

    print(f"Input : {src}  ({size/1048576:.1f} MB, truncated)")
    print(f"Output: {dst}  ({len(repaired)/1048576:.1f} MB)")
    print(f"Recovered {len(data)} complete frames "
          f"(iteration {iters[0]} .. {iters[-1]}).")
    print(f"Discarded {size - len(repaired):,} bytes of incomplete trailing data.")
    return len(data)


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    inplace = "--inplace" in sys.argv

    if not args:
        print(__doc__)
        sys.exit(1)

    src = args[0]
    if inplace:
        dst = src
    elif len(args) > 1:
        dst = args[1]
    else:
        root, ext = os.path.splitext(src)
        dst = root + "_repaired" + ext

    repair(src, dst)


if __name__ == "__main__":
    main()
