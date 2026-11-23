#!/usr/bin/env python3
"""Run scAllele 0.0.9.7 with the pipeline's fixes, without modifying the installed package.

Usage: scallele_patched.py <scAllele arguments>
Run it with the python of the scAllele env (scallele_run in scallele_chrom.sh does this).

The installed scAllele package is copied once per working directory into
.scallele_patched/ and the patches in assets/scallele_patches/ (or $SCALLELE_PATCH_DIR)
are applied to the copy:
  1_write_vcf_SM_names      multi-BAM VCF sample columns
  2_get_sequence_left_edge  indel at the left edge of a read cluster no longer aborts the contig
  3_faidx_per_process       one FASTA index per worker (intermittent IndexNotFoundError)
  4_pool_watchdog           a killed pool worker stops the run (exit 137, SCALLELE_WORKER_KILLED)
                            instead of hanging forever
Each hunk must match the installed source exactly once; a hunk whose added lines are already
present is skipped, so an env that was patched by hand also works. Anything else is an error.
"""
import glob
import hashlib
import importlib.util
import os
import shutil
import sys
import tempfile


def die(msg):
    sys.stderr.write(f"scallele_patched.py: {msg}\n")
    sys.exit(2)


def read_hunks(patch_file):
    """Return [(target_file, old_text, new_text, added_blocks)] from a unified diff."""
    hunks, target, old, new, added = [], None, None, None, None

    def flush():
        if old is not None:
            blocks = ["".join(b) for b in added if b]
            hunks.append((target, "".join(old), "".join(new), blocks))

    with open(patch_file) as fh:
        for line in fh:
            if line.startswith("--- "):
                continue
            if line.startswith("+++ "):
                flush(); old = new = None
                path = line[4:].split("\t")[0].strip()
                target = path.split("/", 1)[1] if path.startswith(("a/", "b/")) else path
                continue
            if line.startswith("@@"):
                flush(); old, new, added = [], [], [[]]
                continue
            if old is None or line.startswith("\\"):
                continue
            body = line[1:] if line[:1] in " +-" else line
            if line.startswith("-"):
                old.append(body)
            elif line.startswith("+"):
                new.append(body); added[-1].append(body)
            else:
                old.append(body); new.append(body)
                if added[-1]:
                    added.append([])
    flush()
    return hunks


def apply_patches(pkg_root, patch_files):
    for pf in patch_files:
        for i, (target, old, new, added) in enumerate(read_hunks(pf), 1):
            fpath = os.path.join(pkg_root, target)
            with open(fpath) as fh:
                src = fh.read()
            if src.count(new) == 1 or (added and all(src.count(b) == 1 for b in added)):
                continue
            n = src.count(old)
            if n != 1:
                die(f"{os.path.basename(pf)} hunk {i}: expected text found {n} times in {target} "
                    "(not scAllele 0.0.9.7?)")
            with open(fpath, "w") as fh:
                fh.write(src.replace(old, new))


def main():
    spec = importlib.util.find_spec("scAllele")
    if spec is None or not spec.submodule_search_locations:
        die("scAllele is not importable with this python; run it with the scAllele env's python")
    installed = list(spec.submodule_search_locations)[0]

    here = os.path.dirname(os.path.realpath(__file__))
    patch_dir = os.environ.get("SCALLELE_PATCH_DIR",
                               os.path.join(here, "..", "assets", "scallele_patches"))
    patch_files = sorted(glob.glob(os.path.join(patch_dir, "*.patch")))
    if not patch_files:
        die(f"no patches found in {patch_dir}")

    key = hashlib.sha1(installed.encode())
    for pf in patch_files:
        with open(pf, "rb") as fh:
            key.update(fh.read())
    cache = os.path.join(os.getcwd(), ".scallele_patched", key.hexdigest()[:12])

    if not os.path.isdir(os.path.join(cache, "scAllele")):
        os.makedirs(os.path.dirname(cache), exist_ok=True)
        build = tempfile.mkdtemp(dir=os.path.dirname(cache))
        shutil.copytree(installed, os.path.join(build, "scAllele"),
                        ignore=shutil.ignore_patterns("__pycache__"))
        apply_patches(build, patch_files)
        try:
            os.rename(build, cache)
        except OSError:
            shutil.rmtree(build, ignore_errors=True)
            if not os.path.isdir(os.path.join(cache, "scAllele")):
                raise

    sys.path.insert(0, cache)
    for mod in [m for m in sys.modules if m == "scAllele" or m.startswith("scAllele.")]:
        del sys.modules[mod]
    from scAllele import scAllele as entry
    if os.path.dirname(os.path.realpath(entry.__file__)) != os.path.realpath(os.path.join(cache, "scAllele")):
        die("patched copy was not imported")
    sys.argv[0] = "scAllele"
    sys.exit(entry.main())


if __name__ == "__main__":
    main()
