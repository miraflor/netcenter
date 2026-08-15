# Replacing an existing GitHub repo with this codebase

This ZIP is laid out as a **repository root**. There is no extra wrapper folder
inside it.

## Existing local repository

1. Open the existing repo folder in File Explorer.
2. Keep the hidden `.git` folder exactly where it is.
3. Delete or replace the old project files **except `.git`**.
4. Extract the replacement ZIP directly into the repo folder.
5. Open GitHub Desktop and review the changed files.
6. Open a terminal in the repo folder and run:

```bash
python -m pip install -e ".[dev]"
python -m pytest -q
python -m netcenter.cli --help
```

A clean v0.3.0 tree should report **69 passed**.

7. Commit only after the tests pass.

Suggested commit message:

```text
Release netcenter v0.3.0: strengthen topology and memory safety
```

Suggested description:

```text
Generalizes safe shared-vertex noding to true T- and X-junctions, removes
orphan nodes after geometry filtering, streams parallel Dijkstra blocks into a
preallocated matrix, restores conservative worker defaults, and aligns package
metadata, documentation, and regression tests.
```

## New repository

Create an empty GitHub repository, clone it, and place these files directly in
the clone root. Do not add another wrapper directory around them.

## Do not commit

The supplied `.gitignore` excludes Python caches, virtual environments, build
artifacts, and common local GIS outputs. Keep source datasets out of the repo
unless they are deliberately small, redistributable fixtures with clear
licensing.
