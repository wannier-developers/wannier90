# Wannier90 documentation

Wannier90 documentation using `mkdocs`.

## Project layout

    mkdocs.yml    # The configuration file.
    docs/
        index.md  # The documentation homepage.
        ...       # Other markdown pages, images and other files.

## Installation

```bash
pip install -r requirements.txt
```

## Commands

* `mkdocs new [dir-name]` - Create a new project.
* `mkdocs serve` - Start the live-reloading docs server.
* `mkdocs build` - Build the documentation site.
* `mkdocs -h` - Print help message and exit.

## Notes on conversion

The original wannier90 latex documentation was converted to markdown using `pandoc`. The resulting markdown files were then manually edited to fit the `mkdocs` format.

```bash
pandoc -s wannier90/doc/user_guide/user_guide.tex -o user_guide.md
```
