# Wannier90 documentation

Wannier90 documentation using `mkdocs`.

## Project layout

```bash
mkdocs.yml    # The configuration file.
docs/
    index.md  # The documentation homepage.
    ...       # Other markdown pages, images and other files.
```

## Installation

```bash
pip install -r requirements.txt
```

## Commands

* `mkdocs new [dir-name]` - Create a new project.
* `mkdocs serve` - Start the live-reloading docs server.
* `mkdocs build` - Build the documentation site.
* `mkdocs -h` - Print help message and exit.

## Versioning

[`mike`](https://github.com/jimporter/mike) is used to manage the versioning of
the documentation.

In general, one can manually deploy the docs when a new version is released, by

```bash
mike deploy --push --update-aliases v[MAJOR].[MINOR] latest
```

Then the docs will be committed and pushed to the `gh-pages` branch.

For development, one can use

```bash
mike serve
```

To list the available versions, use

```bash
mike list
```

### `miki` initialization

For future reference, these are the commands used for initializing the
`gh-pages` branch

```bash
mike delete --all  # clean gh-pages branch
mike deploy --push --update-aliases v3.1.0 latest  # build docs
mike set-default latest  # set default redirect to latest
```

No need to run these commands again!
For future releases, just use `mike deploy`.

### References

* <https://squidfunk.github.io/mkdocs-material/setup/setting-up-versioning>
* <https://github.com/squidfunk/mkdocs-material-example-versioning>
* <https://github.com/jimporter/mike>

## Notes on conversion

The original wannier90 latex documentation was converted to markdown using
`pandoc`. The resulting markdown files were then manually edited to fit the
`mkdocs` format.

```bash
pandoc -s wannier90/doc/user_guide/user_guide.tex -o user_guide.md
```
