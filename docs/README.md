# w90-docs

Repo for prototyping wannier90 documentation using `mkdocs`.

The original w90 latex documentation was converted to markdown using pandoc. The resulting markdown files were then manually edited to fit the mkdocs format.

```bash
pandoc -s wannier90/doc/user_guide/user_guide.tex -o user_guide.md
```

Note pandoc does not support `verb#1#` style latex code, I need to manually replace these by `verb|1|` style code.
