# Meta-documentation for MkDocs setup

## Environment setup

```bash
conda create -n mkdocs python=3.7
pip install mkdocs mkdocs_material
```

## Local development

```bash
cd /path/to/lfads-run-manager
conda activate mkdocs
mkdocs serve
```

## Publish to github

```bash
cd /path/to/lfads-run-manager
conda activate mkdocs
mkdocs gh-deploy
```
