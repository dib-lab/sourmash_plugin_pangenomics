# sourmash_plugin_pangenomics: tools for sourmash-based pangenome analyses

## Installation

```
pip install sourmash_plugin_pangenomics
```

## Usage

more here @CTB

## Support

We suggest filing issues in [the main sourmash issue tracker](https://github.com/dib-lab/sourmash/issues) as that receives more attention!

## Dev docs

`sourmash_plugin_pangenomics` is developed at https://github.com/sourmash-bio/sourmash_plugin_pangenomics.

### Testing

Run:
```
pytest tests
```

### Generating a release

Bump version number in `pyproject.toml` and push.

Make a new release on github.

Then pull, and:

```
python -m build
```

followed by `twine upload dist/...`.
