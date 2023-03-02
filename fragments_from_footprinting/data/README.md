Solvent accessible surface area data. 

Copied from Devany's data located. Need to get original source
/ru-auth/local/home/dwest/scratch/ritracks/2021_12_histone_protect_tests/1kx5_sasa

## Including package data

Modify your package's `pyproject.toml` file.
Update the [tool.setuptools.package_data](https://setuptools.pypa.io/en/latest/userguide/datafiles.html#package-data)
and point it at the correct files.
Paths are relative to `package_dir`.

Package data can be accessed at run time with `importlib.resources` or the `importlib_resources` back port.
See https://setuptools.pypa.io/en/latest/userguide/datafiles.html#accessing-data-files-at-runtime
for suggestions.

If modules within your package will access internal data files using
[the recommended approach](https://setuptools.pypa.io/en/latest/userguide/datafiles.html#accessing-data-files-at-runtime),
you may need to include `importlib_resources` in your package dependencies.
In `pyproject.toml`, include the following in your `[project]` table.
```
dependencies = [
    "importlib-resources;python_version<'3.10'",
]
```

## Manifest

* `look_and_say.dat`: first entries of the "Look and Say" integer series, sequence [A005150](https://oeis.org/A005150)
